/*
 * Copyright 2009 Colin Percival, 2011 ArtForz, 2011-2014 pooler
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 * This file was originally written by Colin Percival as part of the Tarsnap
 * online backup system.
 */

#include "miner.h"

#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <immintrin.h>

static const uint32_t keypad[12] = {
	0x80000000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0x00000280
};
static const uint32_t innerpad[11] = {
	0x80000000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0x000004a0
};
static const uint32_t outerpad[8] = {
	0x80000000, 0, 0, 0, 0, 0, 0, 0x00000300
};
static const uint32_t finalblk[16] = {
	0x00000001, 0x80000000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0x00000620
};

#define N (1024 * 1024)

struct ScryptData
{
    uint32_t X[32];
    uint32_t data[20 + 12];
    uint32_t hash[8];
    uint32_t midstate[8];
    uint32_t tstate[8];
    uint32_t ostate[8];
    unsigned char *scratchbuf;
    uint32_t* V;
    uint32_t* B;
    uint32_t *Bx;
};

struct ScryptDataSet
{
    struct ScryptData sd[4];
};

static inline void HMAC_SHA256_80_init(const uint32_t *key,
	uint32_t *tstate, uint32_t *ostate)
{
	uint32_t ihash[8];
	uint32_t pad[16];
	int i;

	/* tstate is assumed to contain the midstate of key */
	memcpy(pad, key + 16, 16);
	memcpy(pad + 4, keypad, 48);
	sha256_transform(tstate, pad, 0);
	memcpy(ihash, tstate, 32);

	sha256_init(ostate);
	for (i = 0; i < 8; i++)
		pad[i] = ihash[i] ^ 0x5c5c5c5c;
	for (; i < 16; i++)
		pad[i] = 0x5c5c5c5c;
	sha256_transform(ostate, pad, 0);

	sha256_init(tstate);
	for (i = 0; i < 8; i++)
		pad[i] = ihash[i] ^ 0x36363636;
	for (; i < 16; i++)
		pad[i] = 0x36363636;
	sha256_transform(tstate, pad, 0);
}

static inline void PBKDF2_SHA256_80_128(const uint32_t *tstate,
	const uint32_t *ostate, const uint32_t *salt, uint32_t *output)
{
	uint32_t istate[8], ostate2[8];
	uint32_t ibuf[16], obuf[16];
	int i, j;

	memcpy(istate, tstate, 32);
	sha256_transform(istate, salt, 0);
	
	memcpy(ibuf, salt + 16, 16);
	memcpy(ibuf + 5, innerpad, 44);
	memcpy(obuf + 8, outerpad, 32);

	for (i = 0; i < 4; i++) {
		memcpy(obuf, istate, 32);
		ibuf[4] = i + 1;
		sha256_transform(obuf, ibuf, 0);

		memcpy(ostate2, ostate, 32);
		sha256_transform(ostate2, obuf, 0);
		for (j = 0; j < 8; j++)
			output[8 * i + j] = swab32(ostate2[j]);
	}
}

static inline void PBKDF2_SHA256_128_32(uint32_t *tstate, uint32_t *ostate,
	const uint32_t *salt, uint32_t *output)
{
	uint32_t buf[16];
	int i;
	
	sha256_transform(tstate, salt, 1);
	sha256_transform(tstate, salt + 16, 1);
	sha256_transform(tstate, finalblk, 0);
	memcpy(buf, tstate, 32);
	memcpy(buf + 8, outerpad, 32);

	sha256_transform(ostate, buf, 0);
	for (i = 0; i < 8; i++)
		output[i] = swab32(ostate[i]);
}

static inline void salsa_round(__m128i* xv, const int dst, const int src1, const int src2, const int r)
{
    register __m128i _calc = _mm_add_epi32(xv[src1], xv[src2]);
    xv[dst] = _mm_xor_si128(_mm_xor_si128(xv[dst], _mm_slli_epi32(_calc, r)), _mm_srli_epi32(_calc, (32 - r)));
}
/*
*/
static inline void xor_salsa8_SSE(struct ScryptDataSet* ws)
{
    __m128i bv[16];
    __m128i bxv[16];
    __m128i xv[16];
    uint32_t* b = (uint32_t*)bv;
    uint32_t* bx = (uint32_t*)bxv;
    uint32_t* x = (uint32_t*)xv;
    for (int v = 0; v < 16; ++v)
    {
        b[v * 4 + 0] = ws->sd[0].B[v];
        b[v * 4 + 1] = ws->sd[1].B[v];
        b[v * 4 + 2] = ws->sd[2].B[v];
        b[v * 4 + 3] = ws->sd[3].B[v];
        bx[v * 4 + 0] = ws->sd[0].Bx[v];
        bx[v * 4 + 1] = ws->sd[1].Bx[v];
        bx[v * 4 + 2] = ws->sd[2].Bx[v];
        bx[v * 4 + 3] = ws->sd[3].Bx[v];
    }
    for (int v = 0; v < 16; ++v)
    {
        xv[v] = bv[v] = _mm_xor_si128(bv[v], bxv[v]);
    }
    for (int i = 0; i < 8; i += 2) {
        salsa_round(xv, 4, 0, 12, 7);
        salsa_round(xv, 9, 5, 1, 7);
        salsa_round(xv, 14, 10, 6, 7);
        salsa_round(xv, 3, 15, 11, 7);
        salsa_round(xv, 8, 4, 0, 9);
        salsa_round(xv, 13, 9, 5, 9);
        salsa_round(xv, 2, 14, 10, 9);
        salsa_round(xv, 7, 3, 15, 9);
        salsa_round(xv, 12, 8, 4, 13);
        salsa_round(xv, 1, 13, 9, 13);
        salsa_round(xv, 6, 2, 14, 13);
        salsa_round(xv, 11, 7, 3, 13);
        salsa_round(xv, 0, 12, 8, 18);
        salsa_round(xv, 5, 1, 13, 18);
        salsa_round(xv, 10, 6, 2, 18);
        salsa_round(xv, 15, 11, 7, 18);
        salsa_round(xv, 1, 0, 3, 7);
        salsa_round(xv, 6, 5, 4, 7);
        salsa_round(xv, 11, 10, 9, 7);
        salsa_round(xv, 12, 15, 14, 7);
        salsa_round(xv, 2, 1, 0, 9);
        salsa_round(xv, 7, 6, 5, 9);
        salsa_round(xv, 8, 11, 10, 9);
        salsa_round(xv, 13, 12, 15, 9);
        salsa_round(xv, 3, 2, 1, 13);
        salsa_round(xv, 4, 7, 6, 13);
        salsa_round(xv, 9, 8, 11, 13);
        salsa_round(xv, 14, 13, 12, 13);
        salsa_round(xv, 0, 3, 2, 18);
        salsa_round(xv, 5, 4, 7, 18);
        salsa_round(xv, 10, 9, 8, 18);
        salsa_round(xv, 15, 14, 13, 18);
    }
    for (int v = 0; v < 16; ++v)
    {
        bv[v] = _mm_add_epi32(bv[v], xv[v]);
    }
    for (int v = 0; v < 16; ++v)
    {
        ws->sd[0].B[v] = b[v * 4 + 0];
        ws->sd[1].B[v] = b[v * 4 + 1];
        ws->sd[2].B[v] = b[v * 4 + 2];
        ws->sd[3].B[v] = b[v * 4 + 3];
    }
}

static inline void scrypt_core_4(struct ScryptDataSet* ws)
{
    for (uint32_t i = 0; i < N; i++) {
        for (size_t idx = 0; idx < 4; ++idx)
        {
            memcpy(&ws->sd[idx].V[i * 32], ws->sd[idx].X, 128);
            ws->sd[idx].B = ws->sd[idx].X; ws->sd[idx].Bx = ws->sd[idx].X + 16;
        }
        xor_salsa8_SSE(ws);
        for (size_t idx = 0; idx < 4; ++idx)
        {
            ws->sd[idx].B = ws->sd[idx].X + 16; ws->sd[idx].Bx = ws->sd[idx].X;
        }
        xor_salsa8_SSE(ws);
    }
    for (uint32_t i = 0; i < N; i++) {
        for (size_t idx = 0; idx < 4; ++idx)
        {
            uint32_t j = 32 * (ws->sd[idx].X[16] & (N - 1));
            for (uint8_t k = 0; k < 32; k++)
            {
                ws->sd[idx].X[k] ^= ws->sd[idx].V[j + k];
            }
            ws->sd[idx].B = ws->sd[idx].X; ws->sd[idx].Bx = ws->sd[idx].X + 16;
        }
        xor_salsa8_SSE(ws);
        for (size_t idx = 0; idx < 4; ++idx)
        {
            ws->sd[idx].B = ws->sd[idx].X + 16; ws->sd[idx].Bx = ws->sd[idx].X;
        }
        xor_salsa8_SSE(ws);
    }
}

unsigned char *scrypt_buffer_alloc_sse(int nn)
{
    uchar* unalignedPtr = (uchar*)(malloc((size_t)N * 128 * 4 + 64));
    uchar* ptr = (uchar*)(((uintptr_t)unalignedPtr + 64) & ~(64 - 1));
	return ptr;
}

static void scrypt_1024_1_1_256(struct ScryptDataSet* scryptDataSet)
{
    for (size_t i = 0; i < 4; ++i)
    {
        memcpy(scryptDataSet->sd[i].tstate, scryptDataSet->sd[i].midstate, 32);
        HMAC_SHA256_80_init(scryptDataSet->sd[i].data, scryptDataSet->sd[i].tstate, scryptDataSet->sd[i].ostate);
        PBKDF2_SHA256_80_128(scryptDataSet->sd[i].tstate, scryptDataSet->sd[i].ostate, scryptDataSet->sd[i].data, scryptDataSet->sd[i].X);
    }

    scrypt_core_4(scryptDataSet);
    for (size_t i = 0; i < 4; ++i)
    {
        PBKDF2_SHA256_128_32(scryptDataSet->sd[i].tstate, scryptDataSet->sd[i].ostate, scryptDataSet->sd[i].X, scryptDataSet->sd[i].hash);
    }
}

int scanhash_scrypt_sse(int thr_id, struct work *work, uint32_t max_nonce, uint64_t *hashes_done, unsigned char *scratchbuf, uint32_t nn)
{
	uint32_t *pdata = work->data;
	uint32_t *ptarget = work->target;
	uint32_t n = pdata[19] - 1;
	const uint32_t Htarg = ptarget[7];
    struct ScryptDataSet dataSet;
    for (size_t i = 0; i < 4; ++i)
    {
        dataSet.sd[i].scratchbuf = scratchbuf + i * (size_t)N * 128;
        dataSet.sd[i].V = (uint32_t *)dataSet.sd[i].scratchbuf;
        memcpy(dataSet.sd[i].data, pdata, 80);

        sha256_init(dataSet.sd[i].midstate);
        sha256_transform(dataSet.sd[i].midstate, dataSet.sd[i].data, 0);
    }
	
	do {
        for (size_t i = 0; i < 4; ++i)
        {
            dataSet.sd[i].data[19] = ++n;
        }
		
		scrypt_1024_1_1_256(&dataSet);
		
        for (size_t i = 0; i < 4; ++i)
        {
            if (unlikely(dataSet.sd[i].hash[7] <= Htarg && fulltest(dataSet.sd[i].hash, ptarget))) {
                work_set_target_ratio(work, dataSet.sd[i].hash);
                *hashes_done = n - pdata[19] + 1;
                pdata[19] = dataSet.sd[i].data[19];
                printf("Found, idx %u\n", (unsigned int)i);
                return 1;
            }
        }
	} while (likely(n < max_nonce && !work_restart[thr_id].restart));
	
	*hashes_done = n - pdata[19] + 1;
	pdata[19] = n;
	return 0;
}

// simple cpu test (util.c) 
/*
void scrypthash(void *output, const void *input, uint32_t N)
{
	uint32_t midstate[8];
	char *scratchbuf = scrypt_buffer_alloc(N);

	memset(output, 0, 32);
	if (!scratchbuf)
		return;

	sha256_init(midstate);
	sha256_transform(midstate, input, 0);

	scrypt_1024_1_1_256((uint32_t*)input, (uint32_t*)output, midstate, scratchbuf, N);

	free(scratchbuf);
}
*/