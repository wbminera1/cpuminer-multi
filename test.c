#include <stdlib.h>  
#include <stdio.h>
#include <stdint.h>
#include <conio.h>
#include <memory.h>
#include <immintrin.h>
#include <inttypes.h>

#include "miner.h"

#define R(a, b) (((a) << (b)) | ((a) >> (32 - (b))))

struct WorkData
{
    uint32_t B[16]; 
    uint32_t Bx[16];
};

struct WorkSet
{
    struct WorkData WD[4];
};

void xor_salsa8(uint32_t B[16], const uint32_t Bx[16]);
void salsa20_wordtobyte_tr(uint8_t output[64], const uint32_t input[16]);

inline void swap(uint32_t* a, uint32_t* b)
{
    uint32_t tmp = *a;
    *a = *b;
    *b = tmp;
}

inline void transpose(uint32_t* m)
{
/*
    sX[0] = _mm_shuffle_epi32(row2, 0x93);
    row3 = _mm_shuffle_epi32(row3, 0x4e);
    row4 = _mm_shuffle_epi32(row4, 0x39);
*/
    swap(m+ 1, m+ 4);
    swap(m+ 2, m+ 8);
    swap(m+ 3, m+12);
    swap(m+ 6, m+ 9);
    swap(m+ 7, m+13);
    swap(m+11, m+14);
}

void test_transpose()
{
    uint32_t m[16] = { 0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0A, 0x0B, 0x0C, 0x0D, 0x0E, 0x0F };
    transpose(m);
    transpose(m);
}

int test_rotation()
{
    const uint32_t V1 = 0x6f54c89d;
    uint32_t r1_7 = R(V1, 7);
    uint32_t r1_9 = R(V1, 9);
    uint32_t r1_13 = R(V1, 13);
    uint32_t r1_18 = R(V1, 18);
    uint32_t r2_7 = _rotl(V1, 7);
    uint32_t r2_9 = _rotl(V1, 9);
    uint32_t r2_13 = _rotl(V1, 13);
    uint32_t r2_18 = _rotl(V1, 18);
    return (r1_7 == r2_7) && (r1_9 == r2_9) && (r1_13 == r2_13) && (r1_18 == r2_18);
}

/*static inline */void xor_salsa8_t(uint32_t B[16], const uint32_t Bx[16])
{
    //uint32_t x00, x01, x02, x03, x04, x05, x06, x07, x08, x09, x10, x11, x12, x13, x14, x15;
    uint32_t x[16];
    __m128i* sX = (__m128i*)x;
    __m128i* sB = (__m128i*)B;
    __m128i* sBx = (__m128i*)Bx;
    sB[0] = _mm_xor_si128(sB[0], sBx[0]);
    sB[1] = _mm_xor_si128(sB[1], sBx[1]);
    sB[2] = _mm_xor_si128(sB[2], sBx[2]);
    sB[3] = _mm_xor_si128(sB[3], sBx[3]);
    memcpy(sX, sB, sizeof(uint32_t) * 16);

    /*
    x[ 0] = (B[0] ^= Bx[0]);
    x[ 1] = (B[1] ^= Bx[1]);
    x[ 2] = (B[2] ^= Bx[2]);
    x[ 3] = (B[3] ^= Bx[3]);
    x[ 4] = (B[4] ^= Bx[4]);
    x[ 5] = (B[5] ^= Bx[5]);
    x[ 6] = (B[6] ^= Bx[6]);
    x[ 7] = (B[7] ^= Bx[7]);
    x[ 8] = (B[8] ^= Bx[8]);
    x[ 9] = (B[9] ^= Bx[9]);
    x[10] = (B[10] ^= Bx[10]);
    x[11] = (B[11] ^= Bx[11]);
    x[12] = (B[12] ^= Bx[12]);
    x[13] = (B[13] ^= Bx[13]);
    x[14] = (B[14] ^= Bx[14]);
    x[15] = (B[15] ^= Bx[15]);
*/

    for (int i = 0; i < 8; i += 2) {
        // Operate on columns.
        transpose(x);

        x[1] ^= _rotl(x[0] + x[3], 7);
        x[2] ^= _rotl(x[0] + x[1], 9);
        x[3] ^= _rotl(x[2] + x[1], 13);
        x[0] ^= _rotl(x[2] + x[3], 18);

        x[6] ^= _rotl(x[5] + x[4], 7);
        x[7] ^= _rotl(x[5] + x[6], 9);
        x[4] ^= _rotl(x[7] + x[6], 13);
        x[5] ^= _rotl(x[7] + x[4], 18);

        x[11] ^= _rotl(x[10] + x[9], 7);
        x[8] ^= _rotl(x[10] + x[11], 9);
        x[9] ^= _rotl(x[8] + x[11], 13);
        x[10] ^= _rotl(x[8] + x[9], 18);

        x[12] ^= _rotl(x[15] + x[14], 7);
        x[13] ^= _rotl(x[15] + x[12], 9);
        x[14] ^= _rotl(x[13] + x[12], 13);
        x[15] ^= _rotl(x[13]+ x[14], 18);

        transpose(x);

        /* Operate on rows. */
        x[1] ^= _rotl(x[0] + x[3], 7);
        x[2] ^= _rotl(x[1] + x[0], 9);
        x[3] ^= _rotl(x[2] + x[1], 13);
        x[0] ^= _rotl(x[3] + x[2], 18);

        x[6] ^= _rotl(x[5] + x[4], 7);
        x[7] ^= _rotl(x[6] + x[5], 9);
        x[4] ^= _rotl(x[7] + x[6], 13);
        x[5] ^= _rotl(x[4] + x[7], 18);

        x[11] ^= _rotl(x[10] + x[9], 7);
        x[8] ^= _rotl(x[11] + x[10], 9);
        x[9] ^= _rotl(x[8] + x[11], 13);
        x[10] ^= _rotl(x[9] + x[8], 18);

        x[12] ^= _rotl(x[15] + x[14], 7);
        x[13] ^= _rotl(x[12] + x[15], 9);
        x[14] ^= _rotl(x[13] + x[12], 13);
        x[15] ^= _rotl(x[14] + x[13], 18);
    }
    for (int i = 0; i < 16; ++i)
    {
        B[i] += x[i];
    }
}

static int salsa_idx[32][4] = {
    { 4, 0, 12, 7 },{ 9, 5, 1, 7 },{ 14, 10, 6, 7 },{ 3, 15, 11, 7 },
    { 8, 4, 0, 9 },{ 13, 9, 5, 9 },{ 2, 14, 10, 9 },{ 7, 3, 15, 9 },
    { 12, 8, 4, 13 },{ 1, 13, 9, 13 },{ 6, 2, 14, 13 },{ 11, 7, 3, 13 },
    { 0,12,8, 18 },{ 5,1,13, 18 },{ 10,6,2, 18 },{ 15,11,7, 18 },
    { 1,0,3, 7 },{ 6,5,4, 7 },{ 11,10,9, 7 },{ 12,15,14, 7 },
    { 2,1,0, 9 },{ 7,6,5, 9 },{ 8,11,10, 9 },{ 13,12,15, 9 },
    { 3,2,1, 13 },{ 4,7,6, 13 },{ 9,8,11, 13 },{ 14,13,12, 13 },
    { 0,3,2, 18 },{ 5,4,7, 18 },{ 10,9,8, 18 },{ 15,14,13, 18 }
};

void xor_salsa8_wd(struct WorkData* wd)
{
    uint32_t x[16];
    for (int i = 0; i < 16; ++i)
    {
        x[i] = (wd->B[i] ^= wd->Bx[i]);
    }
    

    for (int i = 0; i < 8; i += 2) {
        for (int j = 0; j < 32; ++j)
        {
            x[salsa_idx[j][0]] ^= _rotl(x[salsa_idx[j][1]] + x[salsa_idx[j][2]], salsa_idx[j][3]);
        }
    }
    for (int i = 0; i < 16; ++i)
    {
        wd->B[i] += x[i];
    }
}

void xor_salsa8_way4_ref(struct WorkSet* ws)
{
    xor_salsa8_wd(ws->WD + 0);
    xor_salsa8_wd(ws->WD + 1);
    xor_salsa8_wd(ws->WD + 2);
    xor_salsa8_wd(ws->WD + 3);
}

void xor_salsa8_way4_SSE(struct WorkSet* ws)
{
    __m128i bv[16];
    __m128i bxv[16];
    __m128i xv[16];
    uint32_t* b = (uint32_t*)bv;
    uint32_t* bx = (uint32_t*)bxv;
    uint32_t* x = (uint32_t*)xv;
    for (int v = 0; v < 16; ++v)
    {
        b[v * 4 + 0] = ws->WD[0].B[v]; 
        b[v * 4 + 1] = ws->WD[1].B[v]; 
        b[v * 4 + 2] = ws->WD[2].B[v]; 
        b[v * 4 + 3] = ws->WD[3].B[v];
        bx[v * 4 + 0] = ws->WD[0].Bx[v];
        bx[v * 4 + 1] = ws->WD[1].Bx[v];
        bx[v * 4 + 2] = ws->WD[2].Bx[v];
        bx[v * 4 + 3] = ws->WD[3].Bx[v];
    }
    for (int v = 0; v < 16; ++v)
    {
        xv[v] = bv[v] = _mm_xor_si128(bv[v], bxv[v]);
    }
    for (int i = 0; i < 8; i += 2) {
        for (int j = 0; j < 32; ++j)
        {
            __m128i _calc = _mm_add_epi32(xv[salsa_idx[j][1]], xv[salsa_idx[j][2]]);
            __m128i _shift_left = _mm_slli_epi32(_calc, salsa_idx[j][3]);
            xv[salsa_idx[j][0]] = _mm_xor_si128(xv[salsa_idx[j][0]], _shift_left);
            __m128i _shift_right = _mm_srli_epi32(_calc, (32 - salsa_idx[j][3]));
            xv[salsa_idx[j][0]] = _mm_xor_si128(xv[salsa_idx[j][0]], _shift_right);
        }
    }
    for (int v = 0; v < 16; ++v)
    {
        bv[v] = _mm_add_epi32(bv[v], xv[v]);
    }
    for (int v = 0; v < 16; ++v)
    {
        ws->WD[0].B[v] = b[v * 4 + 0];
        ws->WD[1].B[v] = b[v * 4 + 1];
        ws->WD[2].B[v] = b[v * 4 + 2];
        ws->WD[3].B[v] = b[v * 4 + 3];
    }
}

#pragma intrinsic(__rdtsc)

int scanhash_scrypt(int thr_id, struct work *work, uint32_t max_nonce, uint64_t *hashes_done, unsigned char *scratchbuf, uint32_t N);

int scanhash_scrypt_sse(int thr_id, struct work *work, uint32_t max_nonce, uint64_t *hashes_done, unsigned char *scratchbuf, uint32_t nn);
unsigned char *scrypt_buffer_alloc_sse(int nn);

void test_scanhash_scrypt()
{
    struct work scryptWork;
    uint32_t data[48] = { 0x06000000, 0x3ee4047a, 0x59863d03, 0x841dde00,
        0x335add65, 0x4854661a, 0x53c83f47, 0x20e5a23d,
        0x35dff734, 0xc6766d52, 0x1ff7508f, 0xdfc9cdcb,
        0x86139a08, 0x743c7ebe, 0x4a942414, 0xdc865136,
        0x70d9b711, 0x073d4e59, 0xff01011e, 0x80000006,
        0x80000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000280,
        0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000 };
    uint32_t target[8] = { 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00215534 };
    uint32_t hash[8] = { 0x1ddefac6, 0x37418d4c, 0xc4b82561, 0x2832c9c3,
        0xde93664e, 0x84d703e0, 0xc9501531, 0x0002720b };
    const int N = 1024 * 1024;
    unsigned char * scratchBuff = scrypt_buffer_alloc_sse(N);
    for (int i = 0; i < 4; ++i)
    {
        memset(&scryptWork, 0, sizeof(scryptWork));
        memcpy(scryptWork.data, data, sizeof(scryptWork.data));
        memcpy(scryptWork.target, target, sizeof(scryptWork.target));
        uint64_t hashesDone;
        const long long _Ctr1 = __rdtsc();;
        int res = scanhash_scrypt_sse(0, &scryptWork, 16, &hashesDone, scratchBuff, N);

        const long long _Ctr2 = __rdtsc();;
        long long delta = _Ctr2 - _Ctr1;
        printf("Test res %d, %lld, time %lld\n", res, hashesDone, delta / (1024 * 1024));
    }
}

int main(int argc, char *argv[])
{
	printf("Test started\n");

    test_transpose();
    int res = test_rotation();
    printf("rotation %d\n", res);

    const uint32_t B[16]  = { 0x6f54c89d, 0x5715059b, 0xc2a5624d, 0x0dc7b677, 0x87b2f312, 0xf62cf550, 0x53b282fd, 0xc8ff34d6, 0x1fc77a23, 0x814ecddc, 0x10642477, 0x1c181f5c, 0xdddabca8, 0x3bce9427, 0x32524048, 0xce512434 };
    const uint32_t Bx[16] = { 0x37b0a819, 0xf1deff63, 0x2f04fc79, 0x36997495, 0x26018ae6, 0x8ba55257, 0x595c23d2, 0x880d99c6, 0x9dfff6ce, 0x3504752c, 0x3df27f4d, 0x597aa991, 0xe20a335e, 0x04bae0d1, 0xdda7c4f8, 0x6ae01c71 };

    uint32_t t1_B[16], t2_B[16], t3_B[16];
    uint32_t t1_Bx[16], t2_Bx[16], t3_Bx[16];
    memcpy(t1_B, B, sizeof(B));
    memcpy(t2_B, B, sizeof(B));
    memcpy(t3_B, B, sizeof(B));
    memcpy(t1_Bx, Bx, sizeof(Bx));
    memcpy(t2_Bx, Bx, sizeof(Bx));
    memcpy(t3_Bx, Bx, sizeof(Bx));

    for (int wd = 0; wd < 16; ++wd) {
        printf("0x%04x,", B[wd]);
    }
    printf("\n");

    xor_salsa8(t1_B, t1_Bx);

    xor_salsa8_t(t2_B, t2_Bx);

    struct WorkData wd;
    memcpy(wd.B, B, sizeof(B));
    memcpy(wd.Bx, Bx, sizeof(Bx));
    xor_salsa8_wd(&wd);

    for (int wd = 0; wd < 16; ++wd) {
        printf("0x%04x,", B[wd]);
    }
    printf("\n");

    struct WorkSet ws1;
    memcpy(ws1.WD[0].B, B, sizeof(B)); memcpy(ws1.WD[0].Bx, Bx, sizeof(Bx));
    memcpy(ws1.WD[1].B, B, sizeof(B)); memcpy(ws1.WD[1].Bx, Bx, sizeof(Bx));
    memcpy(ws1.WD[2].B, B, sizeof(B)); memcpy(ws1.WD[2].Bx, Bx, sizeof(Bx));
    memcpy(ws1.WD[3].B, B, sizeof(B)); memcpy(ws1.WD[3].Bx, Bx, sizeof(Bx));
    
    struct WorkSet ws2;
    memcpy(&ws2, &ws1, sizeof(ws1));
    
    xor_salsa8_way4_ref(&ws1);

    if (0 == memcmp(t1_B, t2_B, sizeof(t1_B)) 
        && 0 == memcmp(t1_B, wd.B, sizeof(t1_B))
        && 0 == memcmp(t1_B, ws1.WD[0].B, sizeof(t1_B))
        && 0 == memcmp(t1_B, ws1.WD[1].B, sizeof(t1_B))
        && 0 == memcmp(t1_B, ws1.WD[2].B, sizeof(t1_B))
        && 0 == memcmp(t1_B, ws1.WD[3].B, sizeof(t1_B))
        )
    {
        printf("Test Ok !\n");
    }
    else
    {
        printf("Test Failed !\n");
    }

    xor_salsa8_way4_SSE(&ws2);
    if (0 == memcmp(&ws1, &ws2, sizeof(ws1)))
    {
        printf("SSE Ok !\n");
    }
    else
    {
        printf("SSE Failed !\n");
    }
    test_scanhash_scrypt();
    printf("Test stopped\n");
    _getch();
	return 0;
}
