#include "miner.h"

#include <string.h>
#include <stdint.h>

#include "crypto/blake2s.h"

void blake2s_hash(void *output, const void *input)
{
	unsigned char hash[128] = { 0 };
	blake2s_state blake2_ctx;

	blake2s_init(&blake2_ctx, BLAKE2S_OUTBYTES);
	blake2s_update(&blake2_ctx, input, 80);
	blake2s_final(&blake2_ctx, hash, BLAKE2S_OUTBYTES);

	memcpy(output, hash, 32);
}

int scanhash_blake2s(int thr_id, struct work *work, uint32_t max_nonce, uint64_t *hashes_done)
{
	uint32_t _ALIGN(128) hash32[8];
	uint32_t _ALIGN(128) endiandata[20];
	uint32_t *pdata = work->data;
	uint32_t *ptarget = work->target;

	const uint32_t Htarg = ptarget[7];
	const uint32_t first_nonce = pdata[19];

	uint32_t n = first_nonce;

	for (int i=0; i < 19; i++) {
		be32enc(&endiandata[i], pdata[i]);
	}

	do {
		be32enc(&endiandata[19], n);
		blake2s_hash(hash32, endiandata);
		if (hash32[7] < Htarg && fulltest(hash32, ptarget)) {
			work_set_target_ratio(work, hash32);
			*hashes_done = n - first_nonce + 1;
			pdata[19] = n;
			return 1;
		}
		n++;

	} while (n < max_nonce && !work_restart[thr_id].restart);

	*hashes_done = n - first_nonce + 1;
	pdata[19] = n;

	return 0;
}
