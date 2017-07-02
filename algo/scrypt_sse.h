#ifndef __SCRYPT_SSE
#define __SCRYPT_SSE
#include <inttypes.h>


#define N (1024 * 1024)

struct ScryptData
{
    uint32_t X[4][32];
    uint32_t data[4][20 + 12];
    uint32_t hash[4][8];
    uint32_t midstate[4][8];
    uint32_t tstate[4][8];
    uint32_t ostate[4][8];
    unsigned char* scratchbuf[4];
    uint32_t* V[4];
#ifndef _TEST
    uint32_t* B[4];
    uint32_t* Bx[4];
#endif _TEST
};

struct ScryptDataSSE
{
    __m128i X[32];
    __m128i data[20 + 12];
    __m128i hash[8];
    __m128i midstate[8];
    __m128i tstate[8];
    __m128i ostate[8];
    __m128i* V;
    unsigned char* scratchbuf;
};
#endif // __SCRYPT_SSE
