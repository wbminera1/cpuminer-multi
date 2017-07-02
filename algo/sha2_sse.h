#pragma once
#include <immintrin.h>
#include <stdint.h>

extern const uint32_t sha256_h[8];
extern const uint32_t sha256_k[64];

extern const __m128i sha256_h_sse[8];
extern __m128i sha256_k_sse[64];

inline void sha256_init_sse(__m128i* state)
{
    for (size_t i = 0; i < 8; ++i)
    {
        register int val = (int)sha256_h[i];
        state[i] = _mm_set1_epi32(val);
    }
}

/* Elementary functions used by SHA256 */
inline uint32_t Ch(uint32_t x, uint32_t y, uint32_t z) {
    return ((x & (y ^ z)) ^ z);
}
inline uint32_t Maj(uint32_t x, uint32_t y, uint32_t z) {
    return ((x & (y | z)) | (y & z));
}

inline uint32_t S0(uint32_t x) { return (_rotr(x, 2) ^ _rotr(x, 13) ^ _rotr(x, 22)); }
inline uint32_t S1(uint32_t x) { return (_rotr(x, 6) ^ _rotr(x, 11) ^ _rotr(x, 25)); }
inline uint32_t s0(uint32_t x) { return (_rotr(x, 7) ^ _rotr(x, 18) ^ (x >> 3)); }
inline uint32_t s1(uint32_t x) { return (_rotr(x, 17) ^ _rotr(x, 19) ^ (x >> 10)); }

/* SHA256 round function */
inline void RND(uint32_t a, uint32_t b, uint32_t c, uint32_t* d, uint32_t e, uint32_t f, uint32_t g, uint32_t* h, uint32_t k)
{
    uint32_t t0 = *h + S1(e) + Ch(e, f, g) + k;
    uint32_t t1 = S0(a) + Maj(a, b, c);
    *d += t0;
    *h = t0 + t1;
}

/* Adjusted round function for rotating state */
inline void RNDr(uint32_t* S, uint32_t*W, size_t i) {
    RND(S[(64 - i) % 8], S[(65 - i) % 8], S[(66 - i) % 8], &S[(67 - i) % 8], 
        S[(68 - i) % 8], S[(69 - i) % 8], S[(70 - i) % 8], &S[(71 - i) % 8], 
        W[i] + sha256_k[i]);
}

inline void sha256_transform_t(uint32_t *state, const uint32_t *block, int swap)
{
    uint32_t W[64];
    uint32_t S[8];

    /* 1. Prepare message schedule W. */
    if (swap) {
        for (size_t i = 0; i < 16; i++)
            W[i] = swab32(block[i]);
    }
    else
        memcpy(W, block, sizeof(uint32_t) * 16);
    for (size_t i = 16; i < 64; i += 2) {
        W[i] = s1(W[i - 2]) + W[i - 7] + s0(W[i - 15]) + W[i - 16];
        W[i + 1] = s1(W[i - 1]) + W[i - 6] + s0(W[i - 14]) + W[i - 15];
    }

    /* 2. Initialize working variables. */
    memcpy(S, state, sizeof(S));

    /* 3. Mix. */
    for (size_t i = 0; i < 64; ++i) {
        RNDr(S, W, i);
    }

    /* 4. Mix local working variables into global state */
    for (size_t i = 0; i < 8; i++)
        state[i] += S[i];
}


inline __m128i XOR_sse(__m128i x, __m128i y)
{
    return _mm_xor_si128(x, y);
}

inline __m128i XOR3_sse(__m128i x, __m128i y, __m128i z)
{
    return _mm_xor_si128(_mm_xor_si128(x, y), z);
}

inline __m128i OR_sse(__m128i x, __m128i y)
{
    return _mm_or_si128(x, y);
}

inline __m128i AND_sse(__m128i x, __m128i y)
{
    return _mm_and_si128(x, y);
}

inline __m128i ADD_sse(__m128i x, __m128i y)
{
    return _mm_add_epi32(x, y);
}

inline __m128i ADD4_sse(__m128i a1, __m128i a2, __m128i a3, __m128i a4)
{
    return _mm_add_epi32(_mm_add_epi32(a1, a2), _mm_add_epi32(a3, a4));
}

inline __m128i ROTR_sse(__m128i x, int s)
{
    return _mm_or_si128(_mm_slli_epi32(x, 32 - s), _mm_srli_epi32(x, s));
}

inline __m128i SHR_sse(__m128i x, int s)
{
    return _mm_srli_epi32(x, s);
}

inline __m128i Ch_sse(__m128i x, __m128i y, __m128i z) {
    return XOR_sse(AND_sse(x,XOR_sse(y,z)),z);
}
inline __m128i Maj_sse(__m128i x, __m128i y, __m128i z) {
    return OR_sse(AND_sse(x,OR_sse(y,z)),AND_sse(y,z));
}

inline __m128i S0_sse(__m128i x) { return XOR3_sse(ROTR_sse(x, 2), ROTR_sse(x, 13), ROTR_sse(x, 22)); }
inline __m128i S1_sse(__m128i x) { return XOR3_sse(ROTR_sse(x, 6), ROTR_sse(x, 11), ROTR_sse(x, 25)); }
inline __m128i s0_sse(__m128i x) { return XOR3_sse(ROTR_sse(x, 7), ROTR_sse(x, 18), SHR_sse(x, 3)); }
inline __m128i s1_sse(__m128i x) { return XOR3_sse(ROTR_sse(x, 17), ROTR_sse(x, 19),SHR_sse(x, 10)); }

inline void RND_sse(__m128i a, __m128i b, __m128i c, __m128i* d, __m128i e, __m128i f, __m128i g, __m128i* h, __m128i k)
{
    __m128i t0 = ADD4_sse(*h, S1_sse(e), Ch_sse(e, f, g), k);
    __m128i t1 = ADD_sse(S0_sse(a), Maj_sse(a, b, c));
    *d = ADD_sse(*d, t0);
    *h = ADD_sse(t0,t1);
}

inline void RNDr_sse(__m128i* S, __m128i* W, size_t i) {
    RND_sse(S[(64 - i) % 8], S[(65 - i) % 8], S[(66 - i) % 8], &S[(67 - i) % 8],
        S[(68 - i) % 8], S[(69 - i) % 8], S[(70 - i) % 8], &S[(71 - i) % 8],
        ADD_sse(W[i], sha256_k_sse[i]));
}

static inline __m128i swab32_sse(__m128i v)
{
    uint32_t* v4 = (uint32_t*)&v;
    for (size_t i = 0; i < 4; ++i)
    {
        v4[i] = swab32(v4[i]);
    }
    return v;
}

inline void sha256_transform_sse(__m128i* state, const __m128i* block, int swap)
{
    __m128i W[64];
    __m128i S[8];

    /* 1. Prepare message schedule W. */
    if (swap) {
        for (size_t i = 0; i < 16; i++)
            W[i] = swab32_sse(block[i]);
    }
    else
        memcpy(W, block, sizeof(__m128i) * 16);
    for (size_t i = 16; i < 64; i += 2) {
        W[i] = ADD4_sse(s1_sse(W[i - 2]), W[i - 7], s0_sse(W[i - 15]), W[i - 16]);
        W[i + 1] = ADD4_sse(s1_sse(W[i - 1]), W[i - 6], s0_sse(W[i - 14]), W[i - 15]);
    }

    /* 2. Initialize working variables. */
    memcpy(S, state, sizeof(S));

    /* 3. Mix. */
    for (size_t i = 0; i < 64; ++i) {
        RNDr_sse(S, W, i);
    }

    /* 4. Mix local working variables into global state */
    for (size_t i = 0; i < 8; i++)
        state[i] = ADD_sse(state[i], S[i]);
}
