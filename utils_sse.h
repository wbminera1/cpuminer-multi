#pragma once

// convert uint32_t arr[count * 4] to mi128[count]
inline void uint32_t_to_mi128(__m128i* _dst, const uint32_t* src, size_t count128)
{
    if ((count128 % 4) == 0)
    {
        uint32_t* dst = (uint32_t*)_dst;
        for (size_t i = 0; i < count128; ++i)
        {
            dst[i * 4 + 0] = src[i + 0];
            dst[i * 4 + 1] = src[i + 4];
            dst[i * 4 + 2] = src[i + 8];
            dst[i * 4 + 3] = src[i + 12];
        }
    }
}

// convert uint32_t arr[4][count] to mi128[count]
inline void uint32_t_arr_to_mi128(__m128i* dst, const uint32_t* src, size_t count128)
{
    /*
    uint32_t* dst = (uint32_t*)_dst;
    for (size_t i = 0; i < count128; ++i)
    {
        dst[i * 4 + 0] = src[i + count128 * 0];
        dst[i * 4 + 1] = src[i + count128 * 1];
        dst[i * 4 + 2] = src[i + count128 * 2];
        dst[i * 4 + 3] = src[i + count128 * 3];
    }
    */
    memcpy(dst, src, sizeof(__m128i) * count128);
}

// convert mi128[count] to uint32_t arr[count * 4]
inline void mi128_to_uint32_t(uint32_t* dst, const __m128i* _src, size_t count128)
{
    if ((count128 % 4) == 0)
    {
        uint32_t* src = (uint32_t*)_src;
        for (size_t i = 0; i < count128; ++i)
        {
            dst[i + 0] = src[i * 4 + 0];
            dst[i + 4] = src[i * 4 + 1];
            dst[i + 8] = src[i * 4 + 2];
            dst[i + 12] = src[i * 4 + 3];
        }
    }
}

// convert mi128[count] to uint32_t arr[4][count]
inline void mi128_arr_to_uint32_t(uint32_t* dst, const __m128i* src, size_t count128)
{
    memcpy(dst, src, sizeof(__m128i) * count128);
}

inline void uint32_t_set_to_mi128(__m128i* dst, const uint32_t* src, size_t count128)
{
    for (size_t i = 0; i < count128; ++i)
    {
        dst[i] = _mm_set1_epi32(src[i]);
    }
}

// copy low part of src to dst
inline void mi128_convert_to_uint32_t(uint32_t* dst, const __m128i* src, size_t count128)
{
    for (size_t i = 0; i < count128; ++i)
    {
        dst[i] = (uint32_t)_mm_cvtsi128_si32(src[i]);
    }
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
inline __m128i INC_sse(__m128i x)
{
    return _mm_add_epi32(x, _mm_set1_epi32(1));
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
