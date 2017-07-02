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
inline void uint32_t_arr_to_mi128(__m128i* _dst, const uint32_t* src, size_t count128)
{
    uint32_t* dst = (uint32_t*)_dst;
    for (size_t i = 0; i < count128; ++i)
    {
        dst[i * 4 + 0] = src[i + count128 * 0];
        dst[i * 4 + 1] = src[i + count128 * 1];
        dst[i * 4 + 2] = src[i + count128 * 2];
        dst[i * 4 + 3] = src[i + count128 * 3];
    }
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