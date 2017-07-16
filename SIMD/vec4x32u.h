#ifndef VEC4x32_H_
#define VEC4x32_H_

#include <inttypes.h>

#include <xmmintrin.h>  // SSE
#include <emmintrin.h>  // SSE2
#include <pmmintrin.h>  // SSE3
#include <tmmintrin.h>  // SSSE3
#include <smmintrin.h>  // SSE4.1
#include <nmmintrin.h>  // SSE4.2
#include <immintrin.h>  // AVX

template <uint32_t i0, uint32_t i1, uint32_t i2, uint32_t i3>
static inline __m128i constant4i() {
    static const union {
        uint32_t  i[4];
        __m128i xmm;
    } u = { { i0,i1,i2,i3 } };
    return u.xmm;
}

static inline __m128i select(__m128i const & s, __m128i const & a, __m128i const & b) {
    return _mm_blendv_epi8(a, b, s);
}

class Vec4x32u {
protected:
    __m128i xmm; 
public:
    Vec4x32u() {
    }
    Vec4x32u(uint32_t i) {
        xmm = _mm_set1_epi32(i);
    }
    Vec4x32u(uint32_t i0, uint32_t i1, uint32_t i2, uint32_t i3) {
        xmm = _mm_setr_epi32(i0, i1, i2, i3);
    }
    Vec4x32u(const __m128i& x) {
        xmm = x;
    }
    Vec4x32u& operator=(const __m128i& x) {
        xmm = x;
        return *this;
    }
    operator __m128i() const {
        return xmm;
    }
    Vec4x32u& operator=(const uint32_t& i) {
        xmm = _mm_set1_epi32(i);
        return *this;
    }
    Vec4x32u& load_u(const Vec4x32u* p) {
        xmm = _mm_loadu_si128((const __m128i*)p);
        return *this;
    }
    Vec4x32u& load_a(const Vec4x32u* p) {
        xmm = _mm_load_si128((const __m128i*)p);
        return *this;
    }
    Vec4x32u& load_l(uint32_t i) {
        xmm = _mm_cvtsi32_si128(i);
        return *this;
    }
    template <size_t index>
    Vec4x32u& load_i(uint32_t i) {
        xmm = _mm_insert_epi32(xmm, (uint32_t)i, (const int)index);
        return *this;
    }
    void store_u(uint32_t* p) const {
        _mm_storeu_si128((__m128i*)p, xmm);
    }
    void store_a(uint32_t* p) const {
        _mm_store_si128((__m128i*)p, xmm);
    }
    uint32_t store_l() const {
        return (uint32_t)_mm_cvtsi128_si32(xmm);
    }
    uint32_t store_i(size_t index) const {
        uint32_t x[4];
        store_u(x);
        return x[index & 0x03];
    }
    static size_t size() {
        return sizeof(Vec4x32u);
    }
    Vec4x32u operator + (const Vec4x32u& a) {
        return _mm_add_epi32(xmm, a.xmm);
    }
    Vec4x32u& operator += (const Vec4x32u &a) {
        *this = *this + a;
        return *this;
    }
    Vec4x32u operator++ (int) {
        Vec4x32u a(xmm);
        *this = *this + 1;
        return a;
    }
    Vec4x32u& operator ++ () {
        *this = *this + 1;
        return *this;
    }
    Vec4x32u operator - (Vec4x32u const & a) {
        return _mm_sub_epi32(xmm, a.xmm);
    }
    Vec4x32u operator - () {
        return _mm_sub_epi32(_mm_setzero_si128(), xmm);
    }
    Vec4x32u& operator -= (const Vec4x32u& a) {
        *this = *this - a;
        return *this;
    }
    Vec4x32u operator -- (int) {
        Vec4x32u a(xmm);
        *this = *this - 1;
        return a;
    }
    Vec4x32u & operator -- () {
        *this = *this - 1;
        return *this;
    }
    Vec4x32u operator == (Vec4x32u const & a) {
        return _mm_cmpeq_epi32(xmm, a.xmm);
    }
/*
    Vec4x32u operator != (Vec4x32u const & a) {
        return Vec4x32u(Vec4x32u(~(*this == a)));
    }
*/
    Vec4x32u operator > (const Vec4x32u& a) {
        return _mm_cmpgt_epi32(xmm, a.xmm);
    }
    Vec4x32u operator < (const Vec4x32u& a) {
        return _mm_cmpgt_epi32(a.xmm, xmm);
    }
/*
    Vec4x32u operator >= (Vec4x32u const & a) {
        return Vec4x32u(Vec4x32u(~(b > a)));
    }
    Vec4x32u operator <= (Vec4x32u const & a, Vec4x32u const & b) {
        return b >= a;
    }
*/
    Vec4x32u compare(const Vec4x32u& a) const {
        return _mm_cmpeq_epi32(xmm, a.xmm);
    }

    int test(const Vec4x32u& a) const {
        return _mm_testc_si128(xmm, a);
    }

    int mask() const {
        return _mm_movemask_epi8(xmm);
    }

    bool all_false() const {
        return 0 == mask();
    }

    bool all_true() const {
        return 0xFFFF == mask();
    }
};


#endif // VEC4x32_H_
