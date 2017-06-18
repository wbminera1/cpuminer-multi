#include <stdio.h>
#include <stdint.h>
#include <conio.h>
#include <memory.h>
#include <immintrin.h>
#include <inttypes.h>


void xor_salsa8(uint32_t B[16], const uint32_t Bx[16]);

inline void swap(uint32_t* a, uint32_t* b)
{
    uint32_t tmp = *a;
    *a = *b;
    *b = tmp;
}

void transpose(uint32_t* m)
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

/*static inline */void xor_salsa8_t(uint32_t B[16], const uint32_t Bx[16])
{
    //uint32_t x00, x01, x02, x03, x04, x05, x06, x07, x08, x09, x10, x11, x12, x13, x14, x15;
    uint32_t x[16];
    __m128i* sX = (__m128i*)x;
    __m128i* sB = (__m128i*)B;
    __m128i* sBx = (__m128i*)Bx;
    int i;
    sX[0] = sB[0] = _mm_xor_si128(sB[0], sBx[0]);
    sX[1] = sB[1] = _mm_xor_si128(sB[1], sBx[1]);
    sX[2] = sB[2] = _mm_xor_si128(sB[2], sBx[2]);
    sX[3] = sB[3] = _mm_xor_si128(sB[3], sBx[3]);

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

    for (i = 0; i < 8; i += 2) {
#define R(a, b) (((a) << (b)) | ((a) >> (32 - (b))))
        // Operate on columns.

        transpose(x);

        x[1] ^= R(x[0] + x[3], 7);
        x[2] ^= R(x[0] + x[1], 9);
        x[3] ^= R(x[2] + x[1], 13);
        x[0] ^= R(x[2] + x[3], 18);

        x[6] ^= R(x[5] + x[4], 7);
        x[7] ^= R(x[5] + x[6], 9);
        x[4] ^= R(x[7] + x[6], 13);
        x[5] ^= R(x[7] + x[4], 18);

        x[11] ^= R(x[10] + x[9], 7);
        x[8] ^= R(x[10] + x[11], 9);
        x[9] ^= R(x[8] + x[11], 13);
        x[10] ^= R(x[8] + x[9], 18);

        x[12] ^= R(x[15] + x[14], 7);
        x[13] ^= R(x[15] + x[12], 9);
        x[14] ^= R(x[13] + x[12], 13);
        x[15] ^= R(x[13]+ x[14], 18);

        transpose(x);

        /* Operate on rows. */
        x[1] ^= R(x[0] + x[3], 7);
        x[2] ^= R(x[1] + x[0], 9);
        x[3] ^= R(x[2] + x[1], 13);
        x[0] ^= R(x[3] + x[2], 18);

        x[6] ^= R(x[5] + x[4], 7);
        x[7] ^= R(x[6] + x[5], 9);
        x[4] ^= R(x[7] + x[6], 13);
        x[5] ^= R(x[4] + x[7], 18);

        x[11] ^= R(x[10] + x[9], 7);
        x[8] ^= R(x[11] + x[10], 9);
        x[9] ^= R(x[8] + x[11], 13);
        x[10] ^= R(x[9] + x[8], 18);

        x[12] ^= R(x[15] + x[14], 7);
        x[13] ^= R(x[12] + x[15], 9);
        x[14] ^= R(x[13] + x[12], 13);
        x[15] ^= R(x[14] + x[13], 18);
#undef R
    }
    B[0] += x[ 0];
    B[1] += x[ 1];
    B[2] += x[ 2];
    B[3] += x[ 3];
    B[4] += x[ 4];
    B[5] += x[ 5];
    B[6] += x[ 6];
    B[7] += x[ 7];
    B[8] += x[ 8];
    B[9] += x[ 9];
    B[10] += x[10];
    B[11] += x[11];
    B[12] += x[12];
    B[13] += x[13];
    B[14] += x[14];
    B[15] += x[15];
}

int main(int argc, char *argv[])
{
	printf("Test started\n");

    test_transpose();

    const uint32_t B[16]  = { 0x6f54c89d, 0x5715059b, 0xc2a5624d, 0x0dc7b677, 0x87b2f312, 0xf62cf550, 0x53b282fd, 0xc8ff34d6, 0x1fc77a23, 0x814ecddc, 0x10642477, 0x1c181f5c, 0xdddabca8, 0x3bce9427, 0x32524048, 0xce512434 };
    const uint32_t Bx[16] = { 0x37b0a819, 0xf1deff63, 0x2f04fc79, 0x36997495, 0x26018ae6, 0x8ba55257, 0x595c23d2, 0x880d99c6, 0x9dfff6ce, 0x3504752c, 0x3df27f4d, 0x597aa991, 0xe20a335e, 0x04bae0d1, 0xdda7c4f8, 0x6ae01c71 };

    uint32_t t1_B[16], t2_B[16];
    uint32_t t1_Bx[16], t2_Bx[16];
    memcpy(t1_B, B, sizeof(B));
    memcpy(t2_B, B, sizeof(B));
    memcpy(t1_Bx, Bx, sizeof(Bx));
    memcpy(t2_Bx, Bx, sizeof(Bx));

    for (int wd = 0; wd < 16; ++wd) {
        printf("0x%04x,", B[wd]);
    }
    printf("\n");

    xor_salsa8(t1_B, t1_Bx);

    xor_salsa8_t(t2_B, t2_Bx);

    for (int wd = 0; wd < 16; ++wd) {
        printf("0x%04x,", B[wd]);
    }
    printf("\n");

    if (0 == memcmp(t1_B, t2_B, sizeof(t1_B)))
    {
        printf("Test Ok !");
    }
    else
    {
        printf("Test Failed !");
    }

    _getch();
	printf("Test stopped\n");
	return 0;
}
