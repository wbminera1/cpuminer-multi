#include "../../Catch/single_include/catch.hpp"
#include "../SIMD/vec4x32u.h"

TEST_CASE("SIMD", "[Vec2x32]")
{
    SECTION("Equ")
    {
        bool passed = true;
        {
            Vec4x32u v1(0);
            Vec4x32u v2(0);
            passed &= (v1 == v2).all_true();
        }
        {
            Vec4x32u v1(1,1,1,1);
            Vec4x32u v2(1,1,1,1);
            passed &= (v1 == v2).all_true();
        }
        {
            Vec4x32u v1(1, 2, 3, 4);
            Vec4x32u v2(1, 2, 3, 4);
            passed &= (v1 == v2).all_true();
        }
        {
            Vec4x32u v1(1, 1, 1, 1);
            Vec4x32u v2(0, 0, 0, 0);
            passed &= (v1 == v2).all_false();
        }
        {
            Vec4x32u v1(1, 2, 3, 4);
            Vec4x32u v2(4, 3, 2, 1);
            passed &= (v1 == v2).all_false();
        }
        REQUIRE(passed);
    }


	SECTION("Set")
	{
		bool passed = true;
        Vec4x32u v1(1,2,3,4);
        Vec4x32u v2(0);
        v2.load_i<0>(1);
        v2.load_i<1>(2);
        v2.load_i<2>(3);
        v2.load_i<3>(4);
        passed &= (v1 == v2).all_true();
        REQUIRE(passed);
	}

    SECTION("Add")
    {
        bool passed = true;
        Vec4x32u v1(1, 1, 1, 1);
        Vec4x32u v2(2, 2, 2, 2);
        Vec4x32u v3(3, 3, 3, 3);
        v2 += v1;
        passed &= (v2 == v3).all_true();

        REQUIRE(passed);
    }
}
