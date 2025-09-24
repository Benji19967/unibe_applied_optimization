#include <iostream>
#include <Utils/StopWatch.hh>
#include <Algorithms/GridSearch.hh>

#include "gtest/gtest.h"


using namespace AOPT;

/** Those unit tests basically checks individual parts of your implementation.
 * They generally work by either
 * 1. comparing your results with ours
 * 2. comparing your results with manually obtained ones
 *    (when it's a simple problem solvable by hand)
 *
 * A good way to work with them is to make each test pass successively.
 * For instance, if a test checks the results of a function, its gradient and its hessian,
 * you can:
 * 1. implement the function itself (eval_f())
 * 2. check that the test's assertion passes (e.g. ASSERT_EQ(eval_f(x), expected_value_of_f(x)) passes)
 * 3. implement the gradient
 * 4. check that the gradient-related assertion passes
 * 5. implement the hessian
 * 6. check that the whole test passes
 *
 * Feel free to modify the tests but ONLY to output intermediary values/results BUT
 * do not modify it in any way that would change its result.
 *
 * For more details about googletest, please visit
 * https://github.com/Macaulay2/googletest-1/blob/master/googletest/docs/Primer.md
**/



TEST(FunctionsTest, QuadraticFucntion2d){
    FunctionQuadratic2D func(1.5);
    FunctionQuadratic2D::Vec x(2);
    x<<-1, 2;

    ASSERT_EQ(func.eval_f(x), 3.5);
}


TEST(FunctionsTest, QuadraticFucntionNd){

    const int n(3);

    FunctionQuadraticND::Mat A(n,n);
    A<<     -1, 2, 3,
            0,-3,4.5,
            4,2, -5;

    FunctionQuadraticND::Vec b(n);
    b<<-1,0,-1;

    const double c(-1);


    FunctionQuadraticND func(A,b,c);
    FunctionQuadratic2D::Vec x(n);
    x<<-1,1,-1;

    ASSERT_EQ(func.eval_f(x), -4.25);
}


TEST(FunctionsTest, NonConvexFunction2d){

    FunctionNonConvex2D func;
    FunctionQuadratic2D::Vec x(2);
    x<<0.5, -0.5;

    ASSERT_FLOAT_EQ(func.eval_f(x), 1.1057945);
}



TEST(GridSearchTests, Grid2dFunction2dGridSearch) {

    const int n(5);
    GridSearch grid(n);

    FunctionQuadratic2D func(-0.8);

    FunctionQuadratic2D::Vec x_l(2);
    x_l<<-1, -1;
    FunctionQuadratic2D::Vec x_u(2);
    x_u<<1, 1;

    double f_min;
    grid.grid_search_2d(&func, x_l, x_u, f_min);

    ASSERT_EQ(f_min, -0.38);
}

TEST(GridSearchTests, GridNdFunctionNdGridSearch) {

    const int grid_n(10);
    GridSearch grid(grid_n);

    const int func_n(5);

    FunctionQuadraticND::Mat A(func_n, func_n);
    A<<     -1, 2, 3,     0, -3,
             0,-3, 4.5, 1.5, -0.5,
             4, 2,-5,     0, -2,
            -4,-2, 4,  -3.5,  3,
             2,-4, 2,   6.5,  1;


    FunctionQuadraticND::Vec b(func_n);
    b<< 2, -1, 2, 3, -5;

    const double c(-1);
    FunctionQuadraticND func(A, b, c);

    FunctionQuadraticND::Vec x_l = FunctionQuadraticND::Vec::Zero(func_n);
    FunctionQuadraticND::Vec x_u = FunctionQuadraticND::Vec::Ones(func_n) * 4;

    double f_min;
    grid.grid_search_nd(&func, x_l, x_u, f_min);

    ASSERT_FLOAT_EQ(f_min, -77);
}

int main(int _argc, char** _argv){

    testing::InitGoogleTest(&_argc, _argv);
    return RUN_ALL_TESTS();

}
