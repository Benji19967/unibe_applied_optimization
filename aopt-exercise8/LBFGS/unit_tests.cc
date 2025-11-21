#include <iostream>

#include <Utils/StopWatch.hh>
#include <Utils/RandomNumberGenerator.hh>

#include <Functions/FunctionQuadratic2D.hh>
#include <Functions/FunctionQuadraticND.hh>

#include <Algorithms/LBFGS.hh>


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




TEST(LBFGS, SimpleQuadraticProblem){

    const int dim(5);


    LBFGS::Mat A(dim, dim);
    A <<    246.652,  107.143,  117.078, -125.157,  21.0117,
            107.143,  244.831, -34.6749, -3.39497,  105.168,
            117.078, -34.6749,  114.209, -57.6975,  14.9266,
           -125.157, -3.39497, -57.6975,  129.647,   99.809,
            21.0117,  105.168,  14.9266,   99.809,  222.531;

    LBFGS::Vec b(dim);
    b << 0.832402, -9.81436,  9.98931, -9.74196,  6.84995;

    double c(0);


    FunctionQuadraticND func(A, b, c);

    LBFGS lbfgs(3);
    LBFGS::Vec init_x(dim);
    init_x<<-1, 5, 4,3,5;

    auto result = lbfgs.solve(&func, init_x);
    std::cout<<result.transpose()<<std::endl;

    LBFGS::Vec expected_result(dim);
    expected_result<<0.854001, -0.123255,  -0.33051,   1.18085, -0.560632;

    ASSERT_LT((result - expected_result).norm(), 1e-5);
}


TEST(LBFGS, AnotherSimpleQuadraticProblem){

    const int dim(4);


    LBFGS::Mat A(dim, dim);
    A <<    193.244,  51.9043, -54.0949, -121.182,
            51.9043,  281.143, -64.3568, -1.16595,
           -54.0949, -64.3568,  125.882,  27.6244,
           -121.182, -1.16595,  27.6244,   131.61;

    LBFGS::Vec b(dim);
    b<<  4.71305, -7.69894,  3.97317, -2.8879;


    double c(0);

    FunctionQuadraticND func(A, b, c);

    LBFGS lbfgs(3);
    LBFGS::Vec init_x(dim);
    init_x<<1, -2, 3, -4;

    auto result = lbfgs.solve(&func, init_x);

    std::cout<<"result: "<<std::setprecision(10)<<result.transpose()<<std::endl;

    LBFGS::Vec expected_result(dim);
    expected_result<< -0.05604022015,   0.0294477851, -0.03578660566,    -0.02188488;

    ASSERT_LT((result - expected_result).norm(), 1e-5);
}




int main(int _argc, char** _argv){

    testing::InitGoogleTest(&_argc, _argv);
    return RUN_ALL_TESTS();

}
