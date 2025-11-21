#include <iostream>
#include <Utils/StopWatch.hh>
#include <Algorithms/ConvexityTest.hh>
#include <Functions/FunctionQuadratic2D.hh>
#include <Functions/FunctionQuadraticND.hh>
#include <Functions/FunctionNonConvex2D.hh>

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


TEST(ConvexityTest, Quadratic2DisConvex){
    AOPT::FunctionQuadratic2D function(1.345);
    
    const double min(-10), max(10);
    const int n_evals(20);
    ASSERT_TRUE(AOPT::ConvexityTest::isConvex(&function, min, max, n_evals));
}

TEST(ConvexityTest, QuadraticNDisConvex){
    
    const int n(4);
    AOPT::FunctionQuadraticND function(n);
    
    const double min(-100), max(100);
    const int n_evals(20);
    ASSERT_TRUE(AOPT::ConvexityTest::isConvex(&function, min, max, n_evals));
}


TEST(ConvexityTest, OtherQuadraticNDisConvex){
    
    const int n(6);
    AOPT::FunctionQuadraticND function(n);
    
    const double min(-10), max(10);
    const int n_evals(10);
    ASSERT_TRUE(AOPT::ConvexityTest::isConvex(&function, min, max, n_evals));
}


TEST(ConvexityTest, YetAnotherQuadraticNDisConvex){
    
    const int n(8);
    AOPT::FunctionQuadraticND function(n);
    
    const double min(-20), max(20);
    const int n_evals(10);
    ASSERT_TRUE(AOPT::ConvexityTest::isConvex(&function, min, max, n_evals));
}


TEST(ConvexityTest, NonConvex2DIsConvexOnCertainInterval){
    
    AOPT::FunctionNonConvex2D function;
    
    const double min(1.0), max(1.001);
    const int n_evals(10);
    ASSERT_TRUE(AOPT::ConvexityTest::isConvex(&function, min, max, n_evals));
}


TEST(ConvexityTest, NonConvex2DIsNonConvexOnOtherInterval){
    
    AOPT::FunctionNonConvex2D function;
    
    const double min(-10), max(10);
    const int n_evals(10);
    ASSERT_FALSE(AOPT::ConvexityTest::isConvex(&function, min, max, n_evals));
}


TEST(ConvexityTest, NonConvex2DIsNonConvexOnYetAnotherInterval){
    
    AOPT::FunctionNonConvex2D function;
    
    const double min(10), max(20);
    const int n_evals(10);
    ASSERT_FALSE(AOPT::ConvexityTest::isConvex(&function, min, max, n_evals));
}


int main(int _argc, char** _argv){

    testing::InitGoogleTest(&_argc, _argv);
    return RUN_ALL_TESTS();

}

