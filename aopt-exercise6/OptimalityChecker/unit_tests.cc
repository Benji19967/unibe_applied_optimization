#include <Utils/OptimalityChecker.hh>
#include <Functions/FunctionQuadraticND.hh>
#include <vector>
#include <iostream>
#include "gtest/gtest.h"


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



TEST(QuadraticNDTest, gradientEvaluation){

    const int n(4);
    AOPT::FunctionQuadraticND::Mat A(n,n);
    AOPT::FunctionQuadraticND::Vec b(n);

    A << 1, -4,  0.5, 5,
        -4,  4,  3,  -6,
         2,  0, -2,   0,
         1, -4,  0,   5;
    b << 1,  1, -3,   4;
    AOPT::FunctionQuadraticND obj_func(A, b, 5);

    AOPT::FunctionQuadraticND::Vec x(n);
    x<<1,-1,-1,0;

    AOPT::FunctionQuadraticND::Vec g(n);
    obj_func.eval_gradient(x, g);

    AOPT::FunctionQuadraticND::Vec expected_g(n);
    expected_g << 4.75,  -8.5, -1.25, 12;

    ASSERT_EQ(g, expected_g);

}


/** This "test fixture" allows to have a zero-initiazlied Matrix variable A and
 * a Vector variable b accessible in each test of the fixture, without having
 * to re-declare them at the beginning. This makes the code overall cleaner */
class KKTConditionsTest : public ::testing::Test{

protected:

      void SetUp() override {
          A = AOPT::FunctionQuadraticND::Mat(2,2);
          b = AOPT::FunctionQuadraticND::Vec(2);

          A.setZero();
          b.setZero();
      }

      AOPT::FunctionQuadraticND::Mat A;
      AOPT::FunctionQuadraticND::Vec b;
};



TEST_F(KKTConditionsTest, worksWithEmptyConstraints){

    //set up the optimization problem
    //-------------------------------------------------------------------------------//

    //1. set objective function
    A << 1, 0,
         0, 1;
    b << 1, 1;
    AOPT::FunctionQuadraticND obj_func(A, b, 0);

    //2. inequality constraints
    std::vector<AOPT::FunctionBase*> ineq_cons;

    //3. equality constraints
    std::vector<AOPT::FunctionBase*> eq_cons;


    //4. set lambdas and vs
    AOPT::FunctionQuadraticND::Vec lambda(2), nu(1);

    //5. set query point
    AOPT::OptimalityChecker::Vec x(2);
    x<<-1, -1;


    //-------------------------------------------------------------------------------//

    //check
    AOPT::OptimalityChecker oc;
    ASSERT_TRUE(oc.is_KKT_satisfied(&obj_func, ineq_cons, eq_cons, x, lambda, nu));
}




TEST_F(KKTConditionsTest, falseIfInequalityConstraintsAreViolated){

    //set up the optimization problem
    //-------------------------------------------------------------------------------//

    //1. set objective function
    A << 1, 0,
         0, 1;
    b << 1, 1;
    AOPT::FunctionQuadraticND obj_func(A, b, 0);

    //2. inequality constraints
    std::vector<AOPT::FunctionBase*> ineq_cons;

    A << -1, 1,
         2, -3;
    b << -4, 5;
    AOPT::FunctionQuadraticND qnd(A, b, 0);
    ineq_cons.push_back(&qnd);


    //3. equality constraints
    std::vector<AOPT::FunctionBase*> eq_cons;


    //4. set lambdas and vs
    AOPT::FunctionQuadraticND::Vec lambda(1), nu(0);

    //5. set query point
    AOPT::OptimalityChecker::Vec x(2);
    x<< -1, 1;


    //-------------------------------------------------------------------------------//

    //check
    AOPT::OptimalityChecker oc;
    ASSERT_FALSE(oc.is_KKT_satisfied(&obj_func, ineq_cons, eq_cons, x, lambda, nu));

}



TEST_F(KKTConditionsTest, falseIfEqualityConstraintsAreViolated){

    //set up the optimization problem
    //-------------------------------------------------------------------------------//

    //1. set objective function
    A << 4, 3,
         -3, 5;
    b << 3, -3;
    AOPT::FunctionQuadraticND obj_func(A, b, 0);

    //2. inequality constraints
    std::vector<AOPT::FunctionBase*> ineq_cons;


    //3. equality constraints
    std::vector<AOPT::FunctionBase*> eq_cons;
    A << -4, 3,
         5, -6;
    b << 0, 3;
    AOPT::FunctionQuadraticND qnd3(A, b, 0);
    eq_cons.push_back(&qnd3);

    //4. set lambdas and vs
    AOPT::FunctionQuadraticND::Vec lambda(0), nu(1);

    //5. set query point
    AOPT::OptimalityChecker::Vec x(2);
    x<< -1, 1;


    //-------------------------------------------------------------------------------//

    //check
    AOPT::OptimalityChecker oc;
    ASSERT_FALSE(oc.is_KKT_satisfied(&obj_func, ineq_cons, eq_cons, x, lambda, nu));
}




TEST_F(KKTConditionsTest, falseIfComplementarySlacknessIsViolated){

    //set up the optimization problem
    //-------------------------------------------------------------------------------//

    //1. set objective function
    A << -5, 2,
          0, 4;
    b << 1, -1;
    AOPT::FunctionQuadraticND obj_func(A, b, 0);

    //2. inequality constraints
    std::vector<AOPT::FunctionBase*> ineq_cons;

    A <<  1,  1,
         -2, -3;
    b << -7, 4;
    AOPT::FunctionQuadraticND qnd(A, b, 0);
    ineq_cons.push_back(&qnd);

    A <<  0,  -4,
         -5, 3;
    b << 2, -4;
    AOPT::FunctionQuadraticND qnd2(A, b, 0);
    ineq_cons.push_back(&qnd2);



    //3. equality constraints
    std::vector<AOPT::FunctionBase*> eq_cons;


    //4. set lambdas and vs
    AOPT::FunctionQuadraticND::Vec lambda(2), nu(0);
    lambda[0] = 1;
    lambda[1] = 0.5;

    //5. set query point
    AOPT::OptimalityChecker::Vec x(2);
    x<< 1, 1;


    //-------------------------------------------------------------------------------//

    //check
    AOPT::OptimalityChecker oc;
    ASSERT_FALSE(oc.is_KKT_satisfied(&obj_func, ineq_cons, eq_cons, x, lambda, nu));

}



/** This minimization problem is taken from the lecture's slides.
 * This example was shown to satisfy the KKT conditions and is_KKT_satisfied
 * should thus return true when given this problem. */
TEST_F(KKTConditionsTest, conditionsRespectedWithLectureExample){

    //set up the optimization problem
    //-------------------------------------------------------------------------------//

    //1. set objective function
    A.setZero();
    b << 1, 0;
    AOPT::FunctionQuadraticND obj_func(A, b, 0);

    //2. inequality constraints
    std::vector<AOPT::FunctionBase*> ineq_cons;

    A <<  2, 0,
          0, 2;
    b << 0, 0;
    AOPT::FunctionQuadraticND qnd(A, b, -1);
    ineq_cons.push_back(&qnd);


    A.setZero();
    b << 1, 0;
    AOPT::FunctionQuadraticND qnd2(A, b, 0);
    ineq_cons.push_back(&qnd2);


    //3. equality constraints
    std::vector<AOPT::FunctionBase*> eq_cons;



    //4. set lambdas and vs
    AOPT::FunctionQuadraticND::Vec lambda(2), nu(0);
    lambda[0] = 0.5;
    lambda[1] = 0;

    //5. set query point
    AOPT::OptimalityChecker::Vec x(2);
    x<< -1, 0;


    //-------------------------------------------------------------------------------//

    //check
    AOPT::OptimalityChecker oc;
    ASSERT_TRUE(oc.is_KKT_satisfied(&obj_func, ineq_cons, eq_cons, x, lambda, nu));

}

int main(int _argc, char** _argv){

    testing::InitGoogleTest(&_argc, _argv);
    return RUN_ALL_TESTS();

}

