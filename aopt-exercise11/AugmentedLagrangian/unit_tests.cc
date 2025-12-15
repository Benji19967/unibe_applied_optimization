
#include <iostream>
#include <Utils/StopWatch.hh>
#include <MassSpringSystemT.hh>
#include <Functions/FunctionQuadraticND.hh>
#include <Functions/CircleConstraint2D.hh>
#include <Functions/CircleConstraintSquared2D.hh>
#include <Functions/AugmentedLagrangianProblem.hh>
#include <Algorithms/AugmentedLagrangian.hh>

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

/** simple function for unit-testing */



TEST(CircleConstraint2D, CheckFunctions){
    
    CircleConstraint2D constraint(2,0,1.5,2.5,0.5);

    CircleConstraint2D::Vec x(2);
    x<<-1, 2;

    ASSERT_FLOAT_EQ(constraint.eval_f(x), 6.25);

    CircleConstraint2D::Vec g(2);
    constraint.eval_gradient(x, g);
    std::cout<<" g: "<<g.transpose()<<std::endl;

    CircleConstraint2D::Vec expected_g(2);
    expected_g << -5, -1;
    ASSERT_FLOAT_EQ((g - expected_g).squaredNorm(), 0.0);

    CircleConstraint2D::SMat H(2,2);
    constraint.eval_hessian(x, H);

    CircleConstraint2D::SMat expected_H(2,2);
    expected_H.insert(0, 0) = 2.;
    expected_H.insert(1, 1) = 2.;

    ASSERT_EQ((H - expected_H).norm(), 0.);
}





TEST(CircleConstraint2D, CheckFunctionsWithDifferentIndices){


    const int n(6);

    CircleConstraint2D::Vec x(n*2);

    for(int i(0); i<n; i++){
        x(2*i) = i;
        x(2*i+1) = i;
    }

    for(int i(0); i<n; i++){

        CircleConstraint2D constraint(n, i, i+1, i+1, i);

        ASSERT_EQ(constraint.eval_f(x), 2 - i*i);

        CircleConstraint2D::Vec g(n*2);
        constraint.eval_gradient(x, g);

        CircleConstraint2D::Vec expected_g(n*2);
        expected_g.setZero();
        expected_g(i*2) = -2;
        expected_g(i*2 + 1) = -2;

        ASSERT_EQ(g, expected_g);

        CircleConstraint2D::SMat H(n*2,n*2);
        constraint.eval_hessian(x, H);

        CircleConstraint2D::SMat expected_H(n*2,n*2);
        for(auto j(0); j<n; j++){
            if(i==j){
                expected_H.insert(2*i,2*i) = 2;
                expected_H.insert(2*i+1,2*i+1) = 2;
            }
        }
        ASSERT_EQ((H - expected_H).norm(), 0.);
    }
}





TEST(CircleConstraintSquared2D, CheckFunctions){

    const int n(10);
    const int index(2);

    CircleConstraintSquared2D constraint(n,index,1.5,2.5,0.5);

    CircleConstraintSquared2D::Vec x(n);
    x.setZero();
    x(2*index)     = -1;
    x(2*index + 1) =  2;

    ASSERT_FLOAT_EQ(constraint.eval_f(x), 39.0625);

    CircleConstraintSquared2D::Vec g(n);
    constraint.eval_gradient(x, g);
    std::cout<<" g: "<<std::endl<<g.transpose()<<std::endl;

    CircleConstraintSquared2D::Vec expected_g(n);
    expected_g.setZero();
    expected_g(2*index)   = -62.5;
    expected_g(2*index+1) = -12.5;
    ASSERT_FLOAT_EQ((g - expected_g).squaredNorm(), 0.0);

    CircleConstraintSquared2D::SMat H(n,n);
    constraint.eval_hessian(x, H);

    CircleConstraintSquared2D::SMat expected_H(n,n);
    expected_H.insert(2*index, 2*index) = 75.;
    expected_H.insert(2*index, 2*index+1) = 10.;
    expected_H.insert(2*index+1, 2*index) = 10.;
    expected_H.insert(2*index+1, 2*index+1) = 27.;
    ASSERT_EQ((H - expected_H).norm(), 0.);
}



TEST(CircleConstraintSquared2D, CheckFunctionsWithDifferentIndices){


    const int n(6);

    CircleConstraintSquared2D::Vec x(n*2);

    for(int i(0); i<n; i++){
        x(2*i) = i;
        x(2*i+1) = i;
    }

    for(int i(0); i<n; i++){

        CircleConstraintSquared2D constraint(n, i, i+1, i+1, i);

        double expected_f(2 - i*i);
        ASSERT_EQ(constraint.eval_f(x), expected_f * expected_f);

        CircleConstraintSquared2D::Vec g(n*2);
        constraint.eval_gradient(x, g);

        CircleConstraintSquared2D::Vec expected_g(n*2);
        expected_g.setZero();
        expected_g(i*2) = -4*(2 - i*i);
        expected_g(i*2 + 1) = -4*(2 - i*i);

        ASSERT_EQ(g, expected_g);
    }
}


TEST(AugmentedLagrangianProblem, CheckFunctionsWithBoundaryConstraintsA_lowDimProblem){

    using Vec  = AugmentedLagrangianProblem::Vec;
    using SMat = AugmentedLagrangianProblem::SMat;
    using Mat  = AugmentedLagrangianProblem::Mat;
    const int dim(1);

    //initial energy
    AOPT::MassSpringSystemT<AOPT::MassSpringProblem2DSparse> mss(dim, dim, 0);
    mss.add_boundary_constraints();

    const int n = mss.get_problem()->n_unknowns();

    auto opt_st = std::make_unique<AOPT::OptimizationStatistic>(mss.get_problem().get());

    //generate the start points
    Vec start_pts(n);
    for(int i(0); i<n/2; i++){
        start_pts(2*i)   = i % (dim+1);
        start_pts(2*i+1) = i % (dim+1) + i;
    }

    //set points
    mss.set_spring_graph_points(start_pts);


    Vec nu(mss.get_constraints().size());
    nu.setOnes() * 0.1;
    double mu(0.5);


    //initialize the augmented lagrangian problem for the unconstrained solver
    AugmentedLagrangianProblem problem(mss.get_problem().get(),
                                       mss.get_constraints(),
                                       mss.get_constraints_squared(),nu,mu);


    Vec x(n);
    x.setZero();
    ASSERT_NEAR(problem.eval_f(x), -1.75, 1e-6);
    ASSERT_NEAR(problem.eval_f(start_pts), 66.25, 1e-6);


    Vec g(n);
    problem.eval_gradient(start_pts, g);
    std::cout<<"g: "<<g.transpose()<<std::endl;
    Vec expected_g(n);
    expected_g << -2.75, -8.75,  3.75,  5.25, -3.75,  5.25,  8.75, 55.25;
    ASSERT_NEAR((g-expected_g).squaredNorm(), 0.0, 1e-6);


    SMat H(n,n);
    problem.eval_hessian(x, H);
    std::cout<<"H: "<<H<<std::endl;
    Mat expected_H(n,n);
    expected_H<<5, 0.5, -1, 0, -1, 0, -1, 0,
            0.5, 5, 0, -1, 0, -1, 0, -1,
            -1, 0, 5, 0.5, -1, 0, -1, 0,
            0, -1, 0.5, 5, 0, -1, 0, -1,
            -1, 0, -1, 0, 5, 0.5, -1, 0,
            0, -1, 0, -1, 0.5, 5, 0, -1,
            -1, 0, -1, 0, -1, 0, 5, 0.5,
            0, -1, 0, -1, 0, -1, 0.5, 5;
    ASSERT_NEAR((H - expected_H.sparseView()).squaredNorm(), 0.0, 1e-6);


}




TEST(AugmentedLagrangianProblem, CheckFunctionsWithBoundaryConstraintsA_higherDimProblem){

    using Vec  = AugmentedLagrangianProblem::Vec;
    using SMat = AugmentedLagrangianProblem::SMat;
    const int dim(10);

    //initial energy
    AOPT::MassSpringSystemT<AOPT::MassSpringProblem2DSparse> mss(dim, dim, 0);
    mss.add_boundary_constraints();

    const int n = mss.get_problem()->n_unknowns();

    std::cout<<" n = "<<n<<std::endl;

    auto opt_st = std::make_unique<AOPT::OptimizationStatistic>(mss.get_problem().get());

    //generate the start points
    Vec start_pts(n);
    for(int i(0); i<n/2; i++){
        start_pts(2*i)   = i % (dim+1);
        start_pts(2*i+1) = i % (dim+1) + (i+1)%(dim);
    }

    //set points
    mss.set_spring_graph_points(start_pts);


    Vec nu(mss.get_constraints().size());
    nu.setOnes() * 0.1;
    double mu(0.5);


    //initialize the augmented lagrangian problem for the unconstrained solver
    AugmentedLagrangianProblem problem(mss.get_problem().get(),
                                       mss.get_constraints(),
                                       mss.get_constraints_squared(),nu,mu);


    Vec x(n);
    x.setZero();
    ASSERT_NEAR(problem.eval_f(x), 23000, 1e-6);
    ASSERT_NEAR(problem.eval_f(start_pts), 43479, 1e-6);

    Vec g(n);
    problem.eval_gradient(start_pts, g);
    const double expected_g_norm(10747932);
    ASSERT_NEAR(g.squaredNorm(), expected_g_norm, 1e-6);


    SMat H(n,n);
    problem.eval_hessian(x, H);
    const double expected_H_norm(215776);
    ASSERT_NEAR(H.squaredNorm(), expected_H_norm, 1e-6);

}



int main(int _argc, char** _argv){

    testing::InitGoogleTest(&_argc, _argv);
    return RUN_ALL_TESTS();

}
