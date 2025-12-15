#include <iostream>

#include <Utils/StopWatch.hh>
#include <Utils/RandomNumberGenerator.hh>
#include <Utils/OptimizationStatistic.hh>

#include <Functions/FunctionQuadratic2D.hh>
#include <Functions/FunctionQuadraticND.hh>
#include <Functions/ConstrainedSpringElement2DLeastSquare.hh>
#include <Functions/MassSpringProblem2DLeastSquare.hh>
#include <Functions/MassSpringProblem2DSparse.hh>
#include <Functions/SpringElement2DLeastSquare.hh>
#include <Functions/SpringElement2DWithLengthLeastSquare.hh>
#include <Functions/ConstrainedSpringElement2DLeastSquare.hh>

#include <Algorithms/NewtonMethods.hh>

#include <MassSpringSystemT.hh>

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




TEST(MassSpringSystemWithInfeasibleStart, CheckMinimumWithSpringWithoutLength){

    int n_grid_x(10), n_grid_y(10), func_index(0);

    //initial energy
    AOPT::MassSpringSystemT<AOPT::MassSpringProblem2DSparse> mss(n_grid_x, n_grid_y, func_index);

    AOPT::NewtonMethods::SMat A;
    AOPT::NewtonMethods::Vec b;
    mss.setup_linear_equality_constraints(A, b);

    //statistic instance
    auto opt_st = std::make_unique<AOPT::OptimizationStatistic>(mss.get_problem().get());

    //generate the start points
    AOPT::RandomNumberGenerator rng(-10., 10.);
    auto start_pts = rng.get_random_nd_vector(opt_st->n_unknowns());

    //set points
    mss.set_spring_graph_points(start_pts);
    //initial energy
    auto energy = mss.initial_system_energy();
    std::cout<<"Initial MassSpring system energy is "<<energy<<std::endl;

    int max_iter(10000);
    AOPT::NewtonMethods::Vec x = AOPT::NewtonMethods::solve_equality_constrained_with_infeasible_start(opt_st.get(), start_pts, A, b, 1e-6, 1e-8, max_iter);

    double final_energy = mss.get_problem()->eval_f(x);
    double expected_final_energy(565.10064157713407);

    ASSERT_NEAR(final_energy, expected_final_energy, 1e-5);
}



TEST(MassSpringSystemWithInfeasibleStart, CheckMinimumWithSpringWithLength){

    int n_grid_x(40), n_grid_y(30), func_index(1);

    //initial energy
    AOPT::MassSpringSystemT<AOPT::MassSpringProblem2DSparse> mss(n_grid_x, n_grid_y, func_index);

    AOPT::NewtonMethods::SMat A;
    AOPT::NewtonMethods::Vec b;
    mss.setup_linear_equality_constraints(A, b);

    //statistic instance
    auto opt_st = std::make_unique<AOPT::OptimizationStatistic>(mss.get_problem().get());

    //generate the start points
    AOPT::RandomNumberGenerator rng(-10., 10.);
    auto start_pts = rng.get_random_nd_vector(opt_st->n_unknowns());

    //set points
    mss.set_spring_graph_points(start_pts);
    //initial energy
    auto energy = mss.initial_system_energy();
    std::cout<<"Initial MassSpring system energy is "<<energy<<std::endl;

    int max_iter(10000);
    AOPT::NewtonMethods::Vec x = AOPT::NewtonMethods::solve_equality_constrained_with_infeasible_start(opt_st.get(), start_pts, A, b, 1e-6, 1e-8, max_iter);

    double final_energy = mss.get_problem()->eval_f(x);
    double expected_final_energy(17389.532762494073);

    ASSERT_NEAR(final_energy, expected_final_energy, 1e-5);
}




//----------------------------------------- BONUS UNIT TESTS

TEST(MassSpringSystemWithHybrid, CheckMinimumWithSpringWithoutLength){

    int n_grid_x(25), n_grid_y(55), func_index(0);

    //initial energy
    AOPT::MassSpringSystemT<AOPT::MassSpringProblem2DSparse> mss(n_grid_x, n_grid_y, func_index);

    AOPT::NewtonMethods::SMat A;
    AOPT::NewtonMethods::Vec b;
    mss.setup_linear_equality_constraints(A, b);

    //statistic instance
    auto opt_st = std::make_unique<AOPT::OptimizationStatistic>(mss.get_problem().get());

    //generate the start points
    AOPT::RandomNumberGenerator rng(-10., 10.);
    auto start_pts = rng.get_random_nd_vector(opt_st->n_unknowns());

    //set points
    mss.set_spring_graph_points(start_pts);
    //initial energy
    auto energy = mss.initial_system_energy();
    std::cout<<"Initial MassSpring system energy is "<<energy<<std::endl;

    int max_iter(10000);
    AOPT::NewtonMethods::Vec x = AOPT::NewtonMethods::solve_equality_constrained_hybrid(opt_st.get(), start_pts, A, b, 1e-6, 1e-8, max_iter);

    double final_energy = mss.get_problem()->eval_f(x);
    double expected_final_energy(6127.1767891866275);

    ASSERT_NEAR(final_energy, expected_final_energy, 1e-5);

}


TEST(MassSpringSystemWithHybrid, CheckMinimumWithSpringWithLength){

    int n_grid_x(10), n_grid_y(5), func_index(1);

    //initial energy
    AOPT::MassSpringSystemT<AOPT::MassSpringProblem2DSparse> mss(n_grid_x, n_grid_y, func_index);

    AOPT::NewtonMethods::SMat A;
    AOPT::NewtonMethods::Vec b;
    mss.setup_linear_equality_constraints(A, b);

    //statistic instance
    auto opt_st = std::make_unique<AOPT::OptimizationStatistic>(mss.get_problem().get());

    //generate the start points
    AOPT::RandomNumberGenerator rng(-10., 10.);
    auto start_pts = rng.get_random_nd_vector(opt_st->n_unknowns());

    //set points
    mss.set_spring_graph_points(start_pts);
    //initial energy
    auto energy = mss.initial_system_energy();
    std::cout<<"Initial MassSpring system energy is "<<energy<<std::endl;

    int max_iter(10000);
    AOPT::NewtonMethods::Vec x = AOPT::NewtonMethods::solve_equality_constrained_hybrid(opt_st.get(), start_pts, A, b, 1e-6, 1e-8, max_iter);

    double final_energy = mss.get_problem()->eval_f(x);
    double expected_final_energy(868.2599748938944);

    ASSERT_NEAR(final_energy, expected_final_energy, 1e-5);

}



int main(int _argc, char** _argv){

    testing::InitGoogleTest(&_argc, _argv);
    return RUN_ALL_TESTS();

}
