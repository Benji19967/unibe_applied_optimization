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




TEST(SpringElement2DWithoutLengthLeastSquare, FunctionsEvaluation){

    SpringElement2DLeastSquare spring;

    typedef SpringElement2DLeastSquare::Vec Vec;


    Vec x(2), coeffs(1);
    x<<-0.7, 2;
    coeffs<<3;

    Vec g(2);

    ASSERT_NEAR(spring.eval_f(x, coeffs), -4.67653718, 1e-6);



    spring.eval_gradient(x, coeffs, g);
    //std::cout<<"energy gradient: "<<std::setprecision(10)<<g<<std::endl;

    Vec expected_g(2);
    expected_g <<  1.732050808,  -1.732050808;
    ASSERT_NEAR((g - expected_g).norm(), 0.0, 1e-6);

}



TEST(ConstrainedSpringElement2DLeastSquare, FunctionsEvaluation){

    ConstrainedSpringElement2DLeastSquare spring;

    typedef ConstrainedSpringElement2DLeastSquare::Vec Vec;

    Vec x(1), coeffs(2);
    x<<0.5;
    coeffs<<1000000, 0.1;

    Vec g(1);

    ASSERT_EQ(spring.eval_f(x, coeffs), 400);

    spring.eval_gradient(x, coeffs, g);
    //std::cout<<"energy gradient: "<<std::setprecision(10)<<g<<std::endl;

    Vec expected_g(1);
    expected_g << 1000;

    ASSERT_EQ(g, expected_g);
}



TEST(SpringElement2DWithLengthLeastSquare, FunctionsEvaluation){

    SpringElement2DWithLengthLeastSquare spring;

    typedef SpringElement2DWithLengthLeastSquare::Vec Vec;


    Vec x_a(2), x_b(2), coeffs(2);

    x_a << -1, 2;
    x_b << 3,-4;
    coeffs << 2, 5;

    Vec x(4);
    x<<x_a, x_b;

    Vec g(4);

    ASSERT_NEAR(spring.eval_f(x, coeffs), 38.183766184073569, 1e-6);

    spring.eval_gradient(x, coeffs, g);
    //std::cout<<"energy gradient: "<<std::setprecision(10)<<g<<std::endl;

    Vec expected_g(4);
    expected_g << -11.3137085,  16.97056275,   11.3137085, -16.97056275;
    ASSERT_NEAR((g - expected_g).norm(), 0.0, 1e-6);

}


TEST(SpringElement2DWithLengthLeastSquare, EnergyOfSpringAtRestIsAlwaysZero){

    SpringElement2DWithLengthLeastSquare spring;

    typedef SpringElement2DWithLengthLeastSquare::Vec Vec;


    RandomNumberGenerator rng;

    const int iteration_count(100);

    for(int i(0); i< iteration_count; i++){
        Vec x_a(rng.get_random_nd_vector(2)),
                direction(rng.get_random_nd_vector(2)),
                coeffs(rng.get_random_nd_vector(2));

        //make sure k and l are positive
        coeffs[0] = std::abs(coeffs[0]);
        coeffs[1] = std::abs(coeffs[1]);

        //generate a second point at length (coeffs[1]) distance from x
        direction.normalize();
        Vec x_b = x_a + coeffs[1] * direction;

        Vec x(4);
        x<<x_a, x_b;

        //evaluation of a spring at rest is always zero
        ASSERT_NEAR(spring.eval_f(x, coeffs), 0.0, 1e-7);
    }
}






TEST(MassSpringProblemWithoutLengthLeastSquares, InitialEnergy){

    int grid_x(4), grid_y(8), spring_type(0);

    //construct mass spring system
    AOPT::MassSpringSystemT<AOPT::MassSpringProblem2DLeastSquare> mss(grid_x, grid_y, spring_type, true);
    //attach spring graph nodes to certain positions
    mss.add_constrained_spring_elements();


    //generate the start points
    FunctionBase::Vec points(2*mss.n_grid_points());
    for(int i(0); i<points.size()/2; i++){
        points[2*i]   = (i % grid_x);
        points[2*i+1] = (i % grid_y);
    }

    //set points
    mss.set_spring_graph_points(points);

    //initial energy
    auto energy = mss.initial_system_energy();

    ASSERT_FLOAT_EQ(energy, 27201040);
}





TEST(MassSpringProblemWithLengthLeastSquares, InitialEnergy){

    int grid_x(6), grid_y(7), spring_type(1);

    //construct mass spring system
    AOPT::MassSpringSystemT<AOPT::MassSpringProblem2DLeastSquare> mss(grid_x, grid_y, spring_type, true);
    //attach spring graph nodes to certain positions
    mss.add_constrained_spring_elements();


    //generate the start points
    FunctionBase::Vec points(2*mss.n_grid_points());
    for(int i(0); i<points.size()/2; i++){
        points[2*i]   = (i % grid_x);
        points[2*i+1] = (i % grid_y);
    }

    //set points
    mss.set_spring_graph_points(points);

    //initial energy
    auto energy = mss.initial_system_energy();

    ASSERT_FLOAT_EQ(energy, 28106546);
}


TEST(MassSpringProblemWithLengthLeastSquares, convergesToSolution){

    int grid_x(5), grid_y(3), spring_type(1);

    //construct mass spring system
    AOPT::MassSpringSystemT<AOPT::MassSpringProblem2DLeastSquare> mss(grid_x, grid_y, spring_type, true);
    //attach spring graph nodes to certain positions
    mss.add_constrained_spring_elements();

    auto opt_st = std::make_unique<AOPT::OptimizationStatistic>(mss.get_problem().get());

    //generate the start points
    FunctionBase::Vec points(2*mss.n_grid_points());
    for(int i(0); i<points.size()/2; i++){
        points[2*i]   = (i % (grid_x + 2));
        points[2*i+1] = (i % (grid_y + 3));
    }

    //set points
    mss.set_spring_graph_points(points);
    
    auto init_energy = mss.initial_system_energy();

    auto x_min = AOPT::NewtonMethods::solve(opt_st.get(), points, 1e-5);
    opt_st->print_statistics();

    const double expected_final_energy(324.44395657791779);
    const double final_energy(opt_st.get()->eval_f(x_min));

    ASSERT_NEAR(final_energy, expected_final_energy, 1e-7);


}




int main(int _argc, char** _argv){

    testing::InitGoogleTest(&_argc, _argv);
    return RUN_ALL_TESTS();

}
