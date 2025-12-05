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

using Vec = AOPT::NewtonMethods::Vec;
using Mat = AOPT::NewtonMethods::Mat;
using SMat = AOPT::NewtonMethods::SMat;

TEST(MassSpringSystemWithEqualityConstraints, ProjectOnAffine){

    const int n(4);
    
    Mat Ap(n,n);
    Ap <<  1, 2,  3,  4,
         -1, 0,  2,  3,
         -2, 1, -2,  4,
          0, 2,  3, -1;

    auto A = 0.5*(Ap.transpose() + Ap).sparseView();


    Vec b(n);
    b << -1, 0, 1, 2;
    
    Vec x(n);
    x <<0,0,0,0;
    
    ASSERT_GT((A*x - b).norm(), 1e-7);

    NewtonMethods::project_on_affine(x, A, b);

    ASSERT_LT((A*x - b).norm(), 1e-7);
    
}


TEST(MassSpringSystemWithEqualityConstraints, CheckMinimum){
    

    int n_grid_x(5), n_grid_y(7), func_index(0);

    //initial energy
    AOPT::MassSpringSystemT<AOPT::MassSpringProblem2DSparse> mss(n_grid_x, n_grid_y, func_index);

    SMat A;
    Vec b;
    mss.setup_linear_equality_constraints(A, b);

    //statistic instance
    auto opt_st = std::make_unique<AOPT::OptimizationStatistic>(mss.get_problem().get());

    //generate the start points
    Vec start_pts(opt_st->n_unknowns());
    for(int i(0); i<start_pts.size();i++){
        start_pts[i] = i;
    }

    //set points
    mss.set_spring_graph_points(start_pts);
    //initial energy
    auto energy = mss.initial_system_energy();
    std::cout<<"Initial MassSpring system energy is "<<energy<<std::endl;

    int max_iter(10000);
    AOPT::NewtonMethods::Vec x = AOPT::NewtonMethods::solve_equality_constrained(opt_st.get(), start_pts, A, b, 1e-7, max_iter);
    
    

    auto final_energy = mss.get_problem()->eval_f(x);
    
    std::cout<<" final energy = "<<std::setprecision(10)<<final_energy<<std::endl;
    
    double expected_energy(234.10279931131697);
    
    
    ASSERT_NEAR(final_energy, expected_energy, 1e-7);


}




int main(int _argc, char** _argv){

    testing::InitGoogleTest(&_argc, _argv);
    return RUN_ALL_TESTS();

}
