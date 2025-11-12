#include <iostream>
#include <Utils/StopWatch.hh>
#include <MassSpringSystemT.hh>
#include <Functions/ConstrainedSpringElement2D.hh>
#include <Functions/FunctionQuadratic2D.hh>
#include <Functions/FunctionNonConvex2D.hh>
#include <FunctionBase/ParametricFunctionWrapper.hh>
#include <Utils/RandomNumberGenerator.hh>
#include <Utils/DerivativeChecker.hh>
#include <Algorithms/GradientDescent.hh>

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


/** Simple log function used by one of the tests */
class LogFunction final : public FunctionBase {
public:
    // f(x) = log(x)

    // constructor
    LogFunction(){}

    // number of unknowns
    inline virtual int n_unknowns() { return 1; }

    // funcion evaluation
    inline virtual double eval_f(const Vec &_x) {
        return log(_x[0]);
    }

    // gradient evaluation
    inline virtual void eval_gradient(const Vec &_x, Vec &_g) {
        _g[0] = 1./_x[0];
    }

    // hessian matrix evaluation
    inline virtual void eval_hessian(const Vec &_x, Mat &_H) {}
};



/** checks that the basic functions of the ConstrainedSpringElement
 * give the expected results (computed by hand) */
TEST(ConstrainedSpringElement, CheckFunctions){
    typedef ConstrainedSpringElement2D::Vec Vec;
    typedef ConstrainedSpringElement2D::Mat Mat;

    ConstrainedSpringElement2D cse;


    Vec x(2), coeffs(3), g(2);
    Mat H(2,2);

    x << 1,1;
    coeffs << 1000000, 0, 0;

    ASSERT_EQ(cse.eval_f(x, coeffs), coeffs[0]);

    cse.eval_gradient(x,coeffs, g);
    Vec expected_g(2);
    expected_g << coeffs[0], coeffs[0];
    ASSERT_EQ(g, expected_g);

    cse.eval_hessian(x, coeffs, H);
    Mat expected_hess(2,2);
    expected_hess<<coeffs[0],  0,
            0,  coeffs[0];
    ASSERT_EQ(H, expected_hess);
}



/** This tests compares your implementation's energy computation result
 * with the solution's*/
TEST(MassSpringProblem, CheckFunctionForScenario1){
    MassSpringSystemT<AOPT::MassSpringProblem2DSparse> mss(5, 5, 1);
    //attach spring graph nodes to certain positions
    //it checks the first scenario which has constrained spring
    //elements at four corner nodes
    mss.add_constrained_spring_elements(1);

    //generate points
    const int n_vertices = 36;

    FunctionBase::Vec points(2*n_vertices);
    for(int i=0; i<n_vertices; ++i) {
        points[2*i] = sin(i * 0.3);
        points[2*i+1] = sin(i * 0.1);
    }

    mss.set_spring_graph_points(points);

    auto energy = mss.initial_system_energy();
    ASSERT_FLOAT_EQ(energy, 20207500);
}

TEST(MassSpringProblem, CheckFunctionForScenario2){
    MassSpringSystemT<AOPT::MassSpringProblem2DSparse> mss(5, 5, 1);
    //attach spring graph nodes to certain positions
    //it checks the first scenario which has constrained spring
    //elements at four corner nodes
    mss.add_constrained_spring_elements(2);

    //generate points
    const int n_vertices = 36;

    FunctionBase::Vec points(2*n_vertices);
    for(int i=0; i<n_vertices; ++i) {
        points[2*i] = sin(i * 0.3);
        points[2*i+1] = sin(i * 0.1);
    }

    mss.set_spring_graph_points(points);

    auto energy = mss.initial_system_energy();
    ASSERT_FLOAT_EQ(energy, 36036676);
}



/** Checks that the MSP2DSparse's energy gradient is close enough to the
 * finite difference computation based on the base function.
 * i.e. checks that the analytical expression of the gradient gives the same
 * result as the numerical approach */
TEST(MassSpringProblem, CheckGradient){
    MassSpringSystemT<AOPT::MassSpringProblem2DSparse> mss(5, 5, 1);
    //attach spring graph nodes to certain positions
    mss.add_constrained_spring_elements();

    DerivativeChecker checker;
    ASSERT_TRUE(checker.check_d1(*(mss.get_problem())));
}



/** Checks that your implementation of the back-tracking line search algorithm
works properly for a simple problem*/
TEST(LineSearch, CheckBackTrackingLineSearch){
    using Vec = FunctionQuadratic2D::Vec;

    LogFunction func;

    Vec start_pt(1);
    start_pt << 10;

    Vec g(1);
    func.eval_gradient(start_pt, g);

    double result = LineSearch::backtracking_line_search(&func, start_pt, g, -g, 200.);

    double expected_result = 84.375;
    ASSERT_EQ(result, expected_result);
}


/** Checks that the gradient descent gives the proper result */
TEST(GradientDescent, CheckAlgorithmOnSpringElementWithoutLength){

    using Vec = SpringElement2D::Vec;

    SpringElement2D sel;
    Vec coeffs(1);
    coeffs<<1;

    ParametricFunctionWrapper<SpringElement2D> non_param_sel(sel, coeffs);


    Vec start_pt(4);
    start_pt << 0, 0, 1, 1;

    Vec result = GradientDescent::solve(&non_param_sel, start_pt);

    Vec expected_result(4);
    expected_result << 0.5, 0.5, 0.5, 0.5;
    ASSERT_NEAR((result-expected_result).norm(),0, 1e-4);
}


/** Checks that the gradient descent gives the proper result */
TEST(GradientDescent, CheckAlgorithmOnSpringElementWithLength){

    using Vec = SpringElement2DWithLength::Vec;

    SpringElement2DWithLength selw;
    Vec coeffs(2);
    coeffs<<1, 5;

    ParametricFunctionWrapper<SpringElement2DWithLength> non_param_selw(selw, coeffs);


    Vec start_pt(4);
    start_pt << 10, 0, -10, 0;

    Vec result = GradientDescent::solve(&non_param_selw, start_pt);

    //The energy should be 0 at rest length
    ASSERT_NEAR(non_param_selw.eval_f(result), 0, 1e-10);
}


/** Checks that the gradient descent gives the proper result */
TEST(GradientDescent, CheckAlgorithmOnConstrainedSpringElement){

    using Vec = ConstrainedSpringElement2D::Vec;

    ConstrainedSpringElement2D csel;
    Vec coeffs(3);
    Eigen::VectorXd desired_location(2);
    desired_location<<-1, 1;
    coeffs<<10, desired_location[0], desired_location[1];

    ParametricFunctionWrapper<ConstrainedSpringElement2D> non_param_csel(csel, coeffs);


    Vec start_pt(2);
    start_pt << 10, 10;

    Vec result = GradientDescent::solve(&non_param_csel, start_pt);

    ASSERT_NEAR((result-desired_location).norm(), 0, 1e-4);
}

int main(int _argc, char** _argv){

    testing::InitGoogleTest(&_argc, _argv);
    return RUN_ALL_TESTS();

}
