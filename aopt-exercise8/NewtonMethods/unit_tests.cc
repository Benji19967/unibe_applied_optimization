#include <iostream>
#include <Functions/SpringElement2DWithLengthPSDHess.hh>
#include <Functions/SpringElement2D.hh>
#include <Functions/SpringElement2DWithLength.hh>
#include <Functions/ConstrainedSpringElement2D.hh>
#include <FunctionBase/DenseFunctionWrapper.hh>

#include <Algorithms/NewtonMethods.hh>

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


/** quadratic 2d function class with sparse hessian matrix */
class FunctionQuadraticNDSparse final : public FunctionBaseSparse {
public:
    // constructor
    FunctionQuadraticNDSparse(const SMat& _A, const Vec& _b, const double _c)
            : FunctionBaseSparse(), A_(_A), b_(_b), c_(_c) {
        if(_A.rows() != _A.cols())
            std::cerr << "Warning: matrix not square in FunctionQuadraticND" << std::endl;
        n_ = A_.rows();

        A_ = 0.5*(_A + SMat(_A.transpose()));
    }

    // number of unknowns
    inline virtual int n_unknowns() { return n_; }

    inline virtual double eval_f(const Vec &_x) {
        double v = 0.5 * _x.transpose() * A_ * _x;
        v += b_.transpose() * _x + c_;
        return v;
    }

    inline virtual void eval_gradient(const Vec &_x, Vec &_g) {
        _g = A_*_x + b_;
    }

    inline virtual void eval_hessian(const Vec &_x, SMat &_H) {
        _H = A_;
    }

private:
    // quadratic function 1/2 x^T A x + b^T x + c
    int n_;
    SMat A_;
    Vec b_;
    double c_;
};


/** non convex 2d function class with sparse hessian matrix */
class FunctionNonConvex2DSparse final : public FunctionBaseSparse {
public:
    // f(x,y) = (y-x^2)^2+cos^2(4*y)*(1-x)^2+x^2+y^2

    // constructor
    FunctionNonConvex2DSparse() {}

    // number of unknowns
    inline virtual int n_unknowns() { return 2; }

    inline virtual double eval_f(const Vec &_x) {
        double helper_0 = std::pow(_x[0], 2);

        return helper_0 + std::pow(_x[1], 2) + std::pow(1 - _x[0], 2)*std::pow(std::cos(4*_x[1]), 2) + std::pow(-helper_0 + _x[1], 2);
    }

    inline virtual void eval_gradient(const Vec &_x, Vec &_g) {
        double helper_0 = 2*_x[0];
        double helper_1 = std::pow(_x[0], 2);
        double helper_2 = 4*_x[1];
        double helper_3 = std::cos(helper_2);

        _g[0] = helper_0 + std::pow(helper_3, 2)*(helper_0 - 2) - 4*_x[0]*(-helper_1 + _x[1]);
        _g[1] = -2*helper_1 + helper_2 - 8*helper_3*std::pow(1 - _x[0], 2)*std::sin(helper_2);
    }

    inline virtual void eval_hessian(const Vec &_x, SMat &_H) {
        double helper_0 = 4*_x[1];
        double helper_1 = std::cos(helper_0);
        double helper_2 = std::pow(helper_1, 2);
        double helper_3 = -4*_x[0];
        double helper_4 = std::sin(helper_0);
        double helper_5 = helper_1*helper_4;
        double helper_6 = 32*std::pow(1 - _x[0], 2);

        _H.coeffRef(0,0) = -helper_0 + 2*helper_2 + 12*std::pow(_x[0], 2) + 2;
        _H.coeffRef(0,1) = helper_3 - 8*helper_5*(2*_x[0] - 2);
        _H.coeffRef(1,0) = helper_3 + helper_5*(16 - 16*_x[0]);
        _H.coeffRef(1,1) = -helper_2*helper_6 + std::pow(helper_4, 2)*helper_6 + 4;
    }
};


/** checks that the hessian of SpringElement2DWithLengthPSDHess
 * is positive definite */
TEST(SpringElement2DWithLengthPSDHess, CheckHessian){
    typedef SpringElement2DWithLengthPSDHess::Vec Vec;
    typedef SpringElement2DWithLengthPSDHess::Mat Mat;

    SpringElement2DWithLengthPSDHess sewlph;

    //initialize x and coeffs to test
    Vec x(4), coeffs(2), g(4);
    Mat H(4,4);

    x << 11.9968, 4.07199, 5.96417, 0.0146203;
    coeffs << 1, 1;

    sewlph.eval_hessian(x, coeffs, H);

    ASSERT_TRUE(H.llt().info() == Eigen::Success);
}


/** Checks that the standard newton method gives the proper result */
TEST(StandardNewton, CheckAlgorithm){
    using Vec = FunctionQuadraticNDSparse::Vec;
    using SMat = FunctionQuadraticNDSparse::SMat;

    SMat A(2,2);
    A.coeffRef(0,0) = 1; A.coeffRef(0,1) = 0;
    A.coeffRef(1,0) = 0; A.coeffRef(1,1) = 100;

    Vec b(2);
    b.setZero();

    FunctionQuadraticNDSparse func(A, b, 0);

    Vec start_pt(2);
    start_pt << 20, 20;

    Vec result = NewtonMethods::solve(&func, start_pt);
    std::cerr<<"x: "<<result.transpose()<<std::endl;

    Vec expected_result(2);
    expected_result << 0, 0;
    ASSERT_EQ(result, expected_result);
}


/** Checks that the gradient descent gives the proper result */
TEST(ProjectedNewton, CheckAlgorithmOnSpringElementWithoutLength){

    using Vec = SpringElement2D::Vec;

    SpringElement2D sel;
    Vec coeffs(1);
    coeffs<<1;

    DenseFunctionWrapper<SpringElement2D> sparse_sel(sel, coeffs);

    Vec start_pt(4);
    start_pt << 0, 0, 1, 1;

    Vec result = NewtonMethods::solve_with_projected_hessian(&sparse_sel, start_pt);

    Vec expected_result(4);
    expected_result << 0.5, 0.5, 0.5, 0.5;
    ASSERT_NEAR((result-expected_result).norm(),0, 1e-4);
}



/** Checks that the gradient descent gives the proper result */
TEST(ProjectedNewton, CheckAlgorithmOnSpringElementWithLength){

    using Vec = SpringElement2DWithLength::Vec;

    SpringElement2DWithLength selw;
    Vec coeffs(2);
    coeffs<<1, 5;

    DenseFunctionWrapper<SpringElement2DWithLength> non_param_selw(selw, coeffs);


    Vec start_pt(4);
    start_pt << 10, 0, -10, 0;

    Vec result = NewtonMethods::solve_with_projected_hessian(&non_param_selw, start_pt);

    //The energy should be 0 at rest length
    ASSERT_NEAR(non_param_selw.eval_f(result), 0, 1e-9);
}


/** Checks that the gradient descent gives the proper result */
TEST(ProjectedNewton, CheckAlgorithmOnConstrainedSpringElement){

    using Vec = ConstrainedSpringElement2D::Vec;

    ConstrainedSpringElement2D csel;
    Vec coeffs(3);
    Eigen::VectorXd desired_location(2);
    desired_location<<-1, 1;
    coeffs<<10, desired_location[0], desired_location[1];

    DenseFunctionWrapper<ConstrainedSpringElement2D> non_param_csel(csel, coeffs);


    Vec start_pt(2);
    start_pt << 10, 10;

    Vec result = NewtonMethods::solve_with_projected_hessian(&non_param_csel, start_pt);

    ASSERT_NEAR((result-desired_location).norm(), 0, 1e-4);
}



/** Checks that the projected newton method gives the proper result */
TEST(ProjectedNewton, CheckAlgorithmOnNonConvexFunction){
    using Vec = FunctionQuadraticNDSparse::Vec;

    //set start point
    Vec start_pt(2);
    start_pt << 2, 2;

    FunctionNonConvex2DSparse func;

    Vec result = NewtonMethods::solve_with_projected_hessian(&func, start_pt);
    
    Vec expected_result(2);
    expected_result<<0.118892,0.337706;

    ASSERT_NEAR((result-expected_result).norm(), 0, 1e-5);
}



int main(int _argc, char** _argv){

    testing::InitGoogleTest(&_argc, _argv);
    return RUN_ALL_TESTS();

}
