
#include <iostream>
#include <Utils/StopWatch.hh>
#include <MassSpringSystemT.hh>
#include <Algorithms/InteriorPoint.hh>
#include <Functions/AreaConstraint2D.hh>

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


class FunctionQuadratic2DSparse final : public FunctionBaseSparse {
public:

    // constructor
    FunctionQuadratic2DSparse(int n)
            : FunctionBaseSparse(), n_(n){}

    // number of unknowns
    inline virtual int n_unknowns() { return n_; }

    inline virtual double eval_f(const Vec &_x) {
        return _x.transpose() * _x;
    }

    inline virtual void eval_gradient(const Vec &_x, Vec &_g) {
        _g = _x;
    }

    inline virtual void eval_hessian(const Vec &_x, SMat &_H) {
        _H.reserve(n_ * n_);
        for(int i(0); i<n_; i++){
            _H.insert(i,i) = 1;
        }
    }

private:

    int n_;
};



TEST(AreaConstraint2D, CheckFunctions){

    typedef AreaConstraint2D::Vec Vec;
    const int n(6);

    AreaConstraint2D area_constraint(n, 0, 1, 2);


    Vec x(n);
    x << 0, 0, 1, 0, 0, 1;

    ASSERT_FLOAT_EQ(area_constraint.eval_f(x), -0.5);

    AreaConstraint2D::Vec g(n);
    area_constraint.eval_gradient(x, g);
    std::cout<<" g: "<<g.transpose()<<std::endl;

    AreaConstraint2D::Vec expected_g(n);
    expected_g << 0.5,  0.5, -0.5,  0,  0, -0.5;
    ASSERT_FLOAT_EQ((g - expected_g).squaredNorm(), 0.0);

    AreaConstraint2D::SMat H(n,n);
    area_constraint.eval_hessian(x, H);
    //std::cout<<" H: "<<std::endl<<H<<std::endl;

    AreaConstraint2D::Mat expected_H(n,n);
    expected_H <<
                  0, 0, 0, -0.5, 0, 0.5,
            0, 0, 0.5, 0, -0.5, 0,
            0, 0.5, 0, 0, 0, -0.5,
            -0.5, 0, 0, 0, 0.5, 0,
            0, -0.5, 0, 0.5, 0, 0,
            0.5, 0, -0.5, 0, 0, 0;

    ASSERT_EQ((H - expected_H.sparseView()).norm(), 0.);

}




TEST(InteriorPointProblem, CheckFunctions){

    typedef InteriorPointProblem::Vec Vec;

    const int n(8);

    std::vector<FunctionBaseSparse*> constraints;


    constraints.push_back(new AreaConstraint2D(n, 0, 1, 3));
    constraints.push_back(new AreaConstraint2D(n, 2, 3, 1));
    constraints.push_back(new AreaConstraint2D(n, 3, 0, 2));
    constraints.push_back(new AreaConstraint2D(n, 0, 1, 2));


    Vec x(n);
    x << 0, 0, 1, 0, 1, 1, 0, 1;


    FunctionQuadratic2DSparse obj(n);


    InteriorPointProblem interior_point_problem(&obj, constraints);
    interior_point_problem.t() = 0.1;



    EXPECT_FLOAT_EQ(interior_point_problem.eval_f(x), 31.725887);

    AreaConstraint2D::Vec g(n);
    interior_point_problem.eval_gradient(x, g);
    std::cout<<" g: "<<g.transpose()<<std::endl;

    AreaConstraint2D::Vec expected_g(n);
    expected_g << 20,  20, -19,  20, -19, -19,  20, -19;

    EXPECT_NEAR((g - expected_g).squaredNorm(), 0.0, 1e-9);

    AreaConstraint2D::SMat H(n,n);
    interior_point_problem.eval_hessian(x, H);
    //std::cout<<" H: "<<std::endl<<H<<std::endl;

    AreaConstraint2D::Mat expected_H(n,n);
    expected_H <<  21,  10, -20, -10,   0, -10,   0,  10,
                   10,  21,  10,   0, -10,   0, -10, -20,
                  -20,  10,  21, -10,   0, -10,   0,  10,
                  -10,   0, -10,  21,  10, -20,  10,   0,
                    0, -10,   0,  10,  21,  10, -20, -10,
                  -10,   0, -10, -20,  10,  21,  10,   0,
                    0, -10,   0,  10, -20,  10,  21, -10,
                   10, -20,  10,   0, -10,   0, -10,  21;


    EXPECT_NEAR((H - expected_H.sparseView()).norm(), 0., 1e-7);


    for(int i(0); i<(int)constraints.size(); i++){
        delete constraints[i];
    }
}




TEST(InteriorPointMethod, CheckMinimum){

    const int dim(3);
    typedef MassSpringProblem2DSparse::Vec Vec;

    AOPT::MassSpringSystemT<AOPT::MassSpringProblem2DSparse> mss(dim, dim, 1);

    //fix four corner nodes and one interior node to certain positions
    mss.add_constrained_spring_elements();
    mss.add_constrained_spring_element_for_center_spring_node();
    mss.add_area_constraints();


    const int n(mss.get_problem()->n_unknowns());


    //generate the start points
    /*Vec start_pts(n);
    for(int i(0); i<n; i++){
        start_pts(i) = (i + 1) % ((dim + i) % (dim+1) + 1);
    }
    mss.set_spring_graph_points(start_pts);*/

    auto start_pts = mss.get_spring_graph_points();
    mss.set_spring_graph_points(start_pts);

    const int max_iter(1000);
    //solve
    AOPT::InteriorPoint::Vec x = AOPT::InteriorPoint::solve(mss.get_problem().get(), start_pts, mss.get_constraints(), 1e-4, 10, max_iter);

    //std::cout<<"x = "<<x.transpose()<<std::endl;

    double expected_final_energy(3253.8495082045383);

    double final_energy = mss.get_problem().get()->eval_f(x);


    ASSERT_NEAR(final_energy, expected_final_energy, 1e-6);

}


int main(int _argc, char** _argv){

    testing::InitGoogleTest(&_argc, _argv);
    return RUN_ALL_TESTS();

}
