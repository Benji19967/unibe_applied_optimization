#include <iostream>
#include <Utils/StopWatch.hh>
#include <MassSpringSystemT.hh>
#include <Utils/RandomNumberGenerator.hh>
#include <Utils/DerivativeChecker.hh>

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


/** Checks that your 5x7 grid MSP has the right number of grid points and edges*/
TEST(MassSpringSystem, ProblemSetup){
    MassSpringSystemT<AOPT::MassSpringProblem2DDense> system(5, 7);
    ASSERT_EQ(system.n_grid_points(), 48);
    ASSERT_EQ(system.n_edges(), 152);
}


/** Compares your Spring Element without length's functions (energy, gradient and Hessian)
 * results with the manually computed ones */
TEST(SpringElements, SpringElement2DFunctions){
    typedef SpringElement2D::Vec Vec;
    typedef SpringElement2D::Mat Mat;

    SpringElement2D sewl;


    Vec x(4), coeffs(2), g(4);
    Mat H(4,4);

    x << 0,0,1,1;
    coeffs << 1,1;

    ASSERT_EQ(sewl.eval_f(x, coeffs), 1);

    sewl.eval_gradient(x,coeffs, g);
    Vec expected_g(4);
    expected_g << -1, -1, 1, 1;
    ASSERT_EQ(g, expected_g);

    sewl.eval_hessian(x, coeffs, H);
    Mat expected_hess(4,4);
    expected_hess<<1,  0, -1,  0,
            0,  1,  0, -1,
            -1,  0,  1,  0,
            0, -1,  0,  1;
    ASSERT_EQ(H, expected_hess);
}



/** Compares your Spring Element with length's functions (energy, gradient and Hessian)
 * results with the manually computed ones */
TEST(SpringElements, SpringElement2DWithLengthFunctions){
    typedef SpringElement2DWithLength::Vec Vec;
    typedef SpringElement2DWithLength::Mat Mat;

    SpringElement2DWithLength sewl;


    Vec x(4), coeffs(2), g(4);
    Mat H(4,4);

    x << 0,0,1,1;
    coeffs << 1,1;

    ASSERT_EQ(sewl.eval_f(x, coeffs), 0.5);

    sewl.eval_gradient(x,coeffs, g);
    Vec expected_g(4);
    expected_g << -2, -2, 2, 2;
    ASSERT_EQ(g, expected_g);

    sewl.eval_hessian(x, coeffs, H);
    Mat expected_hess(4,4);
    expected_hess<<6,  4, -6, -4,
            4,  6, -4, -6,
            -6, -4,  6,  4,
            -4, -6,  4,  6;
    ASSERT_EQ(H, expected_hess);
}



/** Compares your MSP Dense's functions (energy, gradient and Hessian)
 * results with the manually computed ones */
TEST(MassSpringProblem, MassSpringProblem2DDenseFunctions){

    typedef MassSpringProblem2DDense::Vec Vec;
    typedef MassSpringProblem2DDense::Mat Mat;

    int n_unknowns(8);

    SpringElement2DWithLength sewl;
    MassSpringProblem2DDense msp(sewl, n_unknowns);

    Vec x(2 * n_unknowns), g;
    Mat H;

    for(int i(0); i<n_unknowns; i++){
        x[2*i] = i;
        x[2*i+1] = -i+1;
    }

    for(int i(0); i<n_unknowns/2; i++){
        msp.add_spring_element(i, (i+1)%(n_unknowns/2));
    }

    ASSERT_EQ(msp.eval_f(x), 146);


    msp.eval_gradient(x, g);
    Vec expected_g(n_unknowns);
    expected_g << -104,104,0,0,0,0,104,-104;
    ASSERT_EQ(g, expected_g);

    msp.eval_hessian(x, H);
    Mat expected_hess(n_unknowns, n_unknowns);
    expected_hess<<76, -40,  -6,   4,   0,   0, -70,  36,
            -40,  76,   4,  -6,   0,   0,  36, -70,
            -6,   4,  12,  -8,  -6,   4,   0,   0,
            4,  -6,  -8,  12,   4,  -6,   0,   0,
            0,   0,  -6,   4,  12,  -8,  -6,   4,
            0,   0,   4,  -6,  -8,  12,   4,  -6,
            -70,  36,   0,   0,  -6,   4,  76, -40,
            36, -70,   0,   0,   4,  -6, -40,  76;
    ASSERT_EQ(H, expected_hess);
}


/** Compares your MSP Sparse's functions (energy, gradient and Hessian)
 * results with the manually computed ones */
TEST(MassSpringProblem, MassSpringProblem2DSparseFunctions){

    typedef MassSpringProblem2DSparse::Vec Vec;
    typedef MassSpringProblem2DSparse::Mat Mat;
    typedef MassSpringProblem2DSparse::SMat SMat;

    int n_unknowns(8);

    SpringElement2DWithLength sewl;
    MassSpringProblem2DSparse msp(sewl, n_unknowns);


    Vec x(2 * n_unknowns), g;
    SMat H;

    for(int i(0); i<n_unknowns; i++){
        x[2*i] = i;
        x[2*i+1] = -i+1;
    }

    for(int i(0); i<n_unknowns/2; i++){
        msp.add_spring_element(i, (i+1)%(n_unknowns/2));
    }


    ASSERT_EQ(msp.eval_f(x), 146);

    msp.eval_gradient(x, g);
    Vec expected_g(n_unknowns);
    expected_g << -104,104,0,0,0,0,104,-104;
    ASSERT_EQ(g, expected_g);

    msp.eval_hessian(x, H);
    Mat expected_hess(n_unknowns, n_unknowns);
    expected_hess<<76, -40,  -6,   4,   0,   0, -70,  36,
            -40,  76,   4,  -6,   0,   0,  36, -70,
            -6,   4,  12,  -8,  -6,   4,   0,   0,
            4,  -6,  -8,  12,   4,  -6,   0,   0,
            0,   0,  -6,   4,  12,  -8,  -6,   4,
            0,   0,   4,  -6,  -8,  12,   4,  -6,
            -70,  36,   0,   0,  -6,   4,  76, -40,
            36, -70,   0,   0,   4,  -6, -40,  76;

    ASSERT_EQ(Mat(H), expected_hess);
}



/** Checks that the MSP's derivative computation satisfies the derivative properties */
TEST(MassSpringProblem, MassSpringProblem2DSparseCheckDerivative){

    int n_unknowns(25);

    SpringElement2DWithLength sewl;
    MassSpringProblem2DSparse msp(sewl, n_unknowns);

    DerivativeChecker checker;
    ASSERT_TRUE(checker.check_all(msp, 0.1, FLT_EPSILON));
}




/** Compares your MSS's energy computation's results with ours */
TEST(MassSpringSystem, EnergyComputation){
    int n_grid_x(20), n_grid_y(20);

    //generate points
    const int n_vertices = (n_grid_x+1)*(n_grid_y+1);

    FunctionBase::Vec points(2*n_vertices);
    for(int i=0; i<n_vertices; ++i) {
        points[2*i] = sin(i * 0.3);
        points[2*i+1] = sin(i * 0.1);
    }

    //initial energy
    AOPT::MassSpringSystemT<AOPT::MassSpringProblem2DSparse> mss(n_grid_x, n_grid_y, 1);
    //set coordinates for graph nodes
    mss.set_spring_graph_points(points);

    auto energy = mss.initial_system_energy();
    ASSERT_FLOAT_EQ(energy, 1033.1414);
}



/** Runs the same problem in both Dense and Sparse form and checks that
 * the Sparse version is (much) faster than the Dense one, as expected */
TEST(MassSpringSystem, DenseVersusSparseSpeedComparison){
    int n_grid_x(100), n_grid_y(100);

    //generate points
    const int n_vertices = (n_grid_x+1)*(n_grid_y+1);
    //uniformly
    FunctionBase::Vec points(2*n_vertices);
    for(int i=0; i<n_vertices; ++i) {
        points[2*i] = i/(n_grid_x+1);
        points[2*i+1] = i%(n_grid_x+1);
    }

    AOPT::StopWatch<> sw;

    //initial energy

    AOPT::MassSpringSystemT<AOPT::MassSpringProblem2DDense> mss_dense(n_grid_x, n_grid_y, 1);
    //set coordinates for graph nodes
    mss_dense.set_spring_graph_points(points);

    std::cout<<"Evaluating the Dense MassSpringSystem..."<<std::endl;
    auto energy = mss_dense.initial_system_energy();
    std::cout<<"MassSpring system energy is "<<energy<<std::endl;

    int n_unknowns = mss_dense.get_problem()->n_unknowns();
    FunctionBase::Vec gradient_dense(n_unknowns);
    mss_dense.get_problem()->eval_gradient(points, gradient_dense);
    std::cout<<"MassSpring system gradient norm is "<<gradient_dense.norm()<<std::endl;

    sw.start();

    AOPT::MassSpringSystemT<AOPT::MassSpringProblem2DDense>::Mat h(n_unknowns, n_unknowns);
    mss_dense.get_problem()->eval_hessian(points, h);
    std::cout<<"MassSpring system hessian norm is "<<h.norm()<<std::endl;

    int duration_dense_ms = sw.stop();
    std::cout<<"Evaluating on DENSE hessian takes: "<<duration_dense_ms/1000.<<"s"<< std::endl;



    AOPT::MassSpringSystemT<AOPT::MassSpringProblem2DSparse> mss_sparse(n_grid_x, n_grid_y, 1);
    mss_sparse.set_spring_graph_points(points);

    std::cout<<"Evaluating the Sparse MassSpringSystem..."<<std::endl;
    energy = mss_sparse.initial_system_energy();
    std::cout<<"MassSpring system energy is "<<energy<<std::endl;

    n_unknowns = mss_sparse.get_problem()->n_unknowns();
    FunctionBase::Vec gradient_sparse(n_unknowns);
    mss_sparse.get_problem()->eval_gradient(points, gradient_sparse);
    std::cout<<"MassSpring system gradient norm is "<<gradient_sparse.norm()<<std::endl;

    sw.start();

    AOPT::MassSpringSystemT<AOPT::MassSpringProblem2DSparse>::SMat sh(n_unknowns, n_unknowns);
    mss_sparse.get_problem()->eval_hessian(points, sh);
    std::cout<<"MassSpring system hessian norm is "<<sh.norm()<<std::endl;

    int duration_sparse_ms = sw.stop();
    std::cout<<"Evaluating on SPARSE hessian takes: "<<duration_sparse_ms/1000.<<"s"<< std::endl;

    std::cout<<" ===> DENSE vs. SPARSE comparison: "<<duration_dense_ms/1000<<"s vs. "<<duration_sparse_ms/1000<<"s"<<std::endl;

    ASSERT_GT(duration_dense_ms, duration_sparse_ms);

}


int main(int _argc, char** _argv){

    testing::InitGoogleTest(&_argc, _argv);
    return RUN_ALL_TESTS();

}


