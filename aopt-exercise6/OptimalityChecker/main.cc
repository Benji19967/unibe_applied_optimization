#include <Utils/OptimalityChecker.hh>
#include <Functions/FunctionQuadraticND.hh>
#include <vector>
#include <iostream>


int main(int _argc, const char* _argv[]) {
    AOPT::FunctionQuadraticND::Mat A(2, 2);
    A.setZero();
    AOPT::FunctionQuadraticND::Vec b(2);
    b.setZero();

    //-------------------------------------------------------------------------------//
    //Todo: set up the optimization problem by populating the various
    //matrices and vectors coefficients
    //1. set objective function
    A << 2, 0,
        0, -4;
    b << 0, 0;
    AOPT::FunctionQuadraticND obj_func(A, b, 0);

    //2. inequality constraints (reset A and b)
    std::vector<AOPT::FunctionBase*> ineq_cons;
    A << 2, 0,
        0, 0;
    b << 8, -1;
    AOPT::FunctionQuadraticND qnd(A, b, 14);
    ineq_cons.push_back(&qnd);

    A.setZero();
    b << -1, 0;
    AOPT::FunctionQuadraticND qnd2(A, b, -10);
    ineq_cons.push_back(&qnd2);

    //3. equality constraints
    std::vector<AOPT::FunctionBase*> eq_cons;
    A.setZero();
    b << 1, -1;
    AOPT::FunctionQuadraticND qnd3(A, b, 4);
    eq_cons.push_back(&qnd3);

    //Todo: set the lambda and v and query point
    //4. set lambdas and vs
    AOPT::FunctionQuadraticND::Vec lambda(2), nu(1);
    lambda[0] = 4;
    lambda[1] = 0;

    nu[0] = -12;

    //5. set query point
    AOPT::OptimalityChecker::Vec x(2);
    x[0] = -2;
    x[1] = 2;

    //-------------------------------------------------------------------------------//

    //check
    AOPT::OptimalityChecker oc;
    if(oc.is_KKT_satisfied(&obj_func, ineq_cons, eq_cons, x, lambda, nu))
        std::cout<<"\nThe query point satisfies the KKT condition."<<std::endl;
    else
        std::cout<<"\nThe query point does NOT satisfy the KKT condition."<<std::endl;


    return 0;
}

