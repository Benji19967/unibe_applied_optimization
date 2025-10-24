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
    // DON'T FORGET TO MULTIPLY BY 2 TO ACCOUNT FOR THE 1/2 IN THE DEF OF FunctionQuadraticND
    A << 2,  0,
         0, -4;
    b << 0,  0;
    AOPT::FunctionQuadraticND obj_func(A, b, 0);
    
    //2. inequality constraints
    std::vector<AOPT::FunctionBase*> ineq_cons;
    // (x1 + 4)^2 - 2 <= x2 
    // Transforms to: x1^2 + 8*x1 - x2 + 14 <= 0
    // DON'T FORGET TO MULTIPLY BY 2 TO ACCOUNT FOR THE 1/2 IN THE DEF OF FunctionQuadraticND
    A << 2,  0,
         0,  0;
    b << 8, -1;
    AOPT::FunctionQuadraticND qnd1(A, b, 14);
    ineq_cons.push_back(&qnd1);
    // x1 >= -10
    // Transforms to: -x1 - 10 <= 0
    A.setZero();
    b << -1, 0;
    AOPT::FunctionQuadraticND qnd2(A, b, -10);
    ineq_cons.push_back(&qnd2);
    
    //3. equality constraints
    std::vector<AOPT::FunctionBase*> eq_cons;
    // x1 - x2 + 4 = 0
    A.setZero();
    b << 1, -1;
    AOPT::FunctionQuadraticND qnd3(A, b, 4);
    eq_cons.push_back(&qnd3);

    //4. set lambdas and vs
    AOPT::FunctionQuadraticND::Vec lambda(2), nu(1);
    lambda << 4, 0;
    nu << -12;

    //5. set query point
    AOPT::OptimalityChecker::Vec x(2);
    x << -2, 2;
    //-------------------------------------------------------------------------------//

    //TODO: uncomment this to test your implementation
    //NOTE: you can change the variables name if you want
    //check
    AOPT::OptimalityChecker oc;
    if(oc.is_KKT_satisfied(&obj_func, ineq_cons, eq_cons, x, lambda, nu))
        std::cout<<"\nThe query point satisfies the KKT condition."<<std::endl;
    else
        std::cout<<"\nThe query point does NOT satisfy the KKT condition."<<std::endl;

    return 0;
}

