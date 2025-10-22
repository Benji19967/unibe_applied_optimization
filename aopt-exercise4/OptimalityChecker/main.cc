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
    
    //2. inequality constraints

    //3. equality constraints

    //4. set lambdas and vs
   
    //5. set query point
    
    //-------------------------------------------------------------------------------//

    //TODO: uncomment this to test your implementation
    //NOTE: you can change the variables name if you want
    //check
    
    /*AOPT::OptimalityChecker oc;
    if(oc.is_KKT_satisfied(&obj_func, ineq_cons, eq_cons, x, lambda, nu))
        std::cout<<"\nThe query point satisfies the KKT condition."<<std::endl;
    else
        std::cout<<"\nThe query point does NOT satisfy the KKT condition."<<std::endl;
    */

    return 0;
}

