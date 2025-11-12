#pragma once

#include <Eigen/Dense>
#include <iostream>

//== NAMESPACES ===============================================================

namespace AOPT {

    //== CLASS DEFINITION =========================================================


    /* This class is very similar to FunctionBase except it allows for additional
     * parameters, which are required for some functions.
     *
     * For instance, one could imagine a set of linear functions a*x + b
     * where the constant a should be allowed to change at each evaluation.
     * In such a case, the constant a would be an argument of either evaluation function.
     * Furthermore, the '_coeffs' argument vectors do not necessarily need to have the
     * same dimension as the input vector _x.
     * One could, for instance, actually give both constants a and b as input,
     * with _coeffs = [a, b], i.e. the dimension of _coeffs being twice the dimension of _x */
    class ParametricFunctionBase {
    public:

        // (dense) vector type
        typedef Eigen::VectorXd Vec;
        // (dense) matrix type
        typedef Eigen::MatrixXd Mat;

        // default constructor
        ParametricFunctionBase() {}

        // defualt destructor
        virtual ~ParametricFunctionBase() {};

        // number of unknowns
        virtual int n_unknowns() = 0;

        // funcion evaluation
        virtual double eval_f(const Vec &_x, const Vec &_coeffs) = 0;

        // gradient evaluation
        virtual void eval_gradient(const Vec &_x, const Vec &_coeffs, Vec &_g) = 0;

        // hessian matrix evaluation
        virtual void eval_hessian(const Vec &_x, const Vec &_coeffs, Mat &_H) = 0;
    };


//=============================================================================
}
