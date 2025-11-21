#pragma once

#include <iostream>
#include <FunctionBase/FunctionBase.hh>

//== NAMESPACES ===============================================================

namespace AOPT {

    //== CLASS DEFINITION =========================================================
    /* This implements a generic non-convex function.
     * Its sole purpose is to test the ConvexitTest class */
    class FunctionNonConvex2D final : public FunctionBase {
    public:
        // f(x,y) = (y-x^2)^2+cos^2(4*y)*(1-x)^2+x^2+y^2

        // constructor
        FunctionNonConvex2D() {}

        // number of unknowns
        inline virtual int n_unknowns() { return 2; }

        /** funcion evaluation
         * \param _x the value at which to evaluate the function.
         *           It should be a 2D vector*/
        inline virtual double eval_f(const Vec &_x) {
            //------------------------------------------------------//
            //Todo: implement the function f(x,y) = (y-x^2)^2+cos^2(4*y)*(1-x)^2+x^2+y^2
            return std::pow(_x[1] - _x[0]*_x[0], 2) +
                   std::pow(cos(4*_x[1]), 2) * (1-_x[0])*(1-_x[0]) +
                   std::pow(_x[0], 2) + std::pow(_x[1], 2);
            //------------------------------------------------------//
        }

        // gradient evaluation. Not necessary for this function
        inline virtual void eval_gradient(const Vec &_x, Vec &_g) {}

        // hessian matrix evaluation. Not necessary for this function
        inline virtual void eval_hessian(const Vec &_x, Mat &_H) {}
    };

//=============================================================================
}


