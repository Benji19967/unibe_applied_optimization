#pragma once

#include <FunctionBase/ParametricFunctionBase.hh>

//== NAMESPACES ===============================================================

namespace AOPT {


//== CLASS DEFINITION =========================================================

    /* This is for the least square version of the spring element with length, which originally
     * was
     *      E'_ab(x) = 1/2 * k * (((x[0] - x[2])^2 + (x[1] - x[3])^2) - l^2)^2.
     *
     * Written in least square format,
     *      E'_ab(x) = 1/2 * rj^2(x).
     * This class implements the rj(x).
     * Note that because of l^2, we cannot decompose the energy in two evaluations
     * as we did with SpringElement2DLeastSquare and we thus need to pass
     * x = [xa, xb] as we did before with the non-least square versions */
    class SpringElement2DWithLengthLeastSquare : public ParametricFunctionBase {
    public:
        // constructor
        SpringElement2DWithLengthLeastSquare() : ParametricFunctionBase() {}

        // number of unknowns
        inline virtual int n_unknowns() override { return 4; }

        /** evaluates the energy of the function rj(_x)
         * \param _x contains x_a and x_b contiguously,
         *           i.e. _x = [x_a, x_b], i.e. _x is of dimension 4
         * \param _coeffs stores the constants k and l,
         *                i.e. _coeffs[0] = k, _coeffs[1] = l */
        inline virtual double eval_f(const Vec &_x, const Vec &_coeffs) override {
            //------------------------------------------------------//
            //Todo: implement the function rj(x) = sqrt(k) * (||xa-xb||^2 - l^2)
            
            //------------------------------------------------------//
        }

        /** evaluates the gradient of the function rj(_x)
           * \param _x contains x_a and x_b contiguously,
           *           i.e. _x = [x_a, x_b], i.e. _x is of dimension 4
           * \param _coeffs stores the constants k and l,
           *                i.e. _coeffs[0] = k, _coeffs[1] = l
           * \param _g the output gradient, which should also be of dimension 4 */
        inline virtual void eval_gradient(const Vec &_x, const Vec &_coeffs, Vec &_g) override {
            //------------------------------------------------------//
            //Todo: implement the gradient and store in _g
            
            //------------------------------------------------------//
        }

        inline virtual void eval_hessian(const Vec &_x, const Vec &_coeffs, Mat &_H) override {
        }

    };

//=============================================================================
}


