#pragma once

#include <FunctionBase/ParametricFunctionBase.hh>

//== NAMESPACES ===============================================================

namespace AOPT {

//== CLASS DEFINITION =========================================================
    /* This is for the least square version of the spring element, which originally
     * was
     *      E_ab(x) = 1/2 * k * ((x[0] - x[2])^2 + (x[1] - x[3])^2).
     * Written in least square format,
     *      E_ab(x) = 1/2 * (rj^2(xa[0] - xb[0]) + rj^2(xa[1] - xb[1])).
     * This class implements the rj(x) = sqrt(k)*(xa[0,1] - xb[0,1]) and should
     * thus be called twice to get the full energy (once for each dimension)
     *
     * It is a Parametric Function because it requires the elastic constant
     * parameter k_ab for the energy computation. */
    class SpringElement2DLeastSquare : public ParametricFunctionBase {
    public:
        // constructor
        SpringElement2DLeastSquare() : ParametricFunctionBase() {}

        // number of unknowns
        inline virtual int n_unknowns() override { return 2; }


        /** evaluates the energy of the function rj(_x)
        * \param _x contains part of the coordinates of spring node a and b,
        *           e.g. _x = [x_a[0], x_b[0]],  _x is of dimension 2
        * \param _coeffs stores the constant k,
        *                i.e. _coeffs[0] = k */
        inline virtual double eval_f(const Vec &_x, const Vec &_coeffs) override {
            //------------------------------------------------------//
            //Todo: implement the function rj(x) = sqrt(_coeffs[0]) * (x[0] - x[1])
            return sqrt(_coeffs[0]) * (_x[0] - _x[1]);
            //------------------------------------------------------//
        }

        /** evaluates the gradient of the function rj(_x)
         * \param _x contains part of the coordinates of spring node a and b,
         *           e.g. _x = [x_a[0], x_b[0]],  _x is of dimension 2
         * \param _coeffs stores the constant k,
         *                i.e. _coeffs[0] = k
         * \param _g the output gradient, which should also be of dimension 2 */
        inline virtual void eval_gradient(const Vec &_x, const Vec &_coeffs, Vec &_g) override {
            //------------------------------------------------------//
            //Todo: implement the gradient and store in _g
            _g[0] = 1;
            _g[1] = -1;
            _g *= sqrt(_coeffs[0]);
            //------------------------------------------------------//
        }

        inline virtual void eval_hessian(const Vec &_x, const Vec &_coeffs, Mat &_H) override {
        }

    };

//=============================================================================
}


