#pragma once

#include <FunctionBase/ParametricFunctionBase.hh>

//== NAMESPACES ===============================================================

namespace AOPT {

    //== CLASS DEFINITION =========================================================
    /* This is for the least square version of the constrained spring element. The concept
     * is the same as explained in the comment of the SpringElement2DLeastSquare.
     * */
    class ConstrainedSpringElement2DLeastSquare : public ParametricFunctionBase {
    public:
        ConstrainedSpringElement2DLeastSquare(): ParametricFunctionBase() {}

        // number of unknowns
        inline virtual int n_unknowns() override { return 1; }

        /** evaluates the energy of the function rj(_x)
        * \param _x contains part of the coordinates of spring node a,
        *           e.g. _x = [x_a[0]],  _x is of dimension 1
        * \param _coeffs stores the penalty factor, and part of the desired point coordinates.
        * \return the energy of the function rj(_x) */
        inline virtual double eval_f(const Vec &_x, const Vec &_coeffs) final {
            //------------------------------------------------------//
            //Todo: implement the function rj(x) = sqrt(k) * (x-p)

            //------------------------------------------------------//
        }

        /** evaluates the gradient of the function rj(_x)
         * \param _x contains part of the coordinates of spring node a,
        *           e.g. _x = [x_a[0]],  _x is of dimension 1
         * \param _coeffs _coeffs[0] is the penalty factor,
         *                _coeffs[1] is part of the desired point coordinates
         * \param _g gradient output */
        inline virtual void eval_gradient(const Vec &_x, const Vec &_coeffs, Vec &_g) final {
            //------------------------------------------------------//
            //Todo: implement the gradient and store in _g

            //------------------------------------------------------//
        }

        inline virtual void eval_hessian(const Vec &_x, const Vec &_coeffs, Mat &_H) final {
        }
    };

    //=============================================================================
}








