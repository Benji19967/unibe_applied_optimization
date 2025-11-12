#pragma once

#include <FunctionBase/ParametricFunctionBase.hh>

//== NAMESPACES ===============================================================

namespace AOPT {

    //== CLASS DEFINITION =========================================================

    class ConstrainedSpringElement2D : public ParametricFunctionBase {
    public:
        ConstrainedSpringElement2D() : ParametricFunctionBase() {}

        // number of unknowns
        inline virtual int n_unknowns() final { return 2; }

        /** evaluates the spring element's energy
         * \param _x the spring's current position
         * \param _coeffs _coeffs[0] is the penalty factor,
         *                _coeffs[1] and _coeffs[2] are the desired point coordinates
         * \return the energy of the spring */
        inline virtual double eval_f(const Vec &_x, const Vec &_coeffs) final {
            //------------------------------------------------------//
            //Todo: implement the function f(x) = 1/2 * penalty * ((x[0] - px)^2 + (x[1] - py)^2)
            double dx = _x[0] - _coeffs[1];
            double dy = _x[1] - _coeffs[2];
            return 0.5 * _coeffs[0] * (dx*dx + dy*dy);
            //------------------------------------------------------//
        }

        /** evaluates the spring element's energy gradient
         * \param _x the spring's current position
         * \param _coeffs _coeffs[0] is the penalty factor,
         *                _coeffs[1] and _coeffs[2] are the desired point coordinates
         * \param _g gradient output */
        inline virtual void eval_gradient(const Vec &_x, const Vec &_coeffs, Vec &_g) final {
            //------------------------------------------------------//
            //Todo: implement the gradient and store in _g
            _g[0] = _x[0] - _coeffs[1];
            _g[1] = _x[1] - _coeffs[2];
            _g *= _coeffs[0];
            //------------------------------------------------------//
        }

        /** evaluates the spring element's energy Hessian
         * \param _x the spring's current position
         * \param _coeffs _coeffs[0] is the penalty factor,
         *                _coeffs[1] and _coeffs[2] are the desired point coordinates
         * \param _H Hessian output */
        inline virtual void eval_hessian(const Vec &_x, const Vec &_coeffs, Mat &_H) final {
            //------------------------------------------------------//
            //Todo: implement the hessian matrix and store in _H
            _H << _coeffs[0], 0,
                    0, _coeffs[0];
            //------------------------------------------------------//
        }
    };

    //=============================================================================
}





