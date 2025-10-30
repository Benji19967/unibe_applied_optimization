#pragma once

#include <FunctionBase/ParametricFunctionBase.hh>

//== NAMESPACES ===============================================================

namespace AOPT {


//== CLASS DEFINITION =========================================================

/* This class evaluates the energy with a spring with length.
 * It is very similar to the SpringElement2D except it has an addition parameter
 * l_ab which represents the length of the spring at rest. */
    class SpringElement2DWithLength : public ParametricFunctionBase {
    public:
        // E'_ab(x) = 1/2 * k * (((x[0] - x[2])^2 + (x[1] - x[3])^2) - l^2)^2
        // constructor
        SpringElement2DWithLength() : ParametricFunctionBase() {}

        // number of unknowns
        inline virtual int n_unknowns() override { return 4; }

        /** evaluates the spring element's energy
         * \param _x contains x_a and x_b contiguously,
         *           i.e. _x = [x_a, x_b], i.e. _x is of dimension 4
         * \param _coeffs stores the constants k and l,
         *                i.e. _coeffs[0] = k, _coeffs[1] = l */
        inline virtual double eval_f(const Vec &_x, const Vec &_coeffs) override {
            //------------------------------------------------------//
            //Todo: implement the function f(x) = 1/2 * k * (((x[0] - x[2])^2 + (x[1] - x[3])^2) - l^2)^2
            double dx = _x[0] - _x[2];
            double dy = _x[1] - _x[3];

            return 0.5 * _coeffs[0] * std::pow((dx*dx + dy*dy) - _coeffs[1]*_coeffs[1], 2);
            //------------------------------------------------------//
        }


        /** evaluates the spring element's energy gradient
         * \param _x contains x_a and x_b contiguously,
         *           i.e. _x = [x_a, x_b], i.e. _x is of dimension 4
         * \param _coeffs stores the constants k and l,
         *                i.e. _coeffs[0] = k, _coeffs[1] = l
         * \param _g the output gradient, which should also be of dimension 4 */
        inline virtual void eval_gradient(const Vec &_x, const Vec &_coeffs, Vec &_g) override {
            //------------------------------------------------------//
            //Todo: implement the gradient and store in _g
            double part = _coeffs[0]*(std::pow((_x[0]-_x[2]),2) + std::pow((_x[1]-_x[3]),2) - _coeffs[1]*_coeffs[1]);
            part *= 2;
            _g[0] = part * (_x[0] - _x[2]);
            _g[1] = part * (_x[1] - _x[3]);
            _g[2] = -_g[0];
            _g[3] = -_g[1];
            
            //------------------------------------------------------//
        }

        /** evaluates the spring element's energy Hessian
         * \param _x contains x_a and x_b contiguously,
         *           i.e. _x = [x_a, x_b], i.e. _x is of dimension 4
         * \param _coeffs stores the constants k and l,
         *                i.e. _coeffs[0] = k, _coeffs[1] = l
         * \param _H the output Hessian, which should be a 4x4 Matrix */
        inline virtual void eval_hessian(const Vec &_x, const Vec &_coeffs, Mat &_H) override {
            //------------------------------------------------------//
            //Todo: implement the hessian matrix and store in _H
            double dx13 = _x[1] - _x[3];
            double dx13sq = std::pow(dx13, 2);
            double dx02 = _x[0] - _x[2];
            double dx02sq = std::pow(dx02, 2);
            double lsq = std::pow(_coeffs[1], 2);
            _H(0,0) = _coeffs[0]*(-2.0*lsq + 6.0*dx02sq + 2.0*dx13sq);
            _H(0,1) = 4.0*_coeffs[0]*dx02*dx13;
            _H(0,2) = -_H(0,0);
            _H(0,3) = -_H(0,1);
            _H(1,0) = _H(0,1);
            _H(1,1) = _coeffs[0]*(-2.0*lsq + 2.0*dx02sq + 6.0*dx13sq);
            _H(1,2) = -_H(0,1);
            _H(1,3) = -_H(1,1);
            _H(2,0) = _H(0,2);
            _H(2,1) = _H(1,2);
            _H(2,2) = _H(0,0);
            _H(2,3) = _H(0,1);
            _H(3,0) = _H(0,3);
            _H(3,1) = _H(1,3);
            _H(3,2) = _H(2,3);
            _H(3,3) = _H(1,1);
            
            //------------------------------------------------------//
        }
    };

//=============================================================================
}


