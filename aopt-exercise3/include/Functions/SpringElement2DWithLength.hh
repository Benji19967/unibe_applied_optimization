#pragma once

#include <FunctionBase/ParametricFunctionBase.hh>

using std::pow;

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
            const double k = _coeffs[0];
            const double l = _coeffs[1];
            return 0.5 * k * pow(pow(_x[0] - _x[2], 2) + pow(_x[1] - _x[3], 2) - pow(l, 2), 2);

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
            double k = _coeffs[0];
            double l = _coeffs[1];

            Vec xa = _x.segment(0, 2);
            Vec xb = _x.segment(2, 2);
            Vec r = xa - xb;
            const double r_norm = r.norm();

            _g.head<2>() = 2 * k * (pow(r_norm, 2) - pow(l, 2)) * r;
            _g.tail<2>() = -2 * k * (pow(r_norm, 2) - pow(l, 2)) * r;
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
            double k = _coeffs[0];
            double l = _coeffs[1];

            Vec xa = _x.segment(0, 2);
            Vec xb = _x.segment(2, 2);
            Vec r = xa - xb;
            const double r_norm = r.norm();

            Eigen::Matrix2d I = Eigen::Matrix2d::Identity();
            Eigen::Matrix2d A = 2 * r * r.transpose() + (r.transpose() * r - pow(l, 2)) * I;
            Eigen::Matrix4d M = Eigen::Matrix4d::Identity();
            M << A, -A, -A, A;
            _H = 2 * k * M;
            //------------------------------------------------------//
        }
    };

//=============================================================================
}


