#pragma once

#include <iostream>
#include <FunctionBase/FunctionBase.hh>

//== NAMESPACES ===============================================================

namespace AOPT {

    //== CLASS DEFINITION =========================================================

    /* Implements a n-dimensional quadratic function of the form
     * f(x) = 1/2 x^T A x + b^T x + c  with x in R^n
     *
     * The matrix A and vectors b and c should be given at construction (see below) */
    class FunctionQuadraticND final : public FunctionBase {
    public:

        // random problem constructor.
        // It only takes the dimension of the function and whether
        // the problem should be convex or not
        FunctionQuadraticND(const int _n, bool _convex = true)
                : FunctionBase(), n_(_n), A_(n_, n_), b_(n_), c_(0) {
            initialize_random_problem(10, _convex);
        }

        /* Note: the given matrix A should be square, as suggested by the check below. */
        FunctionQuadraticND(const Mat& _A, const Vec& _b, const double _c)
                : FunctionBase(), A_(0.5*(_A + _A.transpose())), b_(_b), c_(_c) {
            if(_A.rows() != _A.cols())
                std::cerr << "Warning: matrix not square in FunctionQuadraticND" << std::endl;
            n_ = A_.rows();

        }

        // number of unknowns
        inline virtual int n_unknowns() { return n_; }

        /** funcion evaluation
         * \param _x the value at which to evaluate the function.
         *           It should be a ND vector*/
        inline virtual double eval_f(const Vec &_x) {
            //-------------------------------------------------------------------------------//
            //Todo: implement the function 0.5 * (x^T A x) + b^T x + c
            double v = 0.5 * _x.transpose() * A_ * _x;
            v += b_.transpose() * _x + c_;
            return v;
            //-------------------------------------------------------------------------------//
        }

        /** evaluates the quadratic function's gradient
         * \param _x the point on which to evaluate the function
         * \param _g gradient output */
        inline virtual void eval_gradient(const Vec &_x, Vec &_g) {}

        /** evaluates the quadratic function's Hessian
         * \param _x the point on which to evaluate the Hessian.
         *           Actually useless since the Hessian is constant but the method
         *           should still use the same interface as FunctionBase
         * \param _H Hessian output */
        inline virtual void eval_hessian(const Vec &_x, Mat &_H) {}

    private:
        void initialize_random_problem(double _max_val = 10.0, bool _convex = true, const int _random_index = 0)
        {
            std::cerr << "initialize random QP problem with " << n_ << " unknowns... " << std::endl;
            // create random matrix
            std::srand(_random_index);
            for (int i = 0; i < n_; ++i)
                for (int j = 0; j < n_; ++j)
                    A_(i, j) = _max_val * (2.0 * double(std::rand()) / double(RAND_MAX) - 1.0);

            // symmetrize
            if(_convex)
            {
                Mat B = A_.transpose() * A_;
                A_ = B;
            }
            else
            {
                Mat B = 0.5 * (A_.transpose() + A_);
                A_ = B;
            }

            for (int i = 0; i < n_; ++i) {
                // random linear part
                b_[i] = _max_val * (2.0 * double(std::rand()) / double(RAND_MAX) - 1.0);

                // make A_ spd (smallest eigenvalue will be 1.0)
                A_(i, i) += 1.0;
            }

            c_ = 0.0;

            std::cerr << "done!" << std::endl;
        }

    private:

        // quadratic function 1/2 x^T A x + b^T x + c
        int n_;
        Mat A_;
        Vec b_;
        double c_;
    };

//=============================================================================
}




