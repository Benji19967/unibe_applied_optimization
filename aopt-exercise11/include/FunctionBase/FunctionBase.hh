#pragma once

#include <Eigen/Dense>
#include <iostream>

//== NAMESPACES ===============================================================

namespace AOPT {

    //== CLASS DEFINITION =========================================================

/* Base class for all Functions of this framework.
 * It implements a generic interface consisting of a direct evaluation,
 * an evaluation of its gradient and its hessian.
 *
 * To extend this interface into an actual function, one must declare a child class
 * of FunctionBase and implement all the functions present in the interface below.
 *
 * Using such a generic interface allows to apply generic processes to any function
 * using this interface. The Grid Search method, for instance, consists of evaluating
 * a function at points laid out on a grid. By abstracting the function as it is done
 * here, we allow to use the same GridSearch class for ANY function. */
    class FunctionBase {
    public:

        /* If you are not familiar with Eigen's Matrix types, please refer to
         *  https://eigen.tuxfamily.org/dox/group__TutorialMatrixClass.html
         * or
         *  https://eigen.tuxfamily.org/dox/classEigen_1_1Matrix.html
         * for more details */

        // (dense) vector type
        typedef Eigen::VectorXd Vec;
        // (dense) matrix type
        typedef Eigen::MatrixXd Mat;

        // default constructor
        FunctionBase() {}

        // default destructor
        virtual ~FunctionBase() {};

        // number of unknowns
        virtual int n_unknowns() = 0;

        // funcion evaluation
        virtual double eval_f(const Vec &_x) = 0;

        /** gradient evaluation
         * \param _g output gradient.
         * It should contain the gradient evaluated for each partial derivative.
         * i.e _g[i] = (df/dx_i)(_x)
         * IMPORTANT NOTE: _g should be properly sized at the start of the function */
        virtual void eval_gradient(const Vec &_x, Vec &_g) = 0;

        /** hessian evaluation
         * \param _H output Hessian matrix.
         * It should contain the hessian evaluated for each partial second-order derivative.
         * i.e _H(i,j) = (d^2f/(dx_i dx_j))(_x)
         * IMPORTANT NOTE: _H should be properly sized at the start of the function */
        virtual void eval_hessian(const Vec &_x, Mat &_H) = 0;
    };


//=============================================================================
}

