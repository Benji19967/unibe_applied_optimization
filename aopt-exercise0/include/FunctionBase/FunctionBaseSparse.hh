#pragma once

#include <Eigen/Sparse>
#include <iostream>

//== NAMESPACES ===============================================================

namespace AOPT {

    //== CLASS DEFINITION =========================================================


/* A sparse-matrix-based version of FunctionBase.
 * For more information on sparse matrices, please refer to Eigen's documentation:
 * https://eigen.tuxfamily.org/dox/group__TutorialSparse.html
 *
 * The same arguments made for FunctionBase (include/FunctionBase/FunctionBase.hh)
 * and its functions' arguments apply to FunctionBaseSparse. */
    class FunctionBaseSparse {
    public:



        /* If you are not familiar with Eigen's Matrix types, please refer to
         *  https://eigen.tuxfamily.org/dox/group__TutorialMatrixClass.html
         * or
         *  https://eigen.tuxfamily.org/dox/classEigen_1_1Matrix.html
         * for more details */
        typedef Eigen::VectorXd Vec;///< (dense) vector type
        typedef Eigen::MatrixXd Mat;///< (dense) matrix type

        using SMat = Eigen::SparseMatrix<double>;///< sparse matrix type
        using T = Eigen::Triplet<double>;///< triplet, used to fill the sparse matrices

        // default constructor
        FunctionBaseSparse() {}

        // default destructor
        virtual ~FunctionBaseSparse() {};

        // number of unknowns
        virtual int n_unknowns() = 0;

        // funcion evaluation
        virtual double eval_f(const Vec &_x) = 0;

        // gradient evaluation
        virtual void eval_gradient(const Vec &_x, Vec &_g) = 0;

        // hessian matrix evaluation
        virtual void eval_hessian(const Vec &_x, SMat& _h) = 0;
    };


//=============================================================================
}

