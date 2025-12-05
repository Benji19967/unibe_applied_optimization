#pragma once

#include <FunctionBase/FunctionBaseSparse.hh>


namespace AOPT {

    /**
     * \brief wrapper to use dense parametric functions with non-parametric, sparse methods (e.g. NewtonMethod)*/
    template<class DenseParametricFunction>
    class DenseFunctionWrapper : public FunctionBaseSparse {

    public:

        typedef typename DenseParametricFunction::Vec Vec;
        typedef typename DenseParametricFunction::Mat Mat;
        typedef typename FunctionBaseSparse::SMat SMat;

        DenseFunctionWrapper(const DenseParametricFunction& param_func, const Vec& coeffs, const double non_zero_entry_eps = 1e-7)
        : param_func_(param_func), coeffs_(coeffs), non_zero_entry_eps_(non_zero_entry_eps){}

        virtual int n_unknowns(){
            return param_func_.n_unknowns();
        }

        // funcion evaluation
        virtual double eval_f(const Vec &_x){
            return param_func_.eval_f(_x, coeffs_);
        }

        /** gradient evaluation
         * \param _g output gradient.
         * It should contain the gradient evaluated for each partial derivative.
         * i.e _g[i] = (df/dx_i)(_x)
         * IMPORTANT NOTE: _g should be properly sized at the start of the function */
        virtual void eval_gradient(const Vec &_x, Vec &_g){
            param_func_.eval_gradient(_x, coeffs_, _g);
        }

        /** hessian evaluation
         * \param _H output Hessian matrix.
         * It should contain the hessian evaluated for each partial second-order derivative.
         * i.e _H(i,j) = (d^2f/(dx_i dx_j))(_x)
         * IMPORTANT NOTE: _H should be properly sized at the start of the function */
        virtual void eval_hessian(const Vec &_x, SMat &_H){
            Mat dense_H(_H.rows(), _H.cols());
            param_func_.eval_hessian(_x, coeffs_, dense_H);
            _H = dense_H.sparseView();
        }


    private:
        DenseParametricFunction param_func_;
        Vec coeffs_;
        double non_zero_entry_eps_;

    };
}