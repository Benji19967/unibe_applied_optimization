#pragma once

#include <FunctionBase/FunctionBase.hh>


namespace AOPT {

    template<class ParametricFunction>
    class ParametricFunctionWrapper : public FunctionBase {

    public:

        typedef typename ParametricFunction::Vec Vec;
        // (dense) matrix type
        typedef typename ParametricFunction::Mat Mat;

        ParametricFunctionWrapper(const ParametricFunction& param_func, const Vec& coeffs) : param_func_(param_func), coeffs_(coeffs){}

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
        virtual void eval_hessian(const Vec &_x, Mat &_H){
            param_func_.eval_hessian(_x, coeffs_, _H);
        }


    private:
        ParametricFunction param_func_;
        Vec coeffs_;

    };
}