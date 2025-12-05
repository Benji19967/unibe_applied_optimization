#pragma once

#include <FunctionBase/FunctionBaseSparse.hh>


//== NAMESPACES ===============================================================

namespace AOPT {


//== CLASS DEFINITION =========================================================

    class AugmentedLagrangianProblem : public FunctionBaseSparse {
    public:
        using Vec = FunctionBaseSparse::Vec;
        using Mat = FunctionBaseSparse::Mat;
        using Edge = std::pair<int, int>;
        // sparse matrix type
        using SMat = FunctionBaseSparse::SMat;
        // triplet
        using T = FunctionBaseSparse::T;

        // default constructor
        AugmentedLagrangianProblem(FunctionBaseSparse* _obj, const std::vector<FunctionBaseSparse*>& _constraints,
                const std::vector<FunctionBaseSparse*>& _squared_constraints, const Vec& _nu, double _mu)
        : FunctionBaseSparse(), obj_(_obj), constraints_(_constraints),
          squared_constraints_(_squared_constraints), nu_(_nu), mu_over_2_(_mu/2.)
        {
            n_ = obj_->n_unknowns();
            g_ = Vec(n_);
            h_ = SMat(n_, n_);
        }

        ~AugmentedLagrangianProblem() {}

        virtual int n_unknowns() override {
            return n_;
        }

        virtual double eval_f(const Vec &_x) override {
            double energy(0);

            //------------------------------------------------------//
            //TODO: accumulate function values (objective function + constraint functions (including squared ones))
            
            
            //------------------------------------------------------//

            return energy;
        }

        virtual void eval_gradient(const Vec &_x, Vec &_g) override {
            _g.resize(n_unknowns());
            _g.setZero();

            //------------------------------------------------------//
            //TODO: accumulate gradients (objective function + constraint functions (including squared ones))
            
            
            //------------------------------------------------------//
        }


        //Hessian in sparse matrix
        virtual void eval_hessian(const Vec &_x, SMat& _h) override {
            _h.resize(n_unknowns(), n_unknowns());
            _h.setZero();

            //------------------------------------------------------//
            //TODO: accumulate hessian matrices (objective function + constraint functions (including squared ones))
            
            
            //------------------------------------------------------//
        }

        //compute constraint function values and store in a vector
        void eval_constraints(const Vec &_x, Vec& _vec_h) {
            for(auto i=0u; i<constraints_.size(); ++i)
                _vec_h[i] = constraints_[i]->eval_f(_x);
        }

        //update nu
        void set_nu(const Vec& _nu) {
            nu_ = _nu;
        };

        //update mu_over_2_
        void set_mu(const double _mu) {
            mu_over_2_ = _mu/2.;
        }


    private:
        int n_;

        FunctionBaseSparse* obj_;
        std::vector<FunctionBaseSparse*> constraints_;
        std::vector<FunctionBaseSparse*> squared_constraints_;

        Vec nu_;
        double mu_over_2_;

        // used as a temporary vector when eval gradient of each constraint
        Vec g_;
        // used as a temporary matrix when eval hessian of each constraint
        SMat h_;
    };

//=============================================================================
}







