#pragma once

#include <FunctionBase/FunctionBase.hh>
#include <vector>
//== NAMESPACES ===============================================================

namespace AOPT {

//== CLASS DEFINITION =========================================================
    class OptimalityChecker {
    public:
        using Vec = Eigen::VectorXd;

        OptimalityChecker(const double _epsilon = 1e-13) : eps_(_epsilon) {}


        /** Checks whether a particular optimization problem satisfies the KKT conditions
         *
         * \param _objective pointer to the objective function, which should be any function
         *        inheriting from FunctionBase
         *
         * \param _inequality_constraints an array containing the inequalities,
         *        each representend as a function inheriting from FunctionBase,
         *        such that _inequality_constraints[i].eval_f(x) <= 0 for all i
         *        if x is in the feasible set
         *
         * \param _equality_constraints an array containing the equalities,
         *        each representend as a function inheriting from FunctionBase,
         *        such that _equality_constraints[i].eval_f(x) == 0 for all i
         *        if x is in the feasible set
         *
         * \param _query_point the point at which the KKT conditions should be checked
         *
         * \param _lambda the lambda coefficients of the dual problem.
         *        It should be of the same dimension as _inequality_constraints since
         *        there is one lambda factor per inequality in the dual problem
         * \param _nu the nu coefficients of the dual problem.
         *        It should be of the same dimension as _equality_constraints since
         *        there is one nu factor per equality in the dual problem
         * */
        bool is_KKT_satisfied(FunctionBase *_objective, const std::vector<FunctionBase *>& _inequality_constraints,
                              const std::vector<FunctionBase *>& _equality_constraints,
                              const Vec& _query_point, const Vec& _lambda, const Vec& _nu) {
            if(_equality_constraints.empty() && _inequality_constraints.empty()) {
                Vec g;
                _objective->eval_gradient(_query_point, g);
                if(g.norm() < eps_)
                    return true;
                else
                    return false;
            } else {
                if(_inequality_constraints.size() != _lambda.size() || _equality_constraints.size() != _nu.size()) {
                    std::cout<<"\nDual variable size should be the same as the constraint size!"<<std::endl;
                    return false;
                }

                //1. check inequality constraints
                for(auto ineq_cons : _inequality_constraints){
                    double val = ineq_cons->eval_f(_query_point);
                    if(val > eps_){
                        std::cout<<"\nInequality Constraint violated! "<<std::endl;
                        return false;
                    }
                }

                //2. check equality constraints
                for(auto eq_cons : _equality_constraints) {
                    double val = eq_cons->eval_f(_query_point);
                    if(fabs(val) > eps_){
                        std::cout<<"\nEq Constraint violated! "<<std::endl;
                        return false;
                    }
                }

                //3. check lambda
                for(size_t i=0; i<_lambda.size(); ++i) {
                    if (_lambda[i] < -eps_) {
                        std::cout<<"\nLambda is smaller than 0! "<<std::endl;
                        return false;
                    }
                }

                //4. check complementary slackness
                for(size_t i=0; i<_lambda.size(); ++i) {
                    double val = _inequality_constraints[i]->eval_f(_query_point);
                    if(fabs(val * _lambda[i]) > eps_) {
                        std::cout<<"\nComplementary slackness violated! "<<std::endl;
                        return false;
                    }
                }

                //5. check gradient
                int dim = _objective->n_unknowns();
                Vec g(dim);

                Vec gi(dim);
                _objective->eval_gradient(_query_point, gi);
                g = gi;

                for(int i=0; i<_inequality_constraints.size(); i++) {
                    _inequality_constraints[i]->eval_gradient(_query_point, gi);
                    g += _lambda[i] * gi;
                }

                for(int i=0; i<_equality_constraints.size(); i++) {
                    _equality_constraints[i]->eval_gradient(_query_point, gi);
                    g += _nu[i] * gi;
                }

                if(g.norm() > eps_) {
                    std::cout<<"\nGradient condition violated! "<<g.norm()<<std::endl;
                    return false;
                }


                return true;
            }
        }

    private:
        double eps_;
    };
//=============================================================================
}



