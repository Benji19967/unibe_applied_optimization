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
            //------------------------------------------------------//
            //Todo:
            //1. check only condition 4 in case there are no constraints
            //2. check inequality constraints (cond. 1.)
            //3. check equality constraints (cond. 1.)
            //4. check lambda (cond. 2.)
            //5. check complementary slackness (cond. 3.)
            //6. check gradient (cond. 4.)
            //------------------------------------------------------//
            if (_equality_constraints.empty() && _inequality_constraints.empty()) {
                // In this case the gradient of the Lagrangian reduces to the gradient of the objective function
                Vec g;
                _objective->eval_gradient(_query_point, g);
                // Norm zero implies zero vector,
                // the best we can do is compare with the tolerance
                if (g.norm() < eps_) {
                    return true;
                } else {
                    return false;
                }
            }
            
            // Size of dual variables should match the number of inequalities and equalities
            if (_inequality_constraints.size() != _lambda.size() || _equality_constraints.size() != _nu.size()) {
                return false;
            }

            // 2. check inequality constraints
            for (auto constraint : _inequality_constraints) {
                double val = constraint->eval_f(_query_point);
                if (val > eps_) {
                    return false;
                }
            }

            // 3. check equality constraints
            for (auto constraint : _equality_constraints) {
                double val = constraint->eval_f(_query_point);
                if (fabs(val) > eps_) {
                    return false;
                }
            }

            // 4. check lambda
            for (size_t i = 0; i < _lambda.size(); ++i) {
                // Comparing with -eps because very small negative numbers should still count as zero
                if (_lambda[i] < -eps_) {
                    return false;
                }
            }

            // 5. check complementary slackness
            for (size_t i = 0; i < _lambda.size(); ++i) {
                double val = _inequality_constraints[i]->eval_f(_query_point);
                if (fabs(val * _lambda[i]) > eps_) {
                    return false;
                }
            }

            // 6. check gradient
            size_t dim = _objective->n_unknowns();
            Vec g(dim); // Total gradient
            Vec gi(dim); // Gradient of the current iteration
            _objective->eval_gradient(_query_point, gi);
            // Add the gradient of the objective function
            g = gi;
            // Add the gradient of the inequality constraints
            for (size_t i = 0; i < _inequality_constraints.size(); i++) {
                _inequality_constraints[i]->eval_gradient(_query_point, gi);
                g += _lambda[i] * gi;
            }
            // Add the gradient of the equality constraints
            for (size_t i = 0; i < _equality_constraints.size(); i++) {
                _equality_constraints[i]->eval_gradient(_query_point, gi);
                g += _nu[i] * gi;
            }
            if (g.norm() > eps_) {
                return false;
            }

            return true;
        }             
            

    private:
        double eps_;
    };
//=============================================================================
}



