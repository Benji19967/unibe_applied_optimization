#pragma once

#include <FunctionBase/FunctionBaseSparse.hh>
#include <Functions/AugmentedLagrangianProblem.hh>
#include <Utils/OptimizationStatistic.hh>
#include <Algorithms/NewtonMethods.hh>
#include "LBFGS.hh"


//== NAMESPACES ===============================================================

namespace AOPT {

    //== CLASS DEFINITION =========================================================
    class AugmentedLagrangian {
    public:
        // LA typedefs
        typedef FunctionBaseSparse::Vec Vec;

        static Vec solve(FunctionBaseSparse *_obj, const Vec& _initial_x, const std::vector<FunctionBaseSparse*>& _constraints, const std::vector<FunctionBaseSparse*>& _squared_constraints,
                const double _eta = 1e-4, const double _tau = 1e-4, const int _max_iters = 20) {
            std::cout << "******** Augmented Lagrangian ********" << std::endl;

            double mu = 10,
            tau = 1./mu,
            eta = std::pow(mu, -0.1),
            hnorm = 0.,
            //previous hnorm
            hnormp = std::numeric_limits<double>::max(),
            tau2 = _tau*_tau;

            //vector of nu and vector of constraint value
            Vec nu(_constraints.size()), h(_constraints.size());
            nu.setZero();
            h.setZero();

            //initialize the augmented lagrangian problem for the unconstrained solver
            AugmentedLagrangianProblem problem(_obj, _constraints, _squared_constraints, nu, mu);
            auto opt_st = std::make_unique<AOPT::OptimizationStatistic>(&problem);

            //get starting point
            Vec x = _initial_x;
            //store previous point
            Vec x_p = x;

            //allocate gradient storage
            Vec g(problem.n_unknowns());

            //------------------------------------------------------//
            //TODO: implement the augmented lagrangian method.
            //Hints: 1. Use projected newton method to solve for an approximated x.
            //          If the maximum iteration is reached and if the norm of the constraints
            //          gets larger, one can say it diverges for simplicity.
            //       2. Use set_mu(...) and set_nu(...) functions in AugmentedLagrangianProblem
            //          class to apply the change of nu and mu
            bool converged = true;
            int iter(0);
            do {
                std::cerr<<"\n---------->Augmented Lagrangian iter: "<<iter<<std::endl;
                //approximate x
                x = NewtonMethods::solve_with_projected_hessian(opt_st.get(), converged, x, 10., std::max(tau*tau, tau2), 1000);

                opt_st->eval_gradient(x, g);

                problem.eval_constraints(x, h);
                hnorm = h.norm();

                if(!converged || hnorm > hnormp) {
                    std::cerr<<"\nDiverged! Restore to previous x.";
                    x = x_p;
                }

                std::cerr<<"\ngradient norm: "<<g.norm()<<"constraint violation: "<<hnorm
                    <<" penalty: "<<mu<<" tau: "<<tau<<" eta: "<<eta<<std::endl;

                if(hnorm <= std::max(eta, _eta)) {
                    if(hnorm <= _eta && g.norm() <= _tau)
                        break;

                    nu += mu*h;
                    eta /= std::pow(mu, 0.9);
                    tau /= mu;
                } else {
                    mu *= 100;
                    eta = std::pow(mu, -0.1);
                    tau = 1./mu;
                }

                //update
                x_p = x;
                problem.set_nu(nu);
                problem.set_mu(mu);
                hnormp = hnorm;

                iter++;
            } while (iter < _max_iters);

            //------------------------------------------------------//

            opt_st->print_statistics();

            return x;
        }
    };
    //=============================================================================

}



