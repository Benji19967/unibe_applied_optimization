#pragma once

#include <FunctionBase/FunctionBaseSparse.hh>
#include "LineSearch.hh"

//== NAMESPACES ===============================================================

namespace AOPT {

    /* Performs a gradient descent on a given problem.
     * This can work with any Problem with a FunctionBase-style interface since
     * the gradient descent method is rather generic mathematically */
    class GradientDescent {
    public:
        typedef FunctionBaseSparse::Vec Vec; ///< Eigen::VectorXd


        /**
         * \param _problem a pointer to a specific Problem, which can be any type that
         *        has the same interface as FunctionBase's (i.e. with eval_f, eval_gradient, etc.)
         * \param _initial_x  the x starting point
         * \param _eps the stopping criterion below which we consider the method
         *             to be done
         * \param _max_iters a capping number of iterations in case you would end-up with a
         *             bad configuration where the successive attemps of finding the
         *             minimum kind of oscillate around the actual minimum without
         *             finding it
         *
         * \return the minimum found by the method. */
        template <class Problem>
        static Vec solve(Problem *_problem, const Vec& _initial_x, const double _eps = 1e-4, const int _max_iters = 1000000) {
            std::cout << "******** Gradient Descent ********" << std::endl;

            // squared epsilon for stopping criterion
            double e2 = _eps * _eps;

            // get starting point
            Vec x = _initial_x;

            // allocate gradient storage
            Vec g(_problem->n_unknowns());
            int iter(0);

            //------------------------------------------------------//
            //TODO: implement the gradient descent
            double f_prev = _problem->eval_f(x);

            for(iter = 0; iter < _max_iters; ++iter) {
                _problem->eval_gradient(x, g);

                if(g.squaredNorm() < e2) {
                    std::cout << "Converged: gradient norm below epsilon after " << iter << " iterations." << std::endl;
                    break;
                }

                double t = LineSearch::backtracking_line_search(_problem, x, g, -g, 1.0);

                if(t == 0) {
                    std::cout << "Line search failed. Stopping after " << iter << " iterations." << std::endl;
                    break;
                }

                x -= t * g;

                double f_new = _problem->eval_f(x);

                if(f_new > f_prev + _eps) {
                    std::cout << "Objective increased. Stopping at iteration " << iter << "." << std::endl;
                    break;
                }

                f_prev = f_new;
            }

            if(iter == _max_iters) {
                std::cout << "Reached maximum iterations (" << _max_iters << ")." << std::endl;
            }

            //------------------------------------------------------//

            return x;
        }
    };
}



