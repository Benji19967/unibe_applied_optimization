#pragma once

#include <FunctionBase/FunctionBaseSparse.hh>
#include "LineSearch.hh"

//== NAMESPACES ===============================================================

namespace AOPT {
    /**
    * @brief NewtonMethods is just a list of functions implementing several variations of the
    * newton's method */
    class NewtonMethods {
    public:
        typedef FunctionBaseSparse::Vec Vec;   // dense vector arbitrary size
        typedef FunctionBaseSparse::Mat Mat;   // dense matrix arbitrary size
        typedef FunctionBaseSparse::T T;        //Triplets
        typedef FunctionBaseSparse::SMat SMat;  // sparse matrix arbitrary size

        /**
         * @brief solve
         * \param _problem pointer to any function/problem inheriting from FunctionBaseSparse
         *        on which the basic Newton Method will be applied
         * \param _initial_x starting point of the method
         * \param _eps epsilon under which the method stops
         * \param _max_iters maximum iteration of the method*/
        static Vec solve(FunctionBaseSparse *_problem, const Vec& _initial_x, const double _eps = 1e-4, const int _max_iters = 1000000) {
            std::cout << "******** Newton Method ********" << std::endl;

            // squared epsilon for stopping criterion
            double e2 = 2* _eps * _eps;

            int n = _problem->n_unknowns();

            // get starting point
            Vec x = _initial_x;

            // allocate gradient storage
            Vec g(n);

            // allocate hessian storage
            SMat H(n, n);

            // allocate search direction vector storage
            Vec delta_x(n);
            int iter(0);


            Eigen::SimplicialLLT<SMat> solver;
  
            //------------------------------------------------------//
            //TODO: implement Newton method
            //fp is the function value of the previous iteration
            double fp = std::numeric_limits<double>::max();

            do {
                ++iter;

                // solve for search direction
                _problem->eval_gradient(x, g);
                _problem->eval_hessian(x, H);

                // H dx = -g
                solver.compute(H);
                if(solver.info() == Eigen::NumericalIssue) {
                    std::cerr << "Warning: LLT factorization has numerical issue!" << std::endl;
                    break;
                }

                delta_x = solver.solve(-g);

                // Newton decrement
                double lambda2 = -g.transpose() * delta_x;

                double f = _problem->eval_f(x);

                // print status
                std::cout << "iter: " << iter <<
                          "   obj = " << f <<
                          "   ||lambda||^2 = " << lambda2 << std::endl;

                if (lambda2 <= e2 || fp <= f) break;



                // step size
                double t = LineSearch::backtracking_line_search(_problem, x, g, delta_x, 1.0);
//            t = LineSearch::wolfe_line_search(_problem, x, g, delta_x, t);

                // update
                x += t * delta_x;
                fp = f;

            } while (iter < _max_iters);

            //------------------------------------------------------//

            return x;
        }

        /**
         * @brief solve with the Projected Hessian method
         * \param _problem pointer to any function/problem inheriting from FunctionBaseSparse.
         *        This problem MUST provide a working eval_hession() function for this method to work.
         *
         * \param _initial_x starting point of the method
         * \param _tau_factor the evolution factor of the tau coefficient
         * \param _eps epsilon under which the method stops
         * \param _max_iters maximum iteration of the method*/
        static Vec solve_with_projected_hessian(FunctionBaseSparse *_problem, const Vec& _initial_x, const double _gamma = 10.0,
                                                const double _eps = 1e-4, const int _max_iters = 1000000) {
            bool converged = false;
            return solve_with_projected_hessian(_problem, converged, _initial_x, _gamma, _eps, _max_iters);
        }

        static Vec solve_with_projected_hessian(FunctionBaseSparse *_problem, bool& _converged, const Vec& _initial_x, const double _gamma = 10.0,
                                                const double _eps = 1e-4, const int _max_iters = 1000000) {
            std::cout << "******** Newton Method with projected hessian ********" << std::endl;

            // squared epsilon for stopping criterion
            double e2 = 2*_eps * _eps;

            int n = _problem->n_unknowns();

            // get starting point
            Vec x = _initial_x;

            // allocate gradient storage
            Vec g(n);

            // allocate hessian storage
            SMat H(n, n);

            // allocate search direction vector storage
            Vec delta_x(n);
            int iter(0);

            // identity and scalar to add positive values to the diagonal
            SMat I(n, n);
            I.setIdentity();

            _converged = false;

            Eigen::SimplicialLLT<SMat> solver;

            //------------------------------------------------------//
            //TODO: implement Newton with projected hessian method
            //Hint: if the factorization fails, then add delta * I to the hessian.
            //      repeat until factorization succeeds (make sure to update delta!)
            double fp = std::numeric_limits<double>::max();

            do {
                ++iter;

                // solve for search direction
                _problem->eval_gradient(x, g);
                _problem->eval_hessian(x, H);

                int cnt = 0;
                double delta = 0.;
                
                std::cout<<" H = "<<H<<std::endl;


                solver.compute(H);
                bool is_not_psd = solver.info() == Eigen::NumericalIssue;
                std::cout<<" psd: "<<!is_not_psd<<std::endl;

                while (is_not_psd && cnt < _max_iters) {
                    if (cnt == 0) {
                        delta = 1e-3 * std::abs(H.diagonal().sum()) / double(n);
                    }
                    H += delta * I;

                    solver.compute(H);
                    is_not_psd = solver.info() == Eigen::NumericalIssue;
                    cnt++;
                    delta *= _gamma;
                }
                delta_x = solver.solve(-g);

                // Newton decrement
                double lambda2 = g.transpose() * (-delta_x);

                double f = _problem->eval_f(x);

                // print status
                std::cout << "iter: " << iter
                          << "   obj = " << f
                          << "   ||lambda||^2 = " << lambda2
                          << "   n_projection_steps = " << cnt << std::endl;

                if (lambda2 <= e2 || fp <= f) {
                    _converged = true;
                    break;
                }

                // step size
                double t = LineSearch::backtracking_line_search(_problem, x, g, delta_x, 1.);

                // update
                x += t * delta_x;
                fp = f;
            } while (iter < _max_iters);

            //------------------------------------------------------//


            return x;
        }



        /**
        * @brief solve problem with the linear equality constraints
        * \param _problem pointer to any function/problem inheriting from FunctionBaseSparse.
        *        This problem MUST provide a working eval_hession() function for this method to work.
        *
        * \param _initial_x starting point of the method
        * \param _A the matrix of constraints Ax=b
        * \param _b the vector of constraints Ax=b
        * \param _eps epsilon under which the method stops
        * \param _max_iters maximum iteration of the method*/
        static Vec solve_equality_constrained(FunctionBaseSparse *_problem, const Vec& _initial_x, const SMat &_A, const Vec &_b,
                                              const double _eps = 1e-4, const int _max_iters = 1000) {
            std::cerr << "******** Equality Constrained Newton ********" << std::endl;

            double eps2 = 2.0 *_eps * _eps;

            // get number of unknowns
            int n = _problem->n_unknowns();

            // get number of constraints
            int p = _A.rows();

            // get starting point
            Vec x = _initial_x;

            // starting point satisfies constraints?
            if ((_A * x - _b).squaredNorm() > eps2)
                project_on_affine(x, _A, _b);

            // allocate gradient storage
            Vec g(n);
            // allocate update storage
            Vec dx(n);
            // allocate hessian storage
            SMat H(n, n);
            // allocate KKT storage
            SMat K(n + p, n + p);
            // allocate rhs storage
            Vec rhs(n + p);
            // allocate solution storage
            Vec dxl(n + p);

            // count number of iterations
            int iter(0);

            Eigen::SparseLU<SMat> solver;
            //------------------------------------------------------//
            //TODO: implement the Newton with equality constraints
            //Hint: the function to set up the KKT matrix is
            //      provided below
            double fp = std::numeric_limits<double>::max();

            while (iter < _max_iters) {
                // get gradient
                _problem->eval_gradient(x, g);

                // get hessian
                _problem->eval_hessian(x, H);

                setup_KKT_matrix(H, _A, K);

                rhs.setZero(n + p);
                rhs.head(n) = -g;

                // solve for constrained Newton step
                solver.compute(K);
                dxl = solver.solve(rhs);

                // extract primal variables
                dx = dxl.head(n);

                // compute Newton decrement
                double lambda2 = -g.transpose() * dx;

                // print status
                double f = _problem->eval_f(x);
                std::cerr << "iter: " << iter <<
                          "   obj = " << f <<
                          "   lambda^2 = " << lambda2 << std::endl;

                // check stopping criterion
                if (lambda2 <= eps2 || f >= fp)
                    break;

                // step size
                double t = LineSearch::backtracking_line_search(_problem, x, g, dx, 1.);

                // update
                x += t * dx;
                fp = f;

                ++iter;
            }
            //------------------------------------------------------//


            return x;
        }

        static Vec solve_equality_constrained_with_infeasible_start(FunctionBaseSparse *_problem, const Vec& _initial_x, const SMat &_A, const Vec &_b,
                                                                    const double _eps = 1e-4, const double _eps_constraints = 1e-4, const int _max_iters = 1000) {
            std::cerr << "******** Equality Constrained Newton with Infeasible Start Point********" << std::endl;
            // get number of unknowns
            int n = _problem->n_unknowns();

            // get number of constraints
            int p = _A.rows();

            // get starting point
            Vec x = _initial_x;

            // allocate gradient storage
            Vec g(n);
            // allocate update storage
            Vec dx(n), dnu(p), w(p), nu(p);

            //initialize nu
            nu.setZero();

            // allocate hessian storage
            SMat H(n, n);
            // allocate KKT storage
            SMat K(n + p, n + p);
            // allocate rhs storage
            Vec rhs(n + p);
            // allocate solution storage
            Vec dxl(n + p);

            //A'*nu
            Vec Anu(n);
            //g+Anu
            Vec rdual(n);
            //Ax-b
            Vec rpri(p);

            //norm of the residual
            double res(0);

            Eigen::SparseLU<SMat> solver;
            //------------------------------------------------------//
            //TODO: implement the Newton with equality constraints
            //count number of iterations
            int iter(0);
            double f, fp = std::numeric_limits<double>::max();

            while (iter < _max_iters) {
                //get gradient
                _problem->eval_gradient(x, g);

                //get function value
                f = _problem->eval_f(x);

                //compute the residual of the primal
                rpri = _A * x - _b;
                double violationSq = rpri.squaredNorm();
                double violation = sqrt(violationSq);

                //dual
                Anu = _A.transpose() * nu;
                rdual = g + Anu;

                //residual
                res = rdual.squaredNorm() + violationSq;
                res = sqrt(res);

                // print status
                std::cerr << "iter: " << iter <<
                          "   obj = " << f <<
                          "   residual = " << res <<" constraint violation = "<< violation<< std::endl;

                if((res < _eps && violation < _eps_constraints)) {
                    std::cout<<" fp = "<<fp<<std::endl;
                    break;
                }


                // get hessian
                _problem->eval_hessian(x, H);

                setup_KKT_matrix(H, _A, K);

                //set right hand side
                rhs.setZero(n + p);
                rhs.head(n) = -rdual;

                //right hand side of infeasible start newton
                rhs.tail(p) = -rpri;

                // solve for constrained Newton step
                solver.compute(K);
                dxl = solver.solve(rhs);

                //get dx
                dx = dxl.head(n);

                //get dnu
                dnu = dxl.tail(p);

                // step size
                double t = LineSearch::backtracking_line_search_newton_with_infeasible_start(_problem, _A, _b, x, nu, dx, dnu, res, 1.);
                // update
                nu += t * dnu;
                x += t * dx;
                fp = f;


                ++iter;
            }
            //------------------------------------------------------//


            return x;
        }


        static Vec solve_equality_constrained_hybrid(FunctionBaseSparse *_problem, const Vec& _initial_x, const SMat &_A, const Vec &_b,
                                                                    const double _eps = 1e-4, const double _eps_constraints = 1e-4, const int _max_iters = 1000) {
            std::cerr << "******** Equality Constrained Newton with hybrid method********" << std::endl;
            // epsilon for newton decrement
            double eps2 = 2.0 *_eps * _eps;


            // get number of unknowns
            int n = _problem->n_unknowns();

            // get number of constraints
            int p = _A.rows();

            // get starting point
            Vec x = _initial_x;

            // allocate gradient storage
            Vec g(n);
            // allocate update storage
            Vec dx(n), dnu(p), w(p), nu(p);

            //initialize nu
            nu.setZero();

            // allocate hessian storage
            SMat H(n, n);
            // allocate KKT storage
            SMat K(n + p, n + p);
            // allocate rhs storage
            Vec rhs(n + p);
            // allocate solution storage
            Vec dxl(n + p);

            //A'*nu
            Vec Anu(n);
            //g+Anu
            Vec rdual(n);
            //Ax-b
            Vec rpri(p);

            //norm of the residual
            double res(0);

            Eigen::SparseLU<SMat> solver;
            //------------------------------------------------------//
            //TODO: implement the Newton with equality constraints
            //count number of iterations
            int iter(0);
            double f, fp = std::numeric_limits<double>::max();

            while (iter < _max_iters) {
                //get function value
                f = _problem->eval_f(x);

                //get gradient
                _problem->eval_gradient(x, g);

                //compute the residual of the primal
                rpri = _A * x - _b;
                double violationSq = rpri.squaredNorm();
                double violation = sqrt(violationSq);

                if(violation >= _eps_constraints) { //infeasible start newton
                    Anu = _A.transpose() * nu;
                    rdual = g + Anu;

                    res = rdual.squaredNorm() + violationSq;
                    res = sqrt(res);

                    // print status
                    std::cerr << "iter: " << iter <<
                              "   obj = " << f <<
                              "   residual = " << res <<" constraint violation = "<< violation << std::endl;

                    if((res < _eps && violation < _eps_constraints))
                        break;
                } else {//feasible start newton
                    //compute Newton decrement
                    double lambda2 = -g.transpose() * dx;

                    //print status
                    std::cerr << "iter: " << iter <<
                              "   obj = " << f <<
                              "   lambda^2 = " << lambda2 << " constraint violation = "<< violation <<std::endl;

                    //check stopping criterion
                    if (lambda2 <= eps2)
                        break;
                }

                // get hessian
                _problem->eval_hessian(x, H);

                setup_KKT_matrix(H, _A, K);

                rhs.setZero(n + p);
                rhs.head(n) = -g;

                //right hand side of infeasible start newton
                if(violation >= _eps_constraints)
                    rhs.tail(p) = -rpri;

                // solve for constrained Newton step
                solver.compute(K);
                dxl = solver.solve(rhs);

                dx = dxl.head(n);

                if(violation >= _eps_constraints) {//use infeasible start newton
                    //get dnu
                    dnu = dxl.tail(p) - nu;

                    // step size
                    double t = LineSearch::backtracking_line_search_newton_with_infeasible_start(_problem, _A, _b, x, nu, dx, dnu, res, 1.);
                    // update
                    nu += t * dnu;
                    x += t * dx;
                } else { //use feasible start newton
                    // step size
                    double t = LineSearch::backtracking_line_search(_problem, x, g, dx, 1.);

                    // update
                    x += t * dx;
                    fp = f;
                }

                ++iter;
            }
            //------------------------------------------------------//


            return x;
        }
        

        
        static void project_on_affine(Vec &_x, const SMat &_A, const Vec &_b) {
            int n = _x.rows();
            int p = _A.rows();

            //------------------------------------------------------//
            //TODO: project x to the hyperplane Ax = b
            Eigen::SparseQR<SMat, Eigen::COLAMDOrdering<int>> sqr;

            SMat A = _A;
            sqr.compute(A);
            if(sqr.info() != Eigen::Success) {
                // decomposition failed
                std::cerr << "SparseQR failed!"<<std::endl;
                return;
            }

            Vec dx(n);
            dx = sqr.solve(_b - _A * _x);

            // get projected
            _x += dx;

            //------------------------------------------------------//

            // check result
            std::cerr << "Constraint violation after projection: " << (_A * _x - _b).squaredNorm() << std::endl;
        }

    private:
        

        static void setup_KKT_matrix(const SMat &_H, const SMat &_A, SMat& _K) {
            const int n  = static_cast<int>(_H.cols());
            const int m  = static_cast<int>(_A.rows());
            const int nf = n+m;

            // set up KKT matrix
            // create sparse matrix
            std::vector<T> trips;
            trips.reserve(_H.nonZeros() + 2*_A.nonZeros());

            // add elements of H
            for (int k=0; k<_H.outerSize(); ++k)
                for (SMat::InnerIterator it(_H,k); it; ++it)
                    trips.emplace_back(static_cast<int>(it.row()),static_cast<int>(it.col()),it.value());

            // add elements of _A
            for (int k=0; k<_A.outerSize(); ++k)
                for (SMat::InnerIterator it(_A,k); it; ++it)
                {
                    // insert _A block below
                    trips.emplace_back(static_cast<int>(it.row())+n,static_cast<int>(it.col()),it.value());

                    // insert _A^T block right
                    trips.emplace_back(static_cast<int>(it.col()),static_cast<int>(it.row())+n,it.value());
                }

            // create KKT matrix
            _K.resize(nf,nf);
            _K.setFromTriplets(trips.begin(), trips.end());
        }
    };

} // namespace AOPT
