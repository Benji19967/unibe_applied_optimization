#pragma once

#include <FunctionBase/FunctionBaseSparse.hh>

//== NAMESPACES ===============================================================

namespace AOPT {

    //== CLASS DEFINITION =========================================================

    class LineSearch {
    public:
        typedef FunctionBaseSparse::Vec Vec;
        typedef FunctionBaseSparse::SMat SMat;

        /** Back-tracking line search method
         *
         * \param _problem a pointer to a specific Problem, which can be any type that
         *        has the same interface as FunctionBase's (i.e. with eval_f, eval_gradient, etc.)
         * \param _x starting point of the method. Should be of the same dimension as the Problem's
         * \param _g gradient at the starting point.
         * \param _dx delta x
         * \param _t0 inital step of the method
         * \param _alpha and _tau variation constant, as stated by the method's definition
         * \return the final step t computed by the back-tracking line search */
        template <class Problem>
        static double backtracking_line_search(Problem *_problem,
                                               const Vec &_x,
                                               const Vec &_g,
                                               const Vec &_dx,
                                               const double _t0,
                                               const double _alpha = 0.5,
                                               const double _tau = 0.75) {

            double t(0);

            //------------------------------------------------------//
            //TODO: implement the backtracking line search algorithm
            t = _t0;

            // pre-compute objective
            double fx = _problem->eval_f(_x);

            // pre-compute dot product
            double gtdx = _g.transpose() * _dx;

            // make sure dx points to a descent direction
            if (gtdx > 0) {
                std::cerr << "dx is in the direction that increases the function value. gTdx = "<<gtdx << std::endl;
                return t;
            }

            // backtracking (stable in case of NAN)
            int i = 0;
            while (!(_problem->eval_f(_x + t * _dx) <= fx + _alpha * t * gtdx) && i<1000) {
                t *= _tau;
                i++;
            }

            //------------------------------------------------------//

            return t;
        }



        /** Back-tracking line search for infeasible start Newton's method
        *
        * \param _problem a pointer to a specific Problem, which can be any type that
        *        has the same interface as FunctionBase's (i.e. with eval_f, eval_gradient, etc.)
        * \param _A matrix of the lhs of linear equality constraints
        * \param _b vector of the rhs of linear equality constraints
        * \param _x primal starting point of the method
        * \param _nu dual starting point of the method
        * \param _dx delta x
        * \param _dnu delta nu
        * \param _initial_res initial residual
        * \param _t0 inital step of the method
        * \param _alpha and _beta variation constant, as stated by the method's definition
        * \return the final step t computed by the back-tracking line search */
        template <class Problem>
        static double backtracking_line_search_newton_with_infeasible_start(Problem *_problem,
                                                                            const SMat& _A,
                                                                            const Vec& _b,
                                                                            const Vec &_x,
                                                                            const Vec &_nu,
                                                                            const Vec &_dx,
                                                                            const Vec &_dnu,
                                                                            const double _initial_res,
                                                                            const double _t0,
                                                                            const double _alpha = 0.2,
                                                                            const double _beta = 0.9) {
            //------------------------------------------------------//

            double t = _t0;

                        //TODO: implement the algorithm
                        int n = _problem->n_unknowns();
                        Vec ga(n), Anua(n), rpr_a(n), rdual_a(n);
                        double resa(0);

                        // backtracking
                        while (1) {
                            _problem->eval_gradient(_x + t * _dx, ga);
                            Anua = _A.transpose() * (_nu + t * _dnu);

                            rpr_a = _A * (_x + t * _dx) - _b;
                            rdual_a = ga + Anua;

                            resa = rpr_a.squaredNorm() + rdual_a.squaredNorm();
                            resa = sqrt(resa);

                            if (resa <= (1. - _alpha * t) * _initial_res)
                                break;

                            if(t < 1e-7) {
                                return t;
                            }
                            t *= _beta;
                        }

            //------------------------------------------------------//

            return t;
        }



        template <class Problem>
        static double wolfe_line_search(Problem *_problem,
                                        const Vec &_x,
                                        const Vec &_g,
                                        const Vec &_dx,
                                        double _t0, double _t_max = 100) {
            //------------------------------------------------------//
            //TODO: implement the line search algorithm that satisfies wolfe condition
            // reference: "Numerical Optimization", "Algorithm 3.5 (Line Search Algorithm)".

            double t = _t0;
            // reference: "Numerical Optimization", "Algorithm 3.5 (Line Search Algorithm)".
            // \alpha  ->  t
            // \phi    ->   f
            // \phi'   ->   dg

            // increase rate
            const double inc = 2.;

            // save the function value at the t = 0
            const double fx_init = _problem->eval_f(_x);
            // projection of gradient on the search direction
            const double dg_init = _g.dot(_dx);
            // make sure dx points to a descent direction
            if (dg_init > 0) {
                std::cerr << "dx is in the direction that increases the function value." << std::endl;
                return t;
            }

            const double dg_test = 1e-4 * dg_init,
                    dg_wolfe = -0.9 * dg_init;

            // first stage:
            // begins with a trial estimate t, and keeps increasing it until it finds either
            // an acceptable step length or an interval that brackets the desired step lengths
            Vec g(_problem->n_unknowns());
            int iter = 1;
            double tp = 0, fxp = fx_init, dgp = dg_init;
            do {
                double fx = _problem->eval_f(_x + t * _dx);
                _problem->eval_gradient(_x + t * _dx, g);
                const double dg = g.dot(_dx);

                if (fx - fx_init > t * dg_test || (1 < iter && fx >= fxp)) {
                    t = zoom(_problem, _x, _dx, fx_init, dg_init, fx, fxp, dg, dgp, t, tp);
                    return t;
                }

                if (std::abs(dg) <= dg_wolfe)
                    return t;

                if (dg >= 0) {
                    t = zoom(_problem, _x, _dx, fx_init, dg_init, fxp, fx, dgp, dg, tp, t);
                    return t;
                }

                tp = t;
                fxp = fx;
                dgp = dg;

                // increase t by 2 in (t, t_max)
                if(t*inc > _t_max) {
                    t += _t_max;
                    t /= 2.;
                } else
                    t *= inc;

                iter++;
            } while (iter <= 20);

            //------------------------------------------------------//

            return t;
        }


    private:
        template <class Problem>
        static double zoom(Problem *_problem,
                           const Vec &_x,
                           const Vec &_dx,
                           double _fx_init,
                           double _dg_init,
                           double _fx_hi,
                           double _fx_lo,
                           double _dg_hi,
                           double _dg_lo,
                           double _thi, double _tlo) {
            // second stage:
            // successively decreases the size of the interval until
            // an acceptable step length is identified.
            const double dg_test = 1e-4 * _dg_init,
                    dg_wolfe = -0.9 * _dg_init;

            int iter(0);
            Vec g(_problem->n_unknowns());
            double t(1.);

            double fx_hi, fx_lo, dg_hi, dg_lo;

            do {
                if (iter == 0) {
                    fx_hi = _fx_hi;
                    dg_hi = _dg_hi;
                    fx_lo = _fx_lo;
                    dg_lo = _dg_lo;
                } else {
                    fx_hi = _problem->eval_f(_x + _thi * _dx);
                    fx_lo = _problem->eval_f(_x + _tlo * _dx);
                    _problem->eval_gradient(_x + _tlo * _dx, g);
                    dg_lo = g.dot(_dx);
                }


                // use {fx_lo, fx_hi, dg_lo} to make a quadric interpolation of
                // the function said interpolation is used to estimate the minimum
                //
                // polynomial: p (x) = a*(x - _t)² + b
                // conditions: p (thi) = fx_hi
                //             p (tlo) = fx_lo
                //             p'(tlo) = dg_lo
                t = (fx_hi - fx_lo) * _tlo - (_thi * _thi - _tlo * _tlo) * dg_lo / 2;
                double div = (fx_hi - fx_lo) - (_thi - _tlo) * dg_lo;
                if(div == 0)
                    t = _tlo;
                else
                    t /= div;

                // if interpolation fails, bisection is used
                if (t <= std::min(_tlo, _thi) || t >= std::max(_tlo, _thi))
                    t = (_tlo + _thi) / 2;

                double fx = _problem->eval_f(_x + t * _dx);
                _problem->eval_gradient(_x + t * _dx, g);

                const double dg = g.dot(_dx);

                if (fx - _fx_init > t * dg_test || fx >= fx_lo) {
                    if (t == _thi) {
                        std::cerr << "t equals to thi, possibly due to insufficient numeric precision." << std::endl;
                        return t;
                    }

                    _thi = t;
                    fx_hi = fx;
                    dg_hi = dg;
                } else {
                    if (std::abs(dg) <= dg_wolfe)
                        return t;

                    if (dg * (_thi - _tlo) >= 0) {
                        _thi = _tlo;
                        fx_hi = fx_lo;
                        dg_hi = dg_lo;
                    }

                    if (t == _tlo) {
                        std::cerr << "t equals to tlo, possibly due to insufficient numeric precision." << std::endl;
                        return t;
                    }

                    _tlo = t;
                    fx_lo = fx;
                    dg_lo = dg;
                }

                iter++;
            } while (iter <= 20);

            return t;
        }
    };
//=============================================================================
}



