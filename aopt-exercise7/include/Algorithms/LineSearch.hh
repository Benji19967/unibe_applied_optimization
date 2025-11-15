#pragma once

#include <FunctionBase/FunctionBaseSparse.hh>
#include <cstdlib>

using std::abs;

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






        template <class Problem>
        static double wolfe_line_search(Problem *_problem,
                                        const Vec &_x,
                                        const Vec &_g,
                                        const Vec &_dx,
                                        double _t0, double _t_max = 100) {
            
            double t(0);

            //------------------------------------------------------//
            //TODO: implement the line search algorithm that satisfies wolfe condition
            // reference: "Numerical Optimization", "Algorithm 3.5 (Line Search Algorithm)".
    
            // double c1 = 0.25;
            // double c2 = 0.75;
            //
            // double f_prime_0 = _g.transpose() * _dx;
            //
            // // Notation: I use f instead of phi
            // // tm1 means t at i-1
            // // fm1 means f at i-1
            // double t0 = t;
            // double f0 = _problem->eval_f(_x + t * _dx);
            // double tm1 = t;
            // double fm1 = _problem->eval_f(_x + t * _dx);
            //                         
            //
            // t = 20;
            // size_t i(1);
            // while(true) {
            //     double fi = _problem->eval(_x + t * _dx);
            //     if (fi > (f0 + c1 * t * f_prime_0) || (fi >= f0 && i > 1)) {
            //       return zoom(tm1, t);
            //     }
            //     _problem->eval_gradient(_x + t * _dx, _g);
            //     double f_prime_i = _g.transpose() * _dx;
            //     if (std::abs(f_prime_i) <= -c2 * f_prime_0) {
            //       return t;
            //     }
            //     if (f_prime_i >= 0) {
            //       return zoom(t, tm1);
            //     }
            //
            //     // Maybe use a constant and do t = c * t
            //     double LO = t + 0.000001;
            //     double HI = _t_max;
            //     t = LO + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(HI-LO)));
            //
            //     i++;
            // }
            //------------------------------------------------------//

            return t;
        }
        //
        // template <class Problem>
        // static double zoom(Problem *_problem,
        //                    double _t_low, 
        //                    double _t_high,
        //                    const Vec &_x,
        //                    const Vec &_g,
        //                    const Vec &_dx,
        //                    const double f0,
        //                    const double f_prime_0,
        //                    const double c1,
        //                    const double c2
        //                    ) {
        //   while (true) {
        //     double tj = 999; // TODO
        //     double fj = _problem->eval(_x + tj * _dx);
        //     double f_low = _problem->eval(_x + _t_low * _dx);
        //     if (fj > f0 + c1 * tj * f_prime_0 || fj >= f_low) {
        //       _t_high = tj;
        //     } else {
        //       _problem->eval_gradient(_x + tj * _dx, _g);
        //       double f_prime_j = _g.transpose() * _dx;
        //       if (std::abs(f_prime_j) <= -c2 * f_prime_0) {
        //         return tj;
        //       }
        //       if (f_prime_j * (_t_high - _t_low) >= 0) {
        //         _t_high = _t_low;
        //       }
        //       _t_low = tj;
        //     }
        //   }
        // }

    private:
        
        
    };
//=============================================================================
}



