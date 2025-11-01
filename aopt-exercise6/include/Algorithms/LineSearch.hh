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
    };
//=============================================================================
}



