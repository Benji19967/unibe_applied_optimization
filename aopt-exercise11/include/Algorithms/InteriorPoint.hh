#pragma once

#include <FunctionBase/FunctionBaseSparse.hh>
#include <Functions/InteriorPointProblem.hh>
#include <Algorithms/NewtonMethods.hh>
#include <Utils/OptimizationStatistic.hh>

//== NAMESPACES ===============================================================

namespace AOPT {

    //== CLASS DEFINITION =========================================================
    class InteriorPoint {
    public:
        // LA typedefs
        using Vec = FunctionBaseSparse::Vec;

        static Vec solve(FunctionBaseSparse *_obj, const Vec& _initial_x, const std::vector<FunctionBaseSparse*>& _constraints,
                const double _eps = 1e-4, const double _mu = 10.0, const int _max_iters = 1000) {
            std::cerr << "******** Interior Point ********" << std::endl;

            // construct log-barrier problem
            InteriorPointProblem problem(_obj, _constraints);
            auto opt_st = std::make_unique<AOPT::OptimizationStatistic>(&problem);

            // count number of iterations
            int iter(0);

            // current barrier parameter
            double t = 1.0;

            // number of constraints
            double m = _constraints.size();

            // points
            Vec x = _initial_x;

            //------------------------------------------------------//
            //TODO: implement the interior point method
            //Hint: Use projected newton method to solve for an approximated x.
            
            
            //------------------------------------------------------//

            opt_st->print_statistics();

            return x;
        }
    };
    //=============================================================================

}





