#pragma once

#include "SpringGraph.hh"

#include <Functions/MassSpringProblem2DDense.hh>
#include <Functions/MassSpringProblem2DSparse.hh>

#include <Functions/SpringElement2D.hh>
#include <Functions/SpringElement2DWithLength.hh>


#include <Utils/RandomNumberGenerator.hh>

#include <memory>

//== NAMESPACES ===============================================================

namespace AOPT {

//== CLASS DEFINITION =========================================================

    /* MassSpringSystemT puts together the Mass Spring Problem, the Spring Graph and
     * the Spring Elements and thus represents the full system.
     *
     * Note that it actually is a 'template' and not a class per say.
     * Shortly put, it is a special type of C++ object that is based on other types
     * to become an actual class. In this case, an MSS is defined by the type of
     * its underlying MSP, as shown right above the declaration.
     * The idea is that an MSS would work the same, no matter its internal MSP.
     * So instead of writing a new MSS class for each MSP, we make it a template,
     * which works in a generic fashion.
     *
     * Templates are an advanced feature of the C++ language so feel free to find
     * more information on this subject online.
     * For now, you can simply have a look at the MassSpringSystemT_impl.hh for
     * an example of implementation. */
    template<class MassSpringProblem>
    class MassSpringSystemT {
    public:
        using Point = Eigen::Vector2d;
        using Edge = std::pair<int, int>;
        using Vec = Eigen::VectorXd;
        using Mat = Eigen::MatrixXd;
        using SMat = Eigen::SparseMatrix<double>;

        MassSpringSystemT(const int _n_grid_x, const int _n_grid_y, const int _spring_element_type = 0, const bool _least_square = false) :
                n_grid_x_(_n_grid_x), n_grid_y_(_n_grid_y), n_unknowns_(0), rng_(-1., 1.) {
            setup_spring_graph();
            setup_problem(_spring_element_type, _least_square);
        }

        ~MassSpringSystemT(){}

        enum SpringElementType {WITHOUT_LENGTH, WITH_LENGTH, WITH_LENGTH_PSD_HESS};

        double initial_system_energy() const;

        std::shared_ptr<MassSpringProblem> get_problem() const;

        void set_spring_graph_points(const Vec& _points);

        Vec get_spring_graph_points() const;

        void save_spring_system(const char *_filename) const;

        size_t n_grid_points() const;

        size_t n_edges() const;
        
        //functions of adding constraints
        void add_constrained_spring_elements(const int _scenario = 1);

    private:
        void setup_problem(const int _spring_element_type, const bool _least_square = false);

        void setup_spring_graph();

        int get_grid_index(const int _i, const int _j) const;

    private:
        int n_grid_x_;
        int n_grid_y_;

        int n_unknowns_;

        SpringGraph sg_;
        RandomNumberGenerator rng_;

        SpringElement2D se_;
        SpringElement2DWithLength sewl_;
        

        std::shared_ptr<MassSpringProblem> msp_;
    };

//=============================================================================
}
//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(MASSSPRINGSYSTEM_C)
#define MASSSPRINGSYSTEM_TEMPLATES
#include "MassSpringSystemT_impl.hh"
#endif
//=============================================================================



