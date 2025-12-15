#define MASSSPRINGSYSTEM_C

#include "MassSpringSystemT.hh"

namespace AOPT {


    template<class MassSpringProblem>
    double MassSpringSystemT<MassSpringProblem>::initial_system_energy() const{
        if(msp_ != nullptr) {
            Vec points = get_spring_graph_points();
            return msp_.get()->eval_f(points);
        }

        return -1;
    }

    template<class MassSpringProblem>
    void MassSpringSystemT<MassSpringProblem>::set_spring_graph_points(const Vec& _points) {
        int n_vertices = sg_.n_vertices();

        for(size_t i=0; i<n_vertices; ++i)
            sg_.set_vertex(i, Point(_points[2*i], _points[2*i+1]));
    }

    template<class MassSpringProblem>
    void MassSpringSystemT<MassSpringProblem>::save_spring_system(const char *_filename) const {
        sg_.save_to_files(_filename);
    }

    template<class MassSpringProblem>
    std::shared_ptr<MassSpringProblem> MassSpringSystemT<MassSpringProblem>::get_problem() const {
        return msp_;
    }

    template<class MassSpringProblem>
    void MassSpringSystemT<MassSpringProblem>::setup_problem(const int _spring_element_type, const bool _least_square) {
        //set unknown variable number
        n_unknowns_ = 2 * sg_.n_vertices();

        //initialize the problem pointer
        //for least square problem (Gauss-Newton)
        if(_least_square) {
            if(_spring_element_type == WITH_LENGTH){
                msp_ = std::make_shared<MassSpringProblem>(sewlls_, n_unknowns_);
            } else if(_spring_element_type == WITHOUT_LENGTH){
                msp_ = std::make_shared<MassSpringProblem>(sels_, n_unknowns_);
            }else {
                std::cout << "Error: spring function index should be 0 or 1!" << std::endl;
                return;
            }
        } else { //for normal problem
            if (_spring_element_type == WITH_LENGTH) {
                msp_ = std::make_shared<MassSpringProblem>(sewl_, n_unknowns_);
            } else if (_spring_element_type == WITHOUT_LENGTH) {
                msp_ = std::make_shared<MassSpringProblem>(se_, n_unknowns_);
            } else if (_spring_element_type == WITH_LENGTH_PSD_HESS) {
                msp_ = std::make_shared<MassSpringProblem>(element_wlen_psd_hess_, n_unknowns_);
            } else {
                std::cout << "Error: spring function index should be 0, 1 or 2!" << std::endl;
                return;
            }
        }



        //add spring elements
        for (size_t i = 0; i < sg_.n_edges(); ++i)
            msp_.get()->add_spring_element(sg_.from_vertex(i), sg_.to_vertex(i), sg_.coefficient(i), sg_.length(i));
    }

    template<class MassSpringProblem>
    void MassSpringSystemT<MassSpringProblem>::add_constrained_spring_elements(const int _scenario) {
        //------------------------------------------------------//
        //Todo: add constrained spring elements to the problem
        //implement both scenarios here.
        
        const double w(100000);
        
        if(_scenario == 1) {
            msp_.get()->add_constrained_spring_element(get_grid_index(0,         0),         w, 0,         0);
            msp_.get()->add_constrained_spring_element(get_grid_index(n_grid_x_, n_grid_y_), w, 2 * n_grid_x_, 2 * n_grid_y_);
            msp_.get()->add_constrained_spring_element(get_grid_index(n_grid_x_, 0),         w, 2 * n_grid_x_, 0);
            msp_.get()->add_constrained_spring_element(get_grid_index(0,         n_grid_y_), w, 0,         2 * n_grid_y_);
        } else if (_scenario == 2){
            for (int i = 0; i <= n_grid_x_; ++i) {
                msp_.get()->add_constrained_spring_element(get_grid_index(i, 0),         w, i, 0);
                msp_.get()->add_constrained_spring_element(get_grid_index(i, n_grid_y_), w, i, 2 * n_grid_y_);
            }
        }
        //------------------------------------------------------//
    }

    template<class MassSpringProblem>
    void MassSpringSystemT<MassSpringProblem>::add_constrained_spring_element_for_center_spring_node() {
         //fix center node
         msp_.get()->add_constrained_spring_element(get_grid_index(n_grid_x_/2, n_grid_y_/2), 1e5, 2.5*n_grid_x_, n_grid_y_);
    }



    template<class MassSpringProblem>
    void MassSpringSystemT<MassSpringProblem>::setup_linear_equality_constraints(SMat& _A, Vec& _b) const {
        _A.resize(8, n_unknowns_);
        _A.setZero();

        _b.resize(8);
        _b.setZero();

        _A.insert(0, 2*get_grid_index(0, 0)  ) = 1; _b[0] = 0;
        _A.insert(1, 2*get_grid_index(0, 0)+1) = 1; _b[1] = 0;

        _A.insert(2, 2*get_grid_index(n_grid_x_, n_grid_y_)  ) = 1; _b[2] = 2 * n_grid_x_;
        _A.insert(3, 2*get_grid_index(n_grid_x_, n_grid_y_)+1) = 1; _b[3] = 2 * n_grid_y_;

        _A.insert(4, 2*get_grid_index(n_grid_x_, 0)  ) = 1; _b[4] = 2 * n_grid_x_;
        _A.insert(5, 2*get_grid_index(n_grid_x_, 0)+1) = 1; _b[5] = 0;

        _A.insert(6, 2*get_grid_index(0, n_grid_y_)  ) = 1; _b[6] = 0;
        _A.insert(7, 2*get_grid_index(0, n_grid_y_)+1) = 1; _b[7] = 2 * n_grid_y_;
    }


    template<class MassSpringProblem>
    void MassSpringSystemT<MassSpringProblem>::add_boundary_constraints() {
        double cx = n_grid_x_/2., cy = n_grid_y_/2., radius = (cx+cy);

        //constrain the boundary nodes to the unit circle
        for(int i = 0; i <= n_grid_x_; ++i) {
            circle_constraints_.emplace_back(n_unknowns_, get_grid_index(i, 0), cx, cy, radius);
            circle_constraints_.emplace_back(n_unknowns_, get_grid_index(i, n_grid_y_), cx, cy, radius);

            //quadratic part (for augmented lagrangian method)
            squared_circle_constraints_.emplace_back(n_unknowns_, get_grid_index(i, 0), cx, cy, radius);
            squared_circle_constraints_.emplace_back(n_unknowns_, get_grid_index(i, n_grid_y_), cx, cy, radius);
        }

        for(int j = 1; j < n_grid_y_; ++j) {
            circle_constraints_.emplace_back(n_unknowns_, get_grid_index(0, j), cx, cy, radius);
            circle_constraints_.emplace_back(n_unknowns_, get_grid_index(n_grid_x_, j), cx, cy, radius);

            //quadratic part (for augmented lagrangian method)
            squared_circle_constraints_.emplace_back(n_unknowns_, get_grid_index(0, j), cx, cy, radius);
            squared_circle_constraints_.emplace_back(n_unknowns_, get_grid_index(n_grid_x_, j), cx, cy, radius);
        }


        //store pointers
        for(auto& cst : circle_constraints_)
            constraint_pointers_.push_back(&cst);

        for(auto& sq_cst : squared_circle_constraints_)
            squared_constraint_pointers_.push_back(&sq_cst);
    }


    template<class MassSpringProblem>
    void MassSpringSystemT<MassSpringProblem>::add_area_constraints() {
        //------------------------------------------------------//
        //TODO: implement the constraints to make sure each triangle has positive area
        //convention: node0 - node1 - node2 (ordered counter-clockwise)
        //triangle: 2 o
        //            |  \
        //            |    \
        //          0 o----- o 1
        //
        // warning: make sure your triangles are rotationally consistent.
        //i.e., all nodes indices should be placed as if you rotated the triangle above
        //i.e.
        //                   o 1
        //                /  |
        //              /    |
        //          2 o----- o 0

        
        
        //------------------------------------------------------//


        //store pointers
        for(auto& cst : area_constraints_)
            constraint_pointers_.push_back(&cst);
    }


    template<class MassSpringProblem>
    const std::vector<FunctionBaseSparse*>& MassSpringSystemT<MassSpringProblem>::get_constraints() const {
        return constraint_pointers_;
    }

    template<class MassSpringProblem>
    const std::vector<FunctionBaseSparse*>& MassSpringSystemT<MassSpringProblem>::get_constraints_squared() const {
        return squared_constraint_pointers_;
    }


    template<class MassSpringProblem>
    void MassSpringSystemT<MassSpringProblem>::setup_spring_graph() {

        //------------------------------------------------------//
        //TODO: set up the spring graph of n_grid_x by n_grid_y ()
        //add vertices

        double k(1);
        double l(1);

        for(int j = 0; j <= n_grid_y_; ++j)
            for (int i = 0; i <= n_grid_x_; ++i)
                sg_.add_vertex(Point(i, j));

        //add edges
        for(int j = 0; j < n_grid_y_; ++j) {
            for(int i = 0; i < n_grid_x_; ++i) {
                //horizontal edge
                sg_.add_edge(get_grid_index(i, j), get_grid_index(i+1, j), k, l);
                //vertical edge
                sg_.add_edge(get_grid_index(i, j), get_grid_index(i, j+1), k, l);
                //diagonal edge
                sg_.add_edge(get_grid_index(i, j), get_grid_index(i+1, j+1), k, sqrt(2.) * l);
                //diagonal edge
                sg_.add_edge(get_grid_index(i+1, j), get_grid_index(i, j+1), k, sqrt(2.) * l);
            }
        }

        //add right most
        for(int j = 0; j < n_grid_y_; ++j)
            sg_.add_edge(get_grid_index(n_grid_x_, j), get_grid_index(n_grid_x_, j+1), k, 1.);

        //add top cap
        for(int i = 0; i < n_grid_x_; ++i)
            sg_.add_edge(get_grid_index(i, n_grid_y_), get_grid_index(i+1, n_grid_y_), k, 1.);
        //------------------------------------------------------//
    }

    template<class MassSpringProblem>
    typename MassSpringSystemT<MassSpringProblem>::Vec
    MassSpringSystemT<MassSpringProblem>::get_spring_graph_points() const {
        Vec points(n_unknowns_);
        int n_vertices = sg_.n_vertices();

        for(size_t i=0; i<n_vertices; ++i) {
            points[2*i] = sg_.point(i)[0];
            points[2*i+1] = sg_.point(i)[1];
        }

        return points;
    }

    template<class MassSpringProblem>
    int MassSpringSystemT<MassSpringProblem>::get_grid_index(const int _i, const int _j) const {
        assert(_i<=n_grid_x_ && _j<=n_grid_y_);
        return (n_grid_x_+1)*_j + _i;
    }

    template<class MassSpringProblem>
    size_t MassSpringSystemT<MassSpringProblem>::n_grid_points() const{
        return (n_grid_x_+1) * (n_grid_y_+1);
    }

    template<class MassSpringProblem>
    size_t MassSpringSystemT<MassSpringProblem>::n_edges() const{
        return sg_.n_edges();
    }
}
