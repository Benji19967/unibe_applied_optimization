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
           
        } else { //for normal problem
            if (_spring_element_type == WITH_LENGTH) {
                msp_ = std::make_shared<MassSpringProblem>(sewl_, n_unknowns_);
            } else if (_spring_element_type == WITHOUT_LENGTH) {
                msp_ = std::make_shared<MassSpringProblem>(se_, n_unknowns_);
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
    void MassSpringSystemT<MassSpringProblem>::setup_spring_graph() {
        //------------------------------------------------------//
        //TODO: set up the spring graph of n_grid_x by n_grid_y ()
        for (int i = 0; i < n_grid_x_ + 1; i++) {
            for (int j = 0; j < n_grid_y_ + 1; j++) {
                sg_.add_vertex(Point(i, j));
            }
        }

        // Each vertex has up to 8 neighbours
        // Add neighbours if the indexes are not out of bounds.
        // Avoid adding duplicate edges by ensuring v < u
        int n_rows = n_grid_y_ + 1;
        int n_cols = n_grid_x_ + 1;
        for (int i = 0; i < n_rows; ++i) {
            for (int j = 0; j < n_cols; ++j) {
                int v = i * n_cols + j;

                for (int di = -1; di <= 1; ++di) {
                    for (int dj = -1; dj <= 1; ++dj) {
                        if (di == 0 && dj == 0) continue;

                        int ni = i + di;
                        int nj = j + dj;

                        if (ni >= 0 && ni < n_rows && nj >= 0 && nj < n_cols) {
                            int u = ni * n_cols + nj;

                            if (v < u) {  // add only if v < u to avoid duplicates
                                sg_.add_edge(v, u, 1, sqrt(2));
                            }
                        }
                    }
                }
            }   
        }
    }
    //------------------------------------------------------//



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
