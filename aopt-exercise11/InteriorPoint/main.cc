#include <iostream>
#include <Utils/OptimizationStatistic.hh>
#include <MassSpringSystemT.hh>
#include <Algorithms/InteriorPoint.hh>

int main(int _argc, const char* _argv[]) {
    if(_argc != 7) {
        std::cout << "Usage: input should be 'function index (0: f without length, 1: f with length, 2: f with length(positive hessian)), "
                     "test case (0: without area constraints, 1: with area constraints)"
                     "number of grid in x, number of grid in y, max iteration, filename', e.g. "
                     "./InteriorPoint 0 20 20 10000 /usr/spring" << std::endl;
        return -1;
    }

    //read the input parameters
    int func_index, test_index, n_grid_x, n_grid_y, max_iter;
    func_index = atoi(_argv[1]);
    test_index = atoi(_argv[2]);
    n_grid_x = atoi(_argv[3]);
    n_grid_y = atoi(_argv[4]);
    max_iter = atoi(_argv[5]);

    std::string filename(_argv[6]);

    AOPT::MassSpringSystemT<AOPT::MassSpringProblem2DSparse> mss(n_grid_x, n_grid_y, func_index);

    //fix four corner nodes and one interior node to certain position
    mss.add_constrained_spring_elements();
    mss.add_constrained_spring_element_for_center_spring_node();

    if(test_index == 1) {
        //make sure each triangle has positive area
        mss.add_area_constraints();
    }


    //set start points
    auto start_pts = mss.get_spring_graph_points();
    mss.set_spring_graph_points(start_pts);

    std::cout<<"Saving initial spring graph to "<<filename<<"_*.csv"<<std::endl;
    mss.save_spring_system(filename.c_str());

    //solve
    AOPT::InteriorPoint::Vec x = AOPT::InteriorPoint::solve(mss.get_problem().get(), start_pts, mss.get_constraints(), 1e-4, 10, max_iter);

    mss.set_spring_graph_points(x);
    filename += "_opt";
    std::cout<<"Saving optimized spring graph to "<<filename<<"_*.csv"<<std::endl;
    mss.save_spring_system(filename.c_str());

    return 0;
}

