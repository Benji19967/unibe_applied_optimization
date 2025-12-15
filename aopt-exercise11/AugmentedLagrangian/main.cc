#include <Utils/StopWatch.hh>
#include <iostream>
#include <Utils/OptimizationStatistic.hh>
#include <MassSpringSystemT.hh>
#include <Algorithms/AugmentedLagrangian.hh>

int main(int _argc, const char* _argv[]) {
    if(_argc != 6) {
        std::cout << "Usage: input should be 'function index (0: f without length, 1: f with length, 2: f with length(positive hessian)), "
                     "number of grid in x, number of grid in y, max iteration, filename, e.g. "
                     "./AugmentedLagrangian 0 20 20 10000 /usr/spring" << std::endl;
        return -1;
    }

    //read the input parameters
    int func_index, cnst_index, n_grid_x, n_grid_y, max_iter;
    func_index = atoi(_argv[1]);
    n_grid_x = atoi(_argv[2]);
    n_grid_y = atoi(_argv[3]);
    max_iter = atoi(_argv[4]);

    std::string filename(_argv[5]);

    //initial energy
    AOPT::MassSpringSystemT<AOPT::MassSpringProblem2DSparse> mss(n_grid_x, n_grid_y, func_index);
    mss.add_boundary_constraints();

    auto opt_st = std::make_unique<AOPT::OptimizationStatistic>(mss.get_problem().get());

    //generate the start points
    //for spring element without length, it can be trapped at local
    // minimum where all points lie on same point on the circle
//    AOPT::RandomNumberGenerator rng(-100., 100.);
//    auto start_pts = rng.get_random_nd_vector(opt_st->n_unknowns());

    auto start_pts = mss.get_spring_graph_points();

    //set points
    mss.set_spring_graph_points(start_pts);

    std::cout<<"Saving initial spring graph to "<<filename<<"_*.csv"<<std::endl;
    mss.save_spring_system(filename.c_str());

    filename += "_opt";


    AOPT::AugmentedLagrangian::Vec x = AOPT::AugmentedLagrangian::solve(opt_st.get(),
                                                                        start_pts,
                                                                        mss.get_constraints(),
                                                                        mss.get_constraints_squared(),
                                                                        1e-4, 1e-4, max_iter);

    mss.set_spring_graph_points(x);
    std::cout<<"Saving optimized spring graph to "<<filename<<"_*.csv"<<std::endl;
    mss.save_spring_system(filename.c_str());

    return 0;
}

