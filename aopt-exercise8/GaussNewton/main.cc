#include <Utils/StopWatch.hh>
#include <iostream>
#include <Algorithms/NewtonMethods.hh>
#include <Functions/MassSpringProblem2DLeastSquare.hh>
#include <Utils/OptimizationStatistic.hh>
#include <MassSpringSystemT.hh>

int main(int _argc, const char* _argv[]) {
    if(_argc != 6) {
        std::cout << "Usage: input should be 'function index(0: f without length, 1: f with length), number of grid in x, number of grid in y, max iteration, filename', e.g. "
                     "./GaussNewton 1 20 20 10000 /usr/spring" << std::endl;
        return -1;
    }

    //read the input parameters
    int func_index, n_grid_x, n_grid_y, m, max_iter;
    func_index = atoi(_argv[1]);
    n_grid_x = atoi(_argv[2]);
    n_grid_y = atoi(_argv[3]);
    max_iter = atoi(_argv[4]);

    std::string filename(_argv[5]);

    //construct mass spring system
    AOPT::MassSpringSystemT<AOPT::MassSpringProblem2DLeastSquare> mss(n_grid_x, n_grid_y, func_index, true);
    //attach spring graph nodes to certain positions
    mss.add_constrained_spring_elements();

    //statistic instance
    auto opt_st = std::make_unique<AOPT::OptimizationStatistic>(mss.get_problem().get());

    //generate the start points
    AOPT::RandomNumberGenerator rng(-10., 10.);
    auto start_pts = rng.get_random_nd_vector(opt_st->n_unknowns());

    //set points
    mss.set_spring_graph_points(start_pts);

    //initial energy
    auto energy = mss.initial_system_energy();
    std::cout<<"\nInitial MassSpring system energy is "<<energy<<std::endl;

    std::cout<<"Saving initial spring graph to "<<filename<<"_*.csv"<<std::endl;
    mss.save_spring_system(filename.c_str());

    filename += "_opt";

    AOPT::NewtonMethods::Vec x = AOPT::NewtonMethods::solve(opt_st.get(), start_pts, 1e-4, max_iter);
    opt_st->print_statistics();

    mss.set_spring_graph_points(x);
    std::cout<<"Saving optimized spring graph to "<<filename<<"_*.csv"<<std::endl;
    mss.save_spring_system(filename.c_str());

    return 0;
}

