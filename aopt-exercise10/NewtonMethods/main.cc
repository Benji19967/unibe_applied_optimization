#include <iostream>
#include <string>
#include <Utils/StopWatch.hh>

#include <Algorithms/NewtonMethods.hh>
#include <Utils/OptimizationStatistic.hh>
#include <Utils/RandomNumberGenerator.hh>
#include <MassSpringSystemT.hh>



int main(int _argc, const char* _argv[]) {
    if(_argc != 7) {
        std::cout << "Usage: input should be 'newton's method(0: standard newton, 1: projected hessian),"
                     "function index(0: f without length, 1: f with length, 2: f with length with positive local hessian),"
                     " number of grid in x, number of grid in y, max iteration, filename', e.g. "
                     "./NewtonMethods 0 0 2 2 10000 /usr/spring" << std::endl;
        return -1;
    }

    //read the input parameters
    int newton_index, func_index, n_grid_x, n_grid_y, max_iter;
    newton_index = atoi(_argv[1]);
    func_index = atoi(_argv[2]);
    n_grid_x = atoi(_argv[3]);
    n_grid_y = atoi(_argv[4]);
    max_iter = atoi(_argv[5]);

    std::string filename(_argv[6]);


    //construct mass spring system
    AOPT::MassSpringSystemT<AOPT::MassSpringProblem2DSparse> mss(n_grid_x, n_grid_y, func_index);
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

    //save graph before optimization
    std::cout<<"Saving initial spring graph to "<<filename<<"_*.csv"<<std::endl;
    mss.save_spring_system(filename.c_str());

    AOPT::NewtonMethods::Vec x;

    if(newton_index == 0 || newton_index == 2) {
        opt_st->start_recording();
        x = AOPT::NewtonMethods::solve(opt_st.get(), start_pts, 1e-4, max_iter);
        opt_st->print_statistics();
    } else if(newton_index == 1) {
        opt_st->start_recording();
        x = AOPT::NewtonMethods::solve_with_projected_hessian(opt_st.get(), start_pts, 10, 1e-4, max_iter);
        opt_st->print_statistics();
    }

    //set points after optimization
    mss.set_spring_graph_points(x);

    //save graph after optimization
    filename += "_opt";
    std::cout<<"Saving optimized spring graph to "<<filename<<"_*.csv"<<std::endl;
    mss.save_spring_system(filename.c_str());

    return 0;
}