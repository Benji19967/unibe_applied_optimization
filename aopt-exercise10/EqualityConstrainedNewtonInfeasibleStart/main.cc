#include <Utils/StopWatch.hh>
#include <iostream>
#include <Algorithms/NewtonMethods.hh>
#include <Utils/OptimizationStatistic.hh>
#include <MassSpringSystemT.hh>

int main(int _argc, const char* _argv[]) {
    if(_argc != 7) {
        std::cout << "Usage: input should be 'newton method(0: feasible start 1: infeasible start, 2: hybrid), function index(0: f without length, 1: f with length, 2: f with length(positive hessian)), number of grid in x, number of grid in y, max iteration, filename', e.g. "
                     "./EqualityConstrainedNewton 0 20 20 10000 /usr/spring" << std::endl;
        return -1;
    }

    //read the input parameters
    int newton_method = atoi(_argv[1]);
    int func_index, n_grid_x, n_grid_y, m, max_iter;
    func_index = atoi(_argv[2]);
    if(func_index == 1)
        std::cout<<"Warning: trying to run Newton's method on non-convex function!"<<std::endl;
    n_grid_x = atoi(_argv[3]);
    n_grid_y = atoi(_argv[4]);
    max_iter = atoi(_argv[5]);

    std::string filename(_argv[6]);

    //initial energy
    AOPT::MassSpringSystemT<AOPT::MassSpringProblem2DSparse> mss(n_grid_x, n_grid_y, func_index);

    AOPT::NewtonMethods::SMat A;
    AOPT::NewtonMethods::Vec b;
    mss.setup_linear_equality_constraints(A, b);

    //statistic instance
    auto opt_st = std::make_unique<AOPT::OptimizationStatistic>(mss.get_problem().get());

    //generate the start points
    AOPT::RandomNumberGenerator rng(-10., 10.);
    auto start_pts = rng.get_random_nd_vector(opt_st->n_unknowns());

    //set points
    mss.set_spring_graph_points(start_pts);
    //initial energy
    auto energy = mss.initial_system_energy();
    std::cout<<"Initial MassSpring system energy is "<<energy<<std::endl;

    //save graph before optimization
    std::cout<<"Saving initial spring graph to "<<filename<<"_*.csv"<<std::endl;
    mss.save_spring_system(filename.c_str());

    filename += "_opt";

    AOPT::NewtonMethods::Vec x(opt_st->n_unknowns());
    if(newton_method == 0)
        x = AOPT::NewtonMethods::solve_equality_constrained(opt_st.get(), start_pts, A, b, 1e-4, max_iter);
    else if(newton_method == 1)
        x = AOPT::NewtonMethods::solve_equality_constrained_with_infeasible_start(opt_st.get(), start_pts, A, b, 1e-6, 1e-8, max_iter);
    else if(newton_method == 2)
        x = AOPT::NewtonMethods::solve_equality_constrained_hybrid(opt_st.get(), start_pts, A, b, 1e-4, 1e-8, max_iter);

    opt_st->print_statistics();

    mss.set_spring_graph_points(x);
    std::cout<<"Saving optimized spring graph to "<<filename<<"_*.csv"<<std::endl;
    mss.save_spring_system(filename.c_str());


    return 0;
}

