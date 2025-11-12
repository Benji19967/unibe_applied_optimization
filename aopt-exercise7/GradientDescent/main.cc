#include <iostream>
#include <string>
#include <Utils/StopWatch.hh>

#include <Algorithms/GradientDescent.hh>
#include <Utils/OptimizationStatistic.hh>
#include <Utils/RandomNumberGenerator.hh>
#include <MassSpringSystemT.hh>
#include <Utils/DerivativeChecker.hh>


//initialize different start points
std::vector<AOPT::GradientDescent::Vec> get_start_points(int n_grid_x, int n_grid_y) {
    std::vector<AOPT::GradientDescent::Vec> start_pts;

    //vertex number
    const int n_vertices = (n_grid_x+1)*(n_grid_y+1);

    //generate the first start points
//    AOPT::RandomNumberGenerator rng1(-0.2, 0.2);
//    //
//    AOPT::GradientDescent::Vec points(2*n_vertices);
//    for(int i=0; i<n_vertices; ++i) {
//        points[2*i] = i/(n_grid_x+1);
//        points[2*i+1] = i%(n_grid_x+1);
//    }
//
//    points += rng1.get_random_nd_vector(2*n_vertices);
//
//    start_pts.push_back(points);

    //generate the second start points
    AOPT::RandomNumberGenerator rng2(-10., 10.);
    start_pts.push_back(rng2.get_random_nd_vector(2*n_vertices));


    return start_pts;
}


int main(int _argc, const char* _argv[]) {
    if(_argc != 7) {
        std::cout << "Usage: input should be 'function index(0: f without length, 1: f with length),"
                     "constrained spring scenario  (1 or 2 )"
                     "number of grid in x, number of grid in y, max iteration, filename', e.g. "
                     "./GradientDescent 0 2 2 10000 /usr/spring" << std::endl;
        return -1;
    }

    //read the input parameters
    int func_index, scenario, n_grid_x, n_grid_y, max_iter;
    func_index = atoi(_argv[1]);
    scenario = atoi(_argv[2]);
    n_grid_x = atoi(_argv[3]);
    n_grid_y = atoi(_argv[4]);
    max_iter = atoi(_argv[5]);

    std::string filename(_argv[6]);


    //construct mass spring system
    AOPT::MassSpringSystemT<AOPT::MassSpringProblem2DSparse> mss(n_grid_x, n_grid_y, func_index);
    //attach spring graph nodes to certain positions
    mss.add_constrained_spring_elements(scenario);

    //statistic instance
    auto opt_st = std::make_unique<AOPT::OptimizationStatistic>(mss.get_problem().get());

    //initialize start points
    auto start_points = get_start_points(n_grid_x, n_grid_y);

    //test on two different start points
    for(auto i=0u; i<start_points.size(); ++i) {
        //set points
        mss.set_spring_graph_points(start_points[i]);
        //initial energy
        auto energy = mss.initial_system_energy();
        std::cout<<"\nInitial MassSpring system energy is "<<energy<<std::endl;

        //save graph before optimization
        std::string fn = filename + std::string(std::to_string(i+1));
        std::cout<<"Saving initial spring graph to "<<fn<<"_*.csv"<<std::endl;
        mss.save_spring_system(fn.c_str());

        opt_st->start_recording();
        AOPT::GradientDescent::Vec x = AOPT::GradientDescent::solve(opt_st.get(), start_points[i], 1e-4, max_iter);
        opt_st->print_statistics();

        //set points after optimization
        mss.set_spring_graph_points(x);

        //save graph after optimization
        fn += "_opt";
        std::cout<<"Saving optimized spring graph to "<<fn<<"_*.csv"<<std::endl;
        mss.save_spring_system(fn.c_str());
    }


    return 0;
}
