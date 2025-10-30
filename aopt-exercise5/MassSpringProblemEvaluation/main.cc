#include <iostream>
#include <Utils/StopWatch.hh>
#include <MassSpringSystemT.hh>
#include <Utils/RandomNumberGenerator.hh>
#include <Utils/DerivativeChecker.hh>

using namespace AOPT;



int main(int _argc, const char* _argv[]) {
    if(_argc < 5) {
        std::cout << "Usage: input should be 'function index(0: f without length, 1: f with length), "
                     "sparse hessian (0: dense, 1: sparse), number of grids in x, number of grids in y, filename (optional)', e.g. "
                     "./MassSpringProblemEvaluation 1 1 10 10 /usr/spring" << std::endl;
        return -1;
    }

    //read the input parameters
    int func_index, sparse, n_grid_x, n_grid_y;
    func_index = atoi(_argv[1]);
    sparse = atoi(_argv[2]);
    n_grid_x = atoi(_argv[3]);
    n_grid_y = atoi(_argv[4]);

    //generate points
    const int n_vertices = (n_grid_x+1)*(n_grid_y+1);
    //randomly
//    RandomNumberGenerator rng(-1., 1.);
//    FunctionBase::Vec points = rng.get_random_nd_vector(2*n_vertices);
    //uniformly
    FunctionBase::Vec points(2*n_vertices);
    for(int i=0; i<n_vertices; ++i) {
        points[2*i] = i/(n_grid_x+1);
        points[2*i+1] = i%(n_grid_x+1);
    }

    AOPT::StopWatch<> sw;

    //initial energy
    if(!sparse) {
        AOPT::MassSpringSystemT<AOPT::MassSpringProblem2DDense> mss(n_grid_x, n_grid_y, func_index);
        //set coordinates for graph nodes
        mss.set_spring_graph_points(points);

        std::cout<<"Evaluating the Dense MassSpringSystem..."<<std::endl;
        auto energy = mss.initial_system_energy();
        std::cout<<"MassSpring system energy is "<<energy<<std::endl;

        int n_unknowns = mss.get_problem()->n_unknowns();
        FunctionBase::Vec gradient(n_unknowns);
        mss.get_problem()->eval_gradient(points, gradient);
        std::cout<<"MassSpring system gradient norm is "<<gradient.norm()<<std::endl;

        sw.start();

        AOPT::MassSpringSystemT<AOPT::MassSpringProblem2DDense>::Mat h(n_unknowns, n_unknowns);
        mss.get_problem()->eval_hessian(points, h);
        std::cout<<"MassSpring system hessian norm is "<<h.norm()<<std::endl;

        std::cout<<"Evaluating on DENSE hessian takes: "<<sw.stop()/1000.<<"s"<< std::endl;
    } else {
        AOPT::MassSpringSystemT<AOPT::MassSpringProblem2DSparse> mss(n_grid_x, n_grid_y, func_index);
        mss.set_spring_graph_points(points);

        std::cout<<"Evaluating the Sparse MassSpringSystem..."<<std::endl;
        auto energy = mss.initial_system_energy();
        std::cout<<"MassSpring system energy is "<<energy<<std::endl;

        int n_unknowns = mss.get_problem()->n_unknowns();
        FunctionBase::Vec gradient(n_unknowns);
        mss.get_problem()->eval_gradient(points, gradient);
        std::cout<<"MassSpring system gradient norm is "<<gradient.norm()<<std::endl;

        sw.start();

        AOPT::MassSpringSystemT<AOPT::MassSpringProblem2DSparse>::SMat sh(n_unknowns, n_unknowns);
        mss.get_problem()->eval_hessian(points, sh);
        std::cout<<"MassSpring system hessian norm is "<<sh.norm()<<std::endl;

        std::cout<<"Evaluating on SPARSE hessian takes: "<<sw.stop()/1000.<<"s"<< std::endl;
    }


    //optional: output the spring graph to file
    if(_argc == 6) {
        std::string filename(_argv[5]);

        AOPT::MassSpringSystemT<AOPT::MassSpringProblem2DDense> mss(n_grid_x, n_grid_y, func_index);
        mss.save_spring_system(filename.c_str());
    }

    return 0;
}

