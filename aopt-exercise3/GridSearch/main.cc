#include <iostream>
#include <Utils/StopWatch.hh>
#include <Algorithms/GridSearch.hh>

int main(int _argc, const char* _argv[]) {
    if (_argc < 5) {
        std::cout
                << "Usage: input should be 'function index(0: non-convex, 1: 2d quadratic, 2: nd quadratic), lower bound X, upper bound X, number of grid, dimension n (only needed for function 2)', e.g. "
                   "./grid_search 1 -10 10 10" << std::endl;
        return -1;
    }

    //read the input parameters
    int func_index, n_grid, n;
    double x_l, x_u;
    func_index = atoi(_argv[1]);
    x_l = atof(_argv[2]);
    x_u = atof(_argv[3]);
    n_grid = atoi(_argv[4]);

    //do grid search
    AOPT::GridSearch gs(n_grid);
    if (func_index == 0 || func_index == 1) {
        AOPT::GridSearch::Vec vec_xl(2), vec_xu(2);
        //make the search domain square for simplicity
        vec_xl << x_l, x_l;
        vec_xu << x_u, x_u;

        double f_min;
        AOPT::StopWatch<> sw;
        sw.start();
        if (func_index == 0) {
            AOPT::FunctionNonConvex2D fnc2d;
            gs.grid_search_2d(&fnc2d, vec_xl, vec_xu, f_min);
        }
        else if (func_index == 1) {
            AOPT::FunctionQuadratic2D fq2d(-1.);
            gs.grid_search_2d(&fq2d, vec_xl, vec_xu, f_min);
        }

        std::cout << "Grid search for minimum takes: " << sw.stop() / 1000. << "s" << std::endl;
    } else if (func_index == 2) {
        if (_argc < 6) {
            std::cout
                    << "Input should be 'function index, lower bound X, upper bound X, number of grid, dimension n', e.g. "
                       "./grid_search 2 -10 10 10 5" << std::endl;
            return -1;
        }

        //read the input parameters
        n = atoi(_argv[5]);

        if (n < 0) {
            std::cout << "Error: unknown's dimension should be larger than 0!" << std::endl;
            return -1;
        }

        //make the search domain square for simplicity
        Eigen::VectorXd vec_xl(n), vec_xu(n);
        for (int i = 0; i < n; ++i) {
            vec_xl(i) = x_l;
            vec_xu(i) = x_u;
        }

        // generate ND function
        AOPT::FunctionQuadraticND fqnd(n, false);

        double f_min;
        AOPT::StopWatch<> sw;
        sw.start();
        gs.grid_search_nd(&fqnd, vec_xl, vec_xu, f_min);
        std::cout<<"Grid search took: "<<sw.stop()/1000.<<"s"<< std::endl;


    } else {
        std::cout << "Error: Function index is from 0 to 2." << std::endl;
        return -1;
    }

    return 0;
}
