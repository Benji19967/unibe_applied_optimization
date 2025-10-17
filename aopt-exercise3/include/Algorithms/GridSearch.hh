#pragma once

#include <Functions/FunctionQuadratic2D.hh>
#include <Functions/FunctionQuadraticND.hh>
#include <Functions/FunctionNonConvex2D.hh>
#include <vector>

//== NAMESPACES ===================================================================

namespace AOPT {

    //== CLASS DEFINITION =========================================================
    class GridSearch {
    public:
        using Vec = FunctionBase::Vec;
        using Mat = FunctionBase::Mat;

        GridSearch(const int _grid_dim = 10) : n_grid_(_grid_dim){}
        ~GridSearch() {}

    public:

        /** Evaluation of a 2D function over the whole grid to find its minimum
         *
         * \param _func a pointer to any 2D function inheriting from FunctionBase
         *              (see include/FunctionBase/FunctionBase.hh)
         * \param _x_l the coordinates of the lower corner of the grid
         * \param _x_u the coordinates of the upper corner of the grid.
         *             _x_l and _x_u together define a square in which the grid lies
         * \return 0 if all went well, -1 if not.*/
        int grid_search_2d(FunctionBase* _func, const Vec& _x_l, const Vec& _x_u, double& _f_min) const {
            std::cout<<"Grid searching the minimum of a 2-D function..."<<std::endl;
            double f = 0., fmin = std::numeric_limits<double>::max();
            
            Vec x_min(2);
            
            //------------------------------------------------------//
            //Todo: implement the 2d version of the grid search
            // algorithm to find minimum value of _func between _x_l and _x_u
            //------------------------------------------------------//

            //set-up the evaluation point and delta vector
            Vec x(2), dx(2);
            dx = (_x_u - _x_l) / n_grid_;

            //going through first dimension
            for(int i=0; i <= n_grid_; ++i) {
                x[0] = _x_l[0] + dx[0]*i;

                //and 2nd dimension
                for(int j=0; j <= n_grid_; ++j) {
                    x[1] = _x_l[1] + dx[1]*j;

                    //evaluate function at x
                    f = _func->eval_f(x);

                    //and update the minimum if needed
                    if(f < fmin) {
                        fmin = f;
                        x_min = x;
                    }
                }
            }
            //------------------------------------------------------//
            _f_min = fmin;
            std::cout<<"Minimum value of the function is: "<<fmin<<" at x:\n"<<x_min<<std::endl;
            return 0;
        }



        /** Evaluation of an ND function over the whole grid to find its minimum
         *  using an iterative approach
         *
         * \param _func a pointer to any ND function inheriting from FunctionBase
         *              (see include/FunctionBase/FunctionBase.hh)
         * \param _x_l the coordinates of the lower corner of the grid
         * \param _x_u the coordinates of the upper corner of the grid.
         *             _x_l and _x_u together define an ND cuboid in which the grid lies
         * \return 0 if all went well, -1 if not.*/
        int grid_search_nd(FunctionBase* _func, const Vec& _x_l, const Vec& _x_u, double& _f_min) const {
            int n = _func->n_unknowns();
            if (_x_l.size() != n || _x_u.size() != n) {
                std::cout << "Error: input limits are not of correct dimension!" << std::endl;
                return -1;
            }
            std::cout << "Grid searching the minimum of a " << n << "-D function..." << std::endl;

            double f_min = std::numeric_limits<double>::max();
            Vec x_min(n);
            //------------------------------------------------------//
            //Todo: implement the nd version of the grid search
            // algorithm to find minimum value of a nd quadratic function
            // set f_min with the minimum, which is then stored in the referenced argument _f_min

            // iteratively walk through all grid cells with n-dimensional indices idx
            std::vector<int> idx(n, 0);
            // current grid coordiantes x and grid size dx
            Vec x(n), dx(n);
            dx = (_x_u - _x_l) / double(n_grid_);
            // start from lower limits
            x = _x_l;
            
            unsigned int nt = std::pow(n_grid_ + 1, n);

            for (unsigned int i = 0; i < nt; ++i) {
                // evaluate function
                double f = _func->eval_f(x);
                // if better than best found -> store new solution
                if (f < f_min) {
                    f_min = f;
                    x_min = x;
                }

                // move to next grid cell
                idx[0]++;
                x[0] += dx[0];

                // update grid indices and grid coordinates when necessary
                for (int j = 0; j < n-1; ++j) {
                    if (idx[j] > n_grid_) {
                        // reset current dimension
                        idx[j] = 0;
                        x[j] = _x_l[j];

                        // increase index and coordinate of next dimension
                        idx[j + 1]++;
                        x[j + 1] += dx[j + 1];
                    } else break; // next dimension can only increase if current one increased
                }
            }
            //------------------------------------------------------//
            _f_min = f_min;
            std::cout << "Minimum value of the function is: " << f_min << " at x:\n" << x_min << std::endl;

            return 0;
        }



    private:
        int n_grid_; ///< the number of cells on one side of every dimension
                    ///< e.g., with a ND grid, you would have (n_grid_)^n cells and thus (n_grid_ + 1)^n evaluation points

    };

    //=============================================================================
}





