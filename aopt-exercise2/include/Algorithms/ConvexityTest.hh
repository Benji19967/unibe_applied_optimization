#pragma once

#include <Utils/RandomNumberGenerator.hh>
#include <FunctionBase/FunctionBase.hh>

//== NAMESPACES ===================================================================

namespace AOPT {

    //== CLASS DEFINITION =========================================================

    class ConvexityTest {
    public:
        using Vec = FunctionBase::Vec; ///< Eigen::VectorXd
        using Mat = FunctionBase::Mat; ///< Eigen::MatrixXd

        ConvexityTest() {}

        ~ConvexityTest() {}

    public:

        /** Checks whether the function given as argument is convex or not.
         * If it is not, it should output a point not satisfying the convexity property
         * before returning false.
         * \param _function a function pointer that should be any class inheriting
         * from FunctionBase, e.g. FunctionQuadraticND
         * \param min the minimum value of all tested points' coordinate
         * \param max the maximum value of all tested points' coordinate
         * \param n_evals the number of evaluations/samples tested on the
         *        line between the two points on the function */
        static bool isConvex(FunctionBase* _function, const double min = -1000., const double max = 1000., const int n_evals = 10) {
            const int n = _function->n_unknowns();
            //randomly generate number from min to max
            RandomNumberGenerator rng(min, max);
            
            const int max_sampling_points(1000000);

            const double epsilon = 1e-9; // to avoid false positives due to numerical precision

            if (n_evals < 2) {
                return true; 
            }
            const double delta_t = 1.0 / (n_evals - 1.0);

            for (int i = 0; i < max_sampling_points; ++i) {
                
                Vec x = rng.get_random_nd_vector(n);
                Vec y = rng.get_random_nd_vector(n);
                
                double f_x = _function->eval_f(x);
                double f_y = _function->eval_f(y);
                
                for (int k = 0; k < n_evals; ++k) {
                    
                    double t = k * delta_t;
                    Vec z = (1.0 - t) * y + t * x;
                    // RHS = f((1 - t)y + tx)
                    double f_z = _function->eval_f(z);
                    
                    // RHS = (1 - t)f(y) + t f(x)
                    double rhs = (1.0 - t) * f_y + t * f_x;
                    
                    // Check the convexity condition: f(z) <= (1 - t)f(y) + t f(x)
                    if (f_z > rhs + epsilon) {
                        
                        std::cout << "A counterexample was found" << std::endl;
                        std::cout << "Function is NOT convex." << std::endl;
                        
                        // Output the details using the helper function
                        printPathInfo(y, x, z, t);
                        
                        std::cout << "Values:\n"
                                  << "  LHS (f(z)): " << f_z << "\n"
                                  << "  RHS ((1 - t)f(y) + t f(x)): " << rhs << "\n"
                                  << "Violation: LHS > RHS (" << f_z << " > " << rhs << ")" << std::endl;
                        
                        return false;
                    }
                }
            }
            
            // The loop completes without finding a counterexample, likely convex
            std::cout << "No counterexample found after " << max_sampling_points << " samples." << std::endl;
            std::cout << "Function is likely convex." << std::endl;  
            return true;
        }


    private:
        static void printPathInfo(FunctionBase::Vec p1, FunctionBase::Vec p2, FunctionBase::Vec p, double t) {
            std::cout << "path: p(t) = (1 - t) * p1 + t * p2; \nwith:\n"
                      << "  p1 = (" << p1.transpose() << ")\n"
                      << "  p2 = (" << p2.transpose() << ")\n"
                      << "  p (t = " << t << ") = (" << p.transpose() << ")" << std::endl;
        }

    };




}
