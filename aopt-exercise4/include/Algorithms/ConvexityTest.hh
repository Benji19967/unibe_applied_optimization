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

            //------------------------------------------------------//
            //Todo: Add your code here
            int count = 0;
            double delta = 1.0 / n_evals;
            while (count < max_sampling_points) {
                // generate a pair of points
                Vec p1 = rng.get_random_nd_vector(n);
                Vec p2 = rng.get_random_nd_vector(n);

                // evaluate the function along the convex combination of p1 and p2
                double f1 = _function->eval_f(p1);
                double f2 = _function->eval_f(p2);

                for (int i = 1; i < n_evals; i++) {

                    double t = delta * i;

                    Vec p = (1.0 - t) * p1 + t * p2;
                    double fp = _function->eval_f(p);
                    double sp = (1.0 - t) * f1 + t * f2;

                    // check that the function f(p) <= (1 - t)*f(p1) + t*f(p2)
                    if (fp > sp) {
                        std::cout << "Function non convexity detected: f(p) = " << fp
                                  << " > (1 - t)*f(p1) + t*f(p2) = " << sp << ". Difference is: "<<fp-sp<<"\n";
                        printPathInfo(p1, p2, p, t);
                        return false;
                    }
                }

                count++;

                if(count % 100000 == 0) {
                    std::cout<<"Processed "<<count<<" pairs of points ..."<<std::endl;
                }
            }
            //------------------------------------------------------//
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
