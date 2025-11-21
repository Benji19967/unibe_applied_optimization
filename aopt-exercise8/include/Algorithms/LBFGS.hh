#pragma once

#include <Eigen/Core>
#include <Algorithms/LineSearch.hh>
#include <FunctionBase/FunctionBaseSparse.hh>

//== NAMESPACES ===============================================================

namespace AOPT {

    //== CLASS DEFINITION =========================================================

    class LBFGS {
    public:
        using Vec = FunctionBaseSparse::Vec;
        using Mat = FunctionBaseSparse::Mat;
        using MapVec = Eigen::Map<Vec>;

        LBFGS(const int _m): m_(_m) {}


        /**
        * @brief solve
        * \param _problem pointer to any function/problem inheriting from FunctionBaseSparse
        *        on which the basic Newton Method will be applied
        * \param _initial_x starting point of the method
        * \param _eps epsilon under which the method stops
        * \param _max_iters maximum iteration of the method*/
        template <class Problem>
        Vec solve(Problem *_problem, const Vec& _initial_x, const double _eps = 1e-4, const int _max_iters = 1000000) {
            std::cout << "******** LBFGS ********" << std::endl;

            int n = _problem->n_unknowns();

            //allocate storage
            init_storage(n);

            //get starting point
            Vec x = _initial_x;

            if(m_ < 1) {
                std::cout<<"\nError: m should be larger than 0!"<<std::endl;
                return x;
            }

            //squared epsilon for stopping criterion
            double e2 = _eps * _eps;

            //allocate gradient storage
            Vec g(n), sk(n), yk(n);

            //initialize k
            int k(0);

            //initialize
            double f = _problem->eval_f(x);
            _problem->eval_gradient(x, g);

            xp_ = x;
            gp_ = g;


            do {
                double g2 = g.squaredNorm();

                //print status
                std::cout << "iter: " << k <<
                          "   obj = " << f <<
                          "   ||g||^2 = " << g2<< std::endl;

                if(g2 < e2) {
                    std::cout<<"Gradient norm converges!"<<std::endl;
                    return x;
                }

                //compute r_
                //------------------------------------------------------//
                //TODO: complete the function
                two_loop_recursion(g, sk, yk, k);
                //------------------------------------------------------//

                //compute the step size
                double t = LineSearch::backtracking_line_search(_problem, x, g, -r_, 1.);
//                double t = LineSearch::wolfe_line_search(_problem, x, g, -r_, 1.);

                if(t < 1e-16) {
                    std::cout<<"The step length is too small!"<<std::endl;
                    return x;
                }

                //update previous function value, x and gradient
                fp_ = f;
                xp_ = x;
                gp_ = g;

                //update x
                x -= t * r_;

                //evaluate current f
                f = _problem->eval_f(x);

                if(k > 0 && fp_ <= f) {
                    std::cout<<"Function value converges!"<<std::endl;
                    return x;
                }

                //current gradient
                _problem->eval_gradient(x, g);

                //update storage
                sk = x - xp_;
                yk = g - gp_;

                //------------------------------------------------------//
                //TODO: complete the function
                update_storage(g, sk, yk, k);
                //------------------------------------------------------//

                k++;

            } while (k < _max_iters);


            return x;
        }



    private:
        void two_loop_recursion(const Vec& _g, const Vec& _sk, const Vec& _yk, const int _k) {
            //------------------------------------------------------//
            //TODO: implement the two-loop recursion as described in the lecture slides
            if(_k == 0) {
                r_ = _g;
                return;
            }

            int n = (int)_g.size();

            //check the curvature condition
            double ys = _sk.dot(_yk);
            if(ys < 0) {
                std::cout<<"Curvature condition violated, search in negative gradient direction!"<<std::endl;
                r_ = _g;
                return;
            }

            double yy = _yk.squaredNorm();

            //initialize r_ with gradient
            r_ = _g;

            //loop range
            int range = std::min(m_, _k);

            //current index
            int j = _k;

            for(int i = 0; i < range; ++i) {
                //in case j is out of scope
                j = (j + m_ - 1)%m_;

                alpha_[j] = rho_[j]*(mat_s_.col(j).dot(r_));
                r_.noalias() -= alpha_[j]*mat_y_.col(j);
            }

            //H_0^k: see formula 6.20 in "Numerical Optimization"
            r_ *= (ys / yy);

            for(int i = 0; i < range; ++i) {
                double beta = rho_[j] * (mat_y_.col(j).dot(r_));
                r_.noalias() += mat_s_.col(j) * (alpha_[j] - beta);

                //in case j is out of scope
                j = (j + 1) % m_;
            }
            //------------------------------------------------------//
        }

        void update_storage(const Vec& _g, const Vec& _sk, const Vec& _yk, const int _k) {
            //------------------------------------------------------//
            //TODO: update the si and yi stored in the mat_s_ and mat_y_ respectively
            //update rho_i stored in rho_
            double ys = _sk.dot(_yk);
            if(ys < 0) {
                std::cout<<"Curvature condition violated, skip updating!"<<std::endl;
                return;
            }

            int cur_idx = _k%m_;

            //update s and y in the storage
            mat_s_.col(cur_idx) = _sk;
            mat_y_.col(cur_idx) = _yk;

            //update \rho
            rho_[cur_idx] = 1./ys;
            //------------------------------------------------------//
        }


        void init_storage(const int _n) {
            mat_y_.resize(_n, m_);
            mat_s_.resize(_n, m_);
            xp_.resize(_n);
            gp_.resize(_n);
            r_.resize(_n);
            rho_.resize(m_);
            alpha_.resize(m_);
        }

    private:
        //variables' name convention follow the lecture slides
        //number of most recent pairs of s and y to store
        int m_;
        //store y
        Mat mat_y_;
        //store s
        Mat mat_s_;
        //previous function value
        double fp_;
        //previous x
        Vec xp_;
        //previous gradient
        Vec gp_;
        //move direction
        Vec r_;

        Vec alpha_;
        Vec rho_;

    };

//=============================================================================
}





