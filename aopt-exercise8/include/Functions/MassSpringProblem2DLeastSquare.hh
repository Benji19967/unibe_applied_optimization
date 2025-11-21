#pragma once

#include <FunctionBase/FunctionBaseSparse.hh>
#include <FunctionBase/ParametricFunctionBase.hh>
#include "ConstrainedSpringElement2DLeastSquare.hh"

//== NAMESPACES ===============================================================

namespace AOPT {


//== CLASS DEFINITION =========================================================

    /** This is very similar to what we did with the non-least square version
     * of the MassSpringProblem except we changed all energy computations into
     * a least square form.
     * What this means is that rather than computing
     *      E = sum Eab, for all springs ab
     * We now compute
     *      E = 1/2 * ||r(x)||^2 = 1/2 * sum rj(x)^2,
     * for all springs (identified by their index j)
     *
     * Now, remember that a key notion is that there is one r function per
     * spring with length, but TWO r functions per elemen without, one for
     * each dimension.
     * See eval_r(...) and eval_jacobian(...) */
    class MassSpringProblem2DLeastSquare : public FunctionBaseSparse {
    public:
        using Vec = FunctionBaseSparse::Vec;
        using Mat = Eigen::MatrixXd;

        using Edge = std::pair<int, int>;
        // sparse matrix type
        using SMat = FunctionBaseSparse::SMat;
        // triplet
        using T = FunctionBaseSparse::T;

        enum SpringElementType {WITHOUT_LENGTH, WITH_LENGTH};


        MassSpringProblem2DLeastSquare(ParametricFunctionBase& _spring, const int _n_unknowns) :
                FunctionBaseSparse(),
                n_(_n_unknowns),
                func_(_spring)
        {
            xe_.resize(func_.n_unknowns());
            ge_.resize(func_.n_unknowns());

            cs_xe_.resize(cse_.n_unknowns());
            cs_ge_.resize(cse_.n_unknowns());
            

            if(func_.n_unknowns() == 2)
                spring_type_ = WITHOUT_LENGTH;
            else if(func_.n_unknowns() == 4)
                spring_type_ = WITH_LENGTH;
        }

        ~MassSpringProblem2DLeastSquare() {}

        virtual int n_unknowns() override {
            return n_;
        }



        /** evaluates the spring elements' energy written in the least square form,
         *      f(x) = 1/2*sum(rj^2(x)).
         *
         * \param _x the problem's springs positions.
         *           It should contain the positions of all nodes of the system.
         *           i.e. (_x[2*i], _x[2*i+1]) is the position of the i-th node
         * \return the sum of the function value */
        virtual double eval_f(const Vec &_x) override {
            double energy(0);

            //------------------------------------------------------//
            /**TODO: assemble function values of all the elements
             * f(x) = 1/2*sum(rj^2(x))
             * Hint: implement eval_r to set r, containing all rj, and then use it to compute the energy */


            Vec r;
            eval_r(_x, r);

            energy = 0.5 * r.squaredNorm();

            //------------------------------------------------------//


            return energy;
        }



        /** The gradient of the least square version.
         *
         * \param _x the problem's springs positions.
         *           It should contain the positions of all nodes of the system.
         *           i.e. (_x[2*i], _x[2*i+1]) is the position of the i-th node
         *
         * \param _g gradint of the objective which is J^T*r, where J is the jacobian matrix */
        virtual void eval_gradient(const Vec &_x, Vec &_g) override {
            _g.resize(n_unknowns());


            //------------------------------------------------------//
            /** TODO: approximate the gradient
             *  Hint: implement eval_r(_x, r) to compute r,
             *        then eval_jacobian(_x, J) to compute J,
             * and then compute the gradient with J^T*r */

            Vec r;
            eval_r(_x, r);

            SMat J;
            eval_jacobian(_x, J);

            _g = J.transpose()*r;
            //------------------------------------------------------//
        }


        /**
         * \param _x the problem's springs positions.
         *           It should contain the positions of all nodes of the system.
         *           i.e. (_x[2*i], _x[2*i+1]) is the position of the i-th node
         * \param _h the hessian matrix of the least square problem, approximated as J^T*J.  **/
        virtual void eval_hessian(const Vec &_x, SMat& _h) override {
            //------------------------------------------------------//
            //approximate the hessian with J^T*J
            _h.resize(n_unknowns(), n_unknowns());

            SMat J;
            eval_jacobian(_x, J);

            _h = J.transpose()*J;
            //------------------------------------------------------//
        }


        void add_spring_element(const int _v_idx0, const int _v_idx1, const double _k = 1., const double _l = 1.) {
            if (2 * _v_idx0 > (int) n_ || _v_idx0 < 0 || 2 * _v_idx1 >= (int) n_ || _v_idx1 < 0)
                std::cout << "Warning: invalid spring element was added... " << _v_idx0 << " " << _v_idx1 << std::endl;
            else {
                springs_.emplace_back(_v_idx0, _v_idx1);
                ks_.push_back(_k);
                ls_.push_back(_l);
            }
        }

        void add_constrained_spring_element(const int _v_idx, const double _w = 1., const double _px = 0., const double _py = 0.) {
            if (2 * _v_idx > (int) n_ || _v_idx < 0)
                std::cout << "Warning: invalid node constraint element was added... " << _v_idx << std::endl;
            else {
                attached_node_indices_.push_back(_v_idx);
                weights_.push_back(_w);
                desired_points_.push_back(_px);
                desired_points_.push_back(_py);
            }
        }

    private:

        /** evaluate the  least square expression rj(x) for all springs
         * and then fills the vector _r */
        void eval_r(const Vec &_x, Vec& _r) {

            //set dimension of vector r, depending on the type of spring
            int num_rj;
            if (spring_type_ == WITHOUT_LENGTH)
                num_rj = 2 * springs_.size();
            else if(spring_type_ == WITH_LENGTH)
                num_rj = springs_.size();

            int dim = num_rj + 2 * attached_node_indices_.size();
            _r.resize(dim);


            //------------------------------------------------------//
            /**TODO: assemble function values of all the elements
             *
             * Hint: if spring_type_ == WITHOUT_LENGTH, it is SpringElement2DLeastSquare
             *       and every spring element has two rj(x), where x is a 2D vector
             *       else if spring_type_ == WITH_LENGTH, it is SpringElement2DWithLengthLeastSquare
             *       and every spring element has one rj(x), where x is a 4D vector */

            if(spring_type_ == WITHOUT_LENGTH) {
                Vec coeff(1);

                for(size_t i=0; i<springs_.size(); ++i) {
                    coeff[0] = ks_[i];

                    xe_[0] = _x[2*springs_[i].first];
                    xe_[1] = _x[2*springs_[i].second];
                    _r[2*i] = func_.eval_f(xe_, coeff);

                    xe_[0] = _x[2*springs_[i].first+1];
                    xe_[1] = _x[2*springs_[i].second+1];
                    _r[2*i+1] = func_.eval_f(xe_, coeff);
                }
            } else if(spring_type_ == WITH_LENGTH){
                Vec coeff(2);

                for(size_t i=0; i<springs_.size(); ++i) {
                    coeff[0] = ks_[i];
                    coeff[1] = ls_[i];

                    xe_[0] = _x[2*springs_[i].first];
                    xe_[1] = _x[2*springs_[i].first+1];
                    xe_[2] = _x[2*springs_[i].second];
                    xe_[3] = _x[2*springs_[i].second+1];

                    _r[i] = func_.eval_f(xe_, coeff);
                }
            }

            //------------------------------------------------------//

            /**help: This assembles
             *       rj(x), for all CONSTRAINED springs
             *       Since those are very similar to the springs without length
             *       in their expression, you can use that to get inspired */
            Vec coeff1(2);
            for(int i=0; i<attached_node_indices_.size(); ++i) {
                coeff1[0] = weights_[i];
                coeff1[1] = desired_points_[2*i];

                cs_xe_[0] = _x[2*attached_node_indices_[i]];
                _r[num_rj+2*i] = cse_.eval_f(cs_xe_, coeff1);

                cs_xe_[0] = _x[2*attached_node_indices_[i]+1];
                coeff1[1] = desired_points_[2*i+1];

                _r[num_rj+2*i+1] = cse_.eval_f(cs_xe_, coeff1);
            }
        }





        void eval_jacobian(const Vec &_x, SMat &_J) {

            //get dimension of vector r
            int num_rj;
            if (spring_type_ == WITHOUT_LENGTH)
                num_rj = 2 * springs_.size();
            else if(spring_type_ == WITH_LENGTH)
                num_rj = springs_.size();

            int dim = num_rj + 2 * attached_node_indices_.size();

            //allocate memory
            _J.resize(dim, n_unknowns());

            std::vector<T> triplets;
            triplets.reserve(4*springs_.size() + 2*attached_node_indices_.size());
            //------------------------------------------------------//
            /**TODO: put local gradient vector to the Jacobian matrix
             * Hint: Here you also have to differentiate between
             *      spring_type_ == WITHOUT_LENGTH --> SpringElement2DLeastSquare
             *      spring_type_ == WITH_LENGTH --> SpringElement2DWithLengthLeastSquare
             *
             * Second hint: use triplets to set up the sparse matrix */


            //SpringElement2DLeastSquare
            if (spring_type_ == WITHOUT_LENGTH) {
                Vec coeff(1);
                for (size_t i = 0; i < springs_.size(); ++i) {
                    coeff[0] = ks_[i];

                    xe_[0] = _x[2 * springs_[i].first];
                    xe_[1] = _x[2 * springs_[i].second];
                    // get local gradient
                    func_.eval_gradient(xe_, coeff, ge_);
                    triplets.emplace_back(2 * i, 2 * springs_[i].first, ge_[0]);
                    triplets.emplace_back(2 * i, 2 * springs_[i].second, ge_[1]);

                    xe_[0] = _x[2 * springs_[i].first + 1];
                    xe_[1] = _x[2 * springs_[i].second + 1];
                    // get local gradient
                    func_.eval_gradient(xe_, coeff, ge_);
                    triplets.emplace_back(2 * i + 1, 2 * springs_[i].first + 1, ge_[0]);
                    triplets.emplace_back(2 * i + 1, 2 * springs_[i].second + 1, ge_[1]);
                }
            } else if(spring_type_ == WITH_LENGTH) { //SpringElement2DWithLengthLeastSquare
                Vec coeff(2);
                for (size_t i = 0; i < springs_.size(); ++i) {
                    coeff[0] = ks_[i];
                    coeff[1] = ls_[i];

                    xe_[0] = _x[2 * springs_[i].first];
                    xe_[1] = _x[2 * springs_[i].first + 1];
                    xe_[2] = _x[2 * springs_[i].second];
                    xe_[3] = _x[2 * springs_[i].second + 1];

                    // get local gradient
                    func_.eval_gradient(xe_, coeff, ge_);

                    triplets.emplace_back(i, 2 * springs_[i].first, ge_[0]);
                    triplets.emplace_back(i, 2 * springs_[i].first + 1, ge_[1]);
                    triplets.emplace_back(i, 2 * springs_[i].second, ge_[2]);
                    triplets.emplace_back(i, 2 * springs_[i].second + 1, ge_[3]);
                }
            }


            Vec coeff1(2);
            for (int i = 0; i < attached_node_indices_.size(); ++i) {
                coeff1[0] = weights_[i];

                cs_xe_[0] = _x[2 * attached_node_indices_[i]];
                coeff1[1] = desired_points_[2 * i];
                cse_.eval_gradient(cs_xe_, coeff1, cs_ge_);
                triplets.emplace_back(num_rj + 2 * i, 2 * attached_node_indices_[i], cs_ge_[0]);


                cs_xe_[0] = _x[2 * attached_node_indices_[i] + 1];
                coeff1[1] = desired_points_[2 * i + 1];
                cse_.eval_gradient(cs_xe_, coeff1, cs_ge_);
                triplets.emplace_back(num_rj + 2 * i + 1, 2 * attached_node_indices_[i] + 1, cs_ge_[0]);
            }

            //------------------------------------------------------//

            //set up sparse matrix
            _J.setFromTriplets(triplets.begin(), triplets.end());
        }


    private:
        int n_;
        std::vector<Edge> springs_;

        ParametricFunctionBase& func_;
        int spring_type_;

        //vector of constants
        std::vector<double> ks_;
        std::vector<double> ls_;

        // coordinates of nodes of spring element
        Vec xe_;
        // gradient of each spring element
        Vec ge_;

        std::vector<int> attached_node_indices_;

        ConstrainedSpringElement2DLeastSquare cse_;

        //vector of constants
        std::vector<double> weights_;
        std::vector<double> desired_points_;

        // coordinates of the node of constrained spring element
        Vec cs_xe_;
        // gradient of each attached node
        Vec cs_ge_;
    };

//=============================================================================
}




