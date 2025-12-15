#pragma once

#include <FunctionBase/FunctionBaseSparse.hh>

//== NAMESPACES ===================================================================

namespace AOPT {

    //== CLASS DEFINITION =========================================================
    class AreaConstraint2D : public FunctionBaseSparse {
    public:
        // Area constraint: 1/2*det(v_01 | V02) >= eps
        // f = -1/2*(x1 - x0)(y2 - y0) - (x2 - x0)(y1 - y0)) + eps <= 0
        // constructor
        AreaConstraint2D(const int _n, const int _idx0, const double _idx1, const double _idx2, const double _eps = 1e-10)
                : FunctionBaseSparse(), n_(_n), idx0_(_idx0), idx1_(_idx1), idx2_(_idx2), eps_(_eps) {}

        // number of unknowns
        inline virtual int n_unknowns() override { return n_; }

        // function evaluation
        // _x stores the coordinates of all nodes
        inline virtual double eval_f(const Vec &_x) override {
            //------------------------------------------------------//
            //Todo: implement the constraint function value
            // Define the three points (x, y) of the triangle
            double x0(_x[2*idx0_]), y0(_x[2*idx0_+1]);
            double x1(_x[2*idx1_]), y1(_x[2*idx1_+1]);
            double x2(_x[2*idx2_]), y2(_x[2*idx2_+1]);

            double determinant = (x1 - x0)*(y2 - y0) - (y1 - y0)*(x2 - x0);
            return -0.5 * determinant + eps_;
            //------------------------------------------------------//
        }

        // gradient evaluation
        inline virtual void eval_gradient(const Vec &_x, Vec &_g) override {
            _g.setZero();
            //------------------------------------------------------//
            //Todo: implement the gradient and store in _g
            // Define the three points (x, y) of the triangle
            double x0(_x[2*idx0_]), y0(_x[2*idx0_+1]);
            double x1(_x[2*idx1_]), y1(_x[2*idx1_+1]);
            double x2(_x[2*idx2_]), y2(_x[2*idx2_+1]);

            _g[2*idx0_] = y1-y2;
            _g[2*idx0_+1] = x2-x1;
            _g[2*idx1_] = y2-y0;
            _g[2*idx1_+1] = x0-x2;
            _g[2*idx2_] = y0-y1;
            _g[2*idx2_+1] = x1-x0;

            _g *= -0.5;
            //------------------------------------------------------//
        }

        // hessian matrix evaluation
        inline virtual void eval_hessian(const Vec &_x, SMat &_h) override {
            _h.setZero();
            //------------------------------------------------------//
            //Todo: implement the hessian matrix and store in _h

            // Second derivatives of _g[2*idx0_] with respect to y1 and y2
            _h.insert(2*idx0_, 2*idx1_+1) =  1;
            _h.insert(2*idx0_, 2*idx2_+1) = -1;
            // Second derivatives of _g[2*idx0_+1] with respect to x1 and x2
            _h.insert(2*idx0_+1, 2*idx1_) = -1;
            _h.insert(2*idx0_+1, 2*idx2_) =  1;
            // Second derivatives of _g[2*idx1_] with respect to y0 and y1
            _h.insert(2*idx1_, 2*idx0_+1) = -1;
            _h.insert(2*idx1_, 2*idx2_+1) =  1;
            // Second derivatives of _g[2*idx1_+1] with respect to x0 and x2
            _h.insert(2*idx1_+1, 2*idx0_) =  1;
            _h.insert(2*idx1_+1, 2*idx2_) = -1;
            // Second derivatives of _g[2*idx2_] with respect to y0 and y1
            _h.insert(2*idx2_, 2*idx0_+1) =  1;
            _h.insert(2*idx2_, 2*idx1_+1) = -1;
            // Second derivatives of _g[2*idx2_+1] with respect to x0 and x1
            _h.insert(2*idx2_+1, 2*idx0_) = -1;
            _h.insert(2*idx2_+1, 2*idx1_) =  1;

            _h *= -0.5;
            //------------------------------------------------------//
        }

    private:
        int n_;
        // index of the nodes
        int idx0_;
        int idx1_;
        int idx2_;

        double eps_;
    };

//=============================================================================
}






