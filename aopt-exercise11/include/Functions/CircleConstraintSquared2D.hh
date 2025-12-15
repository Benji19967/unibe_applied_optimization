#pragma once

#include <FunctionBase/FunctionBaseSparse.hh>

//== NAMESPACES ===================================================================

namespace AOPT {

    //== CLASS DEFINITION =========================================================
    class CircleConstraintSquared2D : public FunctionBaseSparse {
    public:
        using Vec = FunctionBaseSparse::Vec;
        using SMat = FunctionBaseSparse::SMat;
        // f(x,y) = ((x[2*idx]- center_x)^2 + (x[2*idx+1] - center_y)^2 - radius^2)^2
        // constructor
        CircleConstraintSquared2D(const int _n, const int _idx, const double _center_x, const double _center_y, const double _radius)
                : FunctionBaseSparse(), n_(_n), idx_(_idx), center_x_(_center_x), center_y_(_center_y), radius_(_radius) {}

        // number of unknowns
        inline virtual int n_unknowns() override { return n_; }

        // funcion evaluation
        // _x stores the coordinates of all nodes
        inline virtual double eval_f(const Vec &_x) override {
            //------------------------------------------------------//
            //Todo: implement the constraint function value
            double f = std::pow(_x[2*idx_] - center_x_, 2) + std::pow(_x[2*idx_+1] - center_y_, 2) - std::pow(radius_, 2);
            return f*f;
            //------------------------------------------------------//
        }

        // gradient evaluation
        // _g stores the gradient of all nodes
        inline virtual void eval_gradient(const Vec &_x, Vec &_g) override {
            _g.setZero();
            //------------------------------------------------------//
            //Todo: implement the gradient and store in _g
            double f = std::pow(_x[2*idx_] - center_x_, 2) + std::pow(_x[2*idx_+1] - center_y_, 2) - std::pow(radius_, 2);
            _g[2*idx_] = 4.0 * (_x[2*idx_] - center_x_)*f;
            _g[2*idx_+1] = 4.0 * (_x[2*idx_+1] - center_y_)*f;
            //------------------------------------------------------//
        }

        // hessian matrix evaluation
        // _h stores the hessian of all nodes
        inline virtual void eval_hessian(const Vec &_x, SMat &_h) override {
            _h.setZero();
            //------------------------------------------------------//
            //Todo: implement the hessian matrix and store in _h
            double difx = std::pow(_x[2*idx_] - center_x_, 2);
            double dify = std::pow(_x[2*idx_+1] - center_y_, 2);

            _h.insert(2*idx_, 2*idx_) = 12*difx + 4*dify - 4*radius_*radius_;
            _h.insert(2*idx_, 2*idx_+1) = 8*(_x[2*idx_] - center_x_)*(_x[2*idx_+1] - center_y_);
            _h.insert(2*idx_+1, 2*idx_) = _h.coeffRef(2*idx_, 2*idx_+1);
            _h.insert(2*idx_+1, 2*idx_+1) = 4*difx + 12*dify - 4*radius_*radius_;
            //------------------------------------------------------//
        }

    private:
        int n_;
        int idx_;
        double center_x_;
        double center_y_;
        double radius_;
    };

//=============================================================================
}





