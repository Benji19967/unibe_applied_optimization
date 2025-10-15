#pragma once

#include <FunctionBase/ParametricFunctionBase.hh>
#include <cmath>

using std::pow

    //== NAMESPACES
    //===============================================================

    namespace AOPT {

  //== CLASS DEFINITION
  //=========================================================

  /* This function evaluates the energy of an ideal spring with no length
   * going from x_a to x_b.
   * It is a Parametric Function because it requires the elastic constant
   * parameter k_ab for the energy computation. */
  class SpringElement2D : public ParametricFunctionBase {
  public:
    // E_ab(x) = 1/2 * k * ((x[0] - x[2])^2 + (x[1] - x[3])^2)
    // constructor
    SpringElement2D() : ParametricFunctionBase() {}

    // number of unknowns
    inline virtual int n_unknowns() override { return 4; }

    /** evaluates the spring element's energy
     * \param _x contains x_a and x_b contiguously,
     *           i.e. _x = [x_a, x_b], i.e. _x is of dimension 4
     * \param _coeffs stores the constant k,
     *                i.e. _coeffs[0] = k */
    inline virtual double eval_f(const Vec &_x, const Vec &_coeffs) override {
      //------------------------------------------------------//
      // Todo: implement the function f(x) = 1/2 * k * ((x[0] - x[2])^2 + (x[1]
      // - x[3])^2)
      return 0.5 * _coeffs[0] * pow(_x[0] - _x[3], 2) + pow(_x[1] - _x[3], 2)

      //------------------------------------------------------//
    }

    /** evaluates the spring element's energy gradient
     * \param _x contains x_a and x_b contiguously,
     *           i.e. _x = [x_a, x_b], i.e. _x is of dimension 4
     * \param _coeffs stores the constant k,
     *                i.e. _coeffs[0] = k
     * \param _g the output gradient, which should also be of dimension 4 */
    inline virtual void eval_gradient(const Vec &_x, const Vec &_coeffs,
                                      Vec &_g) override {
      //------------------------------------------------------//
      // Todo: implement the gradient and store in _g

      //------------------------------------------------------//
    }

    /** evaluates the spring element's energy Hessian
     * \param _x contains x_a and x_b contiguously,
     *           i.e. _x = [x_a, x_b], i.e. _x is of dimension 4
     * \param _coeffs stores the constant k,
     *                i.e. _coeffs[0] = k
     * \param _H the output Hessian, which should be a 4x4 Matrix */
    inline virtual void eval_hessian(const Vec &_x, const Vec &_coeffs,
                                     Mat &_H) override {
      //------------------------------------------------------//
      // Todo: implement the hessian matrix and store in _H

      //------------------------------------------------------//
    }
  };

  //=============================================================================
}
