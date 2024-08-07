//
// Wrapper class for gsl_interp
// Created by M. J. Steil on 2017.01.19.
//

#include "../gsl_wrapper.hpp"

namespace gsl_wrapper {
    using namespace std;

    gsl_spline_obj::gsl_spline_obj() = default;

    gsl_spline_obj::gsl_spline_obj(vector<double> ri, vector<double> fi,size_t length) {

        comp_type = gsl_interp_cspline;
        f_r = gsl_spline_alloc(comp_type, length);
        f_acc = gsl_interp_accel_alloc();

        gsl_spline_init(f_r, ri.data(), fi.data(), length);
    }

    gsl_spline_obj::gsl_spline_obj(vector<double> ri, vector<double> fi,size_t length, const gsl_interp_type *type) {

        comp_type = type;
        f_r = gsl_spline_alloc(comp_type, length);
        f_acc = gsl_interp_accel_alloc();

        gsl_spline_init(f_r, ri.data(), fi.data(), length);
    }

    double gsl_spline_obj::f(double r) {
        return gsl_spline_eval (f_r,r,f_acc);
    }

    double gsl_spline_obj::df(double r) {
        return gsl_spline_eval_deriv (f_r,r,f_acc);
    }

    double gsl_spline_obj::ddf(double r) {
        return gsl_spline_eval_deriv2 (f_r,r,f_acc);
    }

    void gsl_spline_obj::free(){
        gsl_spline_free (f_r);
        gsl_interp_accel_free (f_acc);
    }

}

