//
// Created by M. J. Steil on 2017.02.15.
//

#include "../gsl_wrapper.hpp"

namespace gsl_wrapper {
    using namespace std;

    gsl_hermite_spline_obj::gsl_hermite_spline_obj() = default;

    gsl_hermite_spline_obj::gsl_hermite_spline_obj(vector<double> ri, vector<double> fi, vector<double> dfi, size_t length_in) {
        // Copy data
        length = length_in;
        r_data.assign(ri.begin(),ri.begin()+length);
        f_data.assign(fi.begin(),fi.begin()+length);
        df_data.assign(dfi.begin(),dfi.begin()+length);

        // Setup f hermite spline
        f_acc = gsl_interp_accel_alloc();

        // Setup df cspline
        df_r = gsl_spline_alloc(gsl_interp_cspline, length);
        gsl_spline_init(df_r, r_data.data(), df_data.data(), length);
        df_acc = gsl_interp_accel_alloc();
    }

    double gsl_hermite_spline_obj::interpol_herm(double r) {
        if(r!=fr.first) {
            fr.first = r;

            // Interpolation index lookup - gsl_interp_accel_find
            size_t i = gsl_interp_accel_find(f_acc, r_data.data(), length, r);
            size_t ip1 = i + 1;

            // Cubic hermite interpolation polynomial
            double dx = r_data[ip1] - r_data[i];
            double u = (r - r_data[i]) / dx;
            double u2 = u * u;
            double u3 = u2 * u;

            fr.second = f_data[i] * (2. * u3 - 3. * u2 + 1.)
                        + f_data[ip1] * (3. * u2 - 2. * u3)
                        + df_data[i] * dx * (u3 - 2. * u2 + u)
                        - df_data[ip1] * dx * (u2 - u3);
        }
        return fr.second;
    }

    double gsl_hermite_spline_obj::f(double r) {
        return interpol_herm(r);
    }

    double gsl_hermite_spline_obj::df(double r) {
        return gsl_spline_eval (df_r,r,df_acc);
    }

    void gsl_hermite_spline_obj::free() {
        // Clear internal data
        length = 0;
        r_data = vector<double>();
        f_data = vector<double>();
        df_data = vector<double>();

        // Free Splines and accelerators
        gsl_spline_free (df_r);
        gsl_interp_accel_free (f_acc);
        gsl_interp_accel_free (df_acc);
    }

}