/***
 * gsl_wrapper
 * Created by M. J. Steil on 2017.01.19.
 */


#ifndef GSL_WRAPPER_HPP
#define GSL_WRAPPER_HPP

#include <iostream>
#include <utility>
#include <vector>
#include <functional>
#include <algorithm>
#include <sstream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>

namespace gsl_wrapper {
    using namespace std;

    template< class F > class gsl_function_pp : public gsl_function {
    public:
        /*!
         * gsl_function wrapper
         * @param func
         * @usage
         *  auto ptr = [=](double r)->double{    return 42;  };<br>
         *  gsl_function_pp &lt decltype(ptr)&gt Fp(ptr);<br>
         *  gsl_function *F = static_cast<gsl_function*>(&Fp);
         * @source
         * <a href="http://stackoverflow.com/a/18413206">http://stackoverflow.com/a/18413206, 160726 23:14</a>
         */
        explicit gsl_function_pp(const F& func) : gsl_function_struct(), _func(func) {
            function = &gsl_function_pp::invoke;
            params=this;
        }
    private:
        const F& _func;
        static double invoke(double x, void *params) {
            return static_cast<gsl_function_pp*>(params)->_func(x);
        }

    };

    class gsl_function_pp_fkt : public gsl_function {
    public:
        /*!
         * gsl_function wrapper
         * @param func
         * @usage
         *  gsl_function_pp_fkt Fp([=](double r)->double{return  42;});<br>
         * @source
         * <a href="http://stackoverflow.com/a/18413206">http://stackoverflow.com/a/18413206, 160726 23:14</a>
         */
        explicit gsl_function_pp_fkt(std::function<double(double)> func);
        gsl_function * gsl_fkt();
    private:
        std::function<double(double)> _func;
        static double invoke(double x, void *params);
    };

    // gsl_odeiv2 ODE Wrapper

    class gsl_odeiv2_ode_system{
    public:
        /*!
         * gsl_odeiv2_ode_system for gsl_odeiv2_ode_struct setup
         * @param dim Dimesion of ODE system
         * @usage
         *     gsl_odeiv2_ode_system system(dim);<br>
         *     system.df= &equation; //equation: Lambda function ode <br>
         *     struct gsl_odeiv2_ode_struct system_struct  = {&system};
         */
        explicit gsl_odeiv2_ode_system(int dim);
        int dim;
        function<double (double, const double *, double *, int)> *df;
    };

    struct gsl_odeiv2_ode_struct {gsl_odeiv2_ode_system *ode;};

    class gsl_odeiv2_data{
    public:
        gsl_odeiv2_data(int dim, unsigned long nMax);
        gsl_odeiv2_data();

        int mode, data_dim;
        unsigned long  n ,data_nMax, data_vec_size;
        vector<vector<double>> data;

        void put(unsigned long n, vector<double> rhfi);
        void insert(int data_i,vector<double> *in,int pos=0);

        void scale(double s, int l); // scale vector l with s
        void offset(double a, int l); // offset vector l with a
        void add(int l,int m); // add vectors l and m into vector l

        void add_fkt(int l, int m, function<double(double)> fkt);

        void copy(int data_i,vector<double>* target);
    };

    class gsl_odeiv2_solver{
    public:
        explicit gsl_odeiv2_solver(int dim, int mode = 1,int print = 1);

        int mode, print, solver_status,solver_comp_status;
        size_t  dim;

        double hMin, hMax; // minimal and maximal step size
        double rMax;
        unsigned long nMax;

        double *odeiv2_control_params; // [error_rel, error_abs, error_ay, error_day]
        double *error_scale;

        const gsl_odeiv2_step_type *odeiv2_type;
        gsl_odeiv2_step *odeiv2_step;
        gsl_odeiv2_control *odeiv2_control;
        gsl_odeiv2_evolve *odeiv2_evolve;
        gsl_odeiv2_system odeiv2_system;
        void set_system(struct gsl_odeiv2_ode_struct *, const gsl_odeiv2_step_type *type = gsl_odeiv2_step_rk8pd);

        unsigned int n;

        double r,h;
        double *f, *df;
        gsl_odeiv2_data *data;
        vector<double> rhfi;
        vector<double> ferr_max;

        void set_parameters(double hMin, double hMax, double rMax, unsigned long nMax, double *odeiv2_control_params ,double *error_scale);
        void set_IC(unsigned int &n,double &r, double &h, double *f, double *df, gsl_odeiv2_data *);

        void comp_step();
        void comp();
        void get_nrh(unsigned int &n_in,double &r_in, double &h_in) const;

    };

    int gsl_odeiv2_ode(double, const double *, double *, void *);

    // Interpolation objects

    GSL_VAR const gsl_interp_type * gsl_interp_poly2;

    class gsl_spline_obj{
    private:
        gsl_interp_accel *f_acc;
        gsl_spline *f_r;
        const gsl_interp_type *comp_type;

    public:
        gsl_spline_obj(vector<double> ri ,vector<double> fi, size_t length);
        gsl_spline_obj(vector<double> ri ,vector<double> fi, size_t length, const gsl_interp_type *type);
        gsl_spline_obj();

        double f(double);
        double df(double);
        double ddf(double);

        void free();
    };

    class gsl_hermite_spline_obj {
    private:
        vector<double> r_data, f_data, df_data;
        size_t length;

        gsl_interp_accel *f_acc;
        gsl_interp_accel *df_acc;

        gsl_spline *df_r;
        pair<double, double> fr = make_pair(0, 0);

        double interpol_herm(double r);

    public:
        gsl_hermite_spline_obj(vector<double> ri, vector<double> fi, vector<double> dfi, size_t length_in);

        gsl_hermite_spline_obj();

        double f(double r);

        double df(double r);

        void free();
    };

    // Minimizer and rootfinder

    class gsl_minimizer{
    public:
        gsl_minimizer();
        gsl_minimizer(gsl_function *fkt, int nMax, double abs_error_threshold, double rel_error_threshold, vector<double> x_range, int print = 0);
        gsl_minimizer(gsl_function *fkt, int nMax, double abs_error_threshold, double rel_error_threshold, vector<double> x_range, const gsl_min_fminimizer_type *type, int print = 0);

        gsl_function *fkt{};
        double rel_error_threshold_x{},rel_error_threshold_f{};
        int nMax{};
        double x_lower{}, x_min{}, x_upper{}, x_rel_error{};
        double x_sum{}, x_sum_prev{};
        double f_lower{}, f_min{}, f_upper{}, f_rel_error{};

        int n = 0;
        int x_sum_equal = 0;
        int x_sum_equal_lim = 2;
        int status = GSL_CONTINUE;

        const gsl_min_fminimizer_type *comp_type{};
        gsl_min_fminimizer *comp_minimizer{};


        void set_parameters(gsl_function *fkt, int nMax, double abs_error_threshold, double rel_error_threshold, vector<double> x_range, const gsl_min_fminimizer_type *type = gsl_min_fminimizer_brent);
        void comp(int print = 0);
    };

    class gsl_rootfinder{
    public:
        gsl_rootfinder();
        gsl_rootfinder(gsl_function *fkt, int nMax, double abs_error_threshold, double rel_error_threshold, vector<double> x_range, int print = 0);
        gsl_rootfinder(gsl_function *fkt, int nMax, double abs_error_threshold, double rel_error_threshold, vector<double> x_range, const gsl_root_fsolver_type *type, int print = 0);

        gsl_function *fkt;
        double rel_error_threshold_x;
        int nMax;
        double x_lower, x_min, x_upper, x_rel_error;
        double x_sum, x_sum_prev;

        int n = 0;
        int x_sum_equal = 0;
        int x_sum_equal_lim = 2;
        int status = GSL_CONTINUE;

        const gsl_root_fsolver_type *comp_type;
        gsl_root_fsolver *comp_fsolver;


        void set_parameters(gsl_function *fkt, int nMax, double abs_error_threshold, double rel_error_threshold, vector<double> x_range,const gsl_root_fsolver_type *type = gsl_root_fsolver_brent);
        void comp(int print = 0);
    };

    // Integrator

    class gsl_integration{
    public:
        gsl_integration();
        gsl_integration(gsl_function *fkt, const vector<double>& range, double abs_error_threshold, double rel_error_threshold, const string& method = "qag", int print = 0);
        gsl_integration(gsl_function *fkt, const vector<double>& range, double abs_error_threshold, double rel_error_threshold, const string& method = "qag", size_t size_in = 100, int key_in = 6, int print = 0);

        double rel_error_threshold,abs_error_threshold;
        size_t size;    // Allocated size of  gsl_integration_workspace * comp_Int
        int key;        // Key of gsl_integration_qag of comp_qag()

        double F, F_error, F_rel_error;
        gsl_function *fkt_pp;
        gsl_integration_workspace * comp_Int;

        void set_parameters(gsl_function *fkt, double abs_error_threshold_in, double rel_error_threshold_in, size_t size_in = 100, int key_in = 6);
        void comp(const string& method,const vector<double>& range, int print = 0);

        void comp_qng(vector<double> range, int print = 0); // QUADPACK non-adaptive Gauss-Kronrod

        void comp_qag(vector<double> range, int print = 0); // QUADPACK adaptive Gauss-Kronrod
        void comp_qags(vector<double> range, int print = 0); // QUADPACK adaptive Gauss-Kronrod with singularities

        void comp_qagi(const vector<double>& range, int print = 0);  // QUADPACK adaptive Gauss-Kronrod on infinite interval (-inf,+inf)
        void comp_qagiu(vector<double> range, int print = 0); // QUADPACK adaptive Gauss-Kronrod on semi-infinite interval (a,+inf)
        void comp_qagil(vector<double> range, int print = 0); // QUADPACK adaptive Gauss-Kronrod on semi-infinite interval (-inf,b)

        void comp_odeiv2(vector<double> range, int print = 0); // ODEIV2 integrator

        void comp_trap(vector<double> range, int print = 0);    // Trapezoidal method with step doubling for error estimate

        void free() const;

    };
}
#endif //GSL_WRAPPER_HPP