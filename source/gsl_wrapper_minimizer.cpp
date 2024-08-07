//
// Created by M. J. Steil on 2017.01.19.
//

#include "../gsl_wrapper.hpp"

namespace gsl_wrapper {
    using namespace std;

    gsl_minimizer::gsl_minimizer()= default;

    gsl_minimizer::gsl_minimizer(gsl_function *fkt, int nMax, double rel_error_threshold_x, double rel_error_threshold_f, vector<double> x_range,int print){
        set_parameters(fkt,nMax,rel_error_threshold_x, rel_error_threshold_f,x_range);
        comp(print);
    }

    gsl_minimizer::gsl_minimizer(gsl_function *fkt, int nMax, double rel_error_threshold_x, double rel_error_threshold_f, vector<double> x_range,const gsl_min_fminimizer_type *type,int print){
        set_parameters(fkt,nMax,rel_error_threshold_x, rel_error_threshold_f,x_range,type);
        comp(print);
    }

    void gsl_minimizer::set_parameters(gsl_function *fkt, int nMax, double rel_error_threshold_x, double rel_error_threshold_f, vector<double> x_range, const gsl_min_fminimizer_type *type) {
        this->fkt = fkt;

        this->nMax = nMax;
        this->rel_error_threshold_x = rel_error_threshold_x;
        this->rel_error_threshold_f = rel_error_threshold_f;
        this->x_lower = x_range[0];
        this->x_min   = x_range[1];
        this->x_upper = x_range[2];

        comp_type = type;
        comp_minimizer = gsl_min_fminimizer_alloc (comp_type);
    }

    void gsl_minimizer::comp(int print) {
        //region Print gsl_minimizer::comp header
        if(print){
            printf("|--: gsl_minimizer searching minimum with %s method\n",gsl_min_fminimizer_name(comp_minimizer));
            printf("|--: in x_range={%.3E,%.3E,%.3E} -> f_range={%.3E,%.3E,%.3E}\nn",
                   x_lower,x_min,x_upper,
                   GSL_FN_EVAL(fkt,x_lower),GSL_FN_EVAL(fkt,x_min),GSL_FN_EVAL(fkt,x_upper));
            printf("|---: n\t \t x_min \t \t x_rel_error \t \t f_min \t \t f_rel_error\n");
        }
        //endregion

        gsl_min_fminimizer_set (comp_minimizer, fkt, x_min, x_lower, x_upper);

        // Main gsl_min_fminimizer_iterate loop
        while (status == GSL_CONTINUE && n < nMax){
            n++;
            x_sum_prev=x_lower+x_upper+x_upper;
            gsl_min_fminimizer_iterate (comp_minimizer);

            x_min   = gsl_min_fminimizer_x_minimum (comp_minimizer);
            x_lower = gsl_min_fminimizer_x_lower (comp_minimizer);
            x_upper = gsl_min_fminimizer_x_upper (comp_minimizer);

            x_sum=x_lower+x_upper+x_upper;
            x_rel_error = abs((x_upper-x_lower)/(x_lower+x_upper)*2);

            f_min = gsl_min_fminimizer_f_minimum(comp_minimizer);
            f_lower = gsl_min_fminimizer_f_lower(comp_minimizer);
            f_upper = gsl_min_fminimizer_f_upper(comp_minimizer);

            f_rel_error = abs((f_upper-f_lower)/(f_lower+f_upper)*2);

            //region x_sum==x_sum_prev break
            if(x_sum==x_sum_prev){
                x_sum_equal++;
                if(x_sum_equal>x_sum_equal_lim){
                    status = GSL_SUCCESS;
                    if(print){printf("|--: gsl_minimizer: x_sum did not change in %d consecutive steps!\n",x_sum_equal_lim);}
                }
            }else{
                x_sum_equal = 0;
            }
            //endregion

            //region x_rel_error<rel_error_threshold_x||f_rel_error<rel_error_threshold_f break
            if(x_rel_error<rel_error_threshold_x||f_rel_error<rel_error_threshold_f){
                status = GSL_SUCCESS;
            }
            //endregion

            if(print&&n%5-1<0){printf("|---: %.d \t %.4E \t %.4E \t %.4E \t %.4E\n",n,x_min,x_rel_error,f_min,f_rel_error);}
        }

        //region Print result
        if(print) {
            if (n == nMax || status != GSL_SUCCESS) {
                printf("|--: gsl_minimizer GSL_ERROR: %s | Stopped at Step: n=%.d - x_min=%.8E (rel_err:~%.8E) - f_min=%.8E (rel_err:~%.8E)\n", gsl_strerror(status), n, x_min, x_rel_error,f_min, f_rel_error);
            } else {
                printf("|--: gsl_minimizer %s | Stopped at Step: n=%.d - x_min=%.8E (rel_err:~%.8E) - f_min=%.8E (rel_err:~%.8E) \n", gsl_strerror(status), n, x_min, x_rel_error,f_min, f_rel_error);
            }
        }
        //endregion

        gsl_min_fminimizer_free (comp_minimizer);
    }

}