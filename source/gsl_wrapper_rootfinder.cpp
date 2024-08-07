//
// Created by M. J. Steil on 2017.01.19.
//

#include "../gsl_wrapper.hpp"

namespace gsl_wrapper {
    using namespace std;

    gsl_rootfinder::gsl_rootfinder()= default;

    gsl_rootfinder::gsl_rootfinder(gsl_function *fkt, int nMax, double rel_error_threshold_x, double rel_error_threshold_f, vector<double> x_range,int print){
        set_parameters(fkt,nMax,rel_error_threshold_x, rel_error_threshold_f,x_range);
        comp(print);
    }

    gsl_rootfinder::gsl_rootfinder(gsl_function *fkt, int nMax, double rel_error_threshold_x, double rel_error_threshold_f, vector<double> x_range, const gsl_root_fsolver_type *type, int print){
        set_parameters(fkt,nMax,rel_error_threshold_x, rel_error_threshold_f,x_range,type);
        comp(print);
    }

    void gsl_rootfinder::set_parameters(gsl_function *fkt, int nMax, double rel_error_threshold_x, double rel_error_threshold_f, vector<double> x_range, const gsl_root_fsolver_type *type) {
        this->fkt = fkt;

        this->nMax = nMax;
        this->rel_error_threshold_x = rel_error_threshold_x;
        this->x_lower = x_range[0];
        this->x_min   = x_range[1];
        this->x_upper = x_range[2];

        comp_type = type;
        comp_fsolver = gsl_root_fsolver_alloc (comp_type);
    }

    void gsl_rootfinder::comp(int print) {
        //region Print gsl_rootfinder::comp_qag header
        if(print){
            printf("|--: gsl_rootfinder searching root with %s method\n", gsl_root_fsolver_name (comp_fsolver));
            printf("|--: in x_range={%.3E,-,%.3E} -> f_range={%.3E,-,%.3E}\n",
                   x_lower,x_upper,
                   GSL_FN_EVAL(fkt,x_lower),GSL_FN_EVAL(fkt,x_upper));
            printf("|---: n\t \t x_min \t \t x_rel_error\n");
        }
        //endregion

        gsl_root_fsolver_set (comp_fsolver, fkt, x_lower, x_upper);

        // Main gsl_root_fsolver_iterate loop
        while (status == GSL_CONTINUE && n < nMax){
            n++;
            x_sum_prev=x_lower+x_upper+x_upper;

            gsl_root_fsolver_iterate (comp_fsolver);

            x_min   = gsl_root_fsolver_root (comp_fsolver);
            x_lower = gsl_root_fsolver_x_lower (comp_fsolver);
            x_upper = gsl_root_fsolver_x_upper (comp_fsolver);

            x_sum=x_lower+x_upper+x_upper;
            x_rel_error = abs((x_upper-x_lower)/(x_lower+x_upper)*2);

            //region x_sum==x_sum_prev break
            if(x_sum==x_sum_prev){
                x_sum_equal++;
                if(x_sum_equal>x_sum_equal_lim){
                    status = GSL_SUCCESS;
                    if(print){printf("|--: gsl_rootfinder: x_sum did not change in %d consecutive steps!\n",x_sum_equal_lim);}
                }
            }else{
                x_sum_equal = 0;
            }
            //endregion

            //region x_rel_error<rel_error_threshold_x break
            if(x_rel_error<rel_error_threshold_x){
                status = GSL_SUCCESS;
            }
            //endregion

            if(print&&n%5-1<0){printf("|---: %.d \t %.4E \t %.4E \n",n,x_min,x_rel_error);}
        }

        //region Print result
        if(print) {
            if (n == nMax || status != GSL_SUCCESS) {
                printf("|--: gsl_rootfinder GSL_ERROR: %s | Stopped at Step: n=%.d - x_min=%.8E (rel_err:~%.8E)\n", gsl_strerror(status), n, x_min, x_rel_error);
            } else {
                printf("|--: gsl_rootfinder %s | Stopped at Step: n=%.d - x_min=%.8E (rel_err:~%.8E)\n", gsl_strerror(status), n, x_min, x_rel_error);
            }
        }
        //endregion

        gsl_root_fsolver_free (comp_fsolver);
    }

}