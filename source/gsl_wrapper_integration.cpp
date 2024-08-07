//
// Created by M. J. Steil on 2017.01.19.
//

#include "../gsl_wrapper.hpp"

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wreturn-stack-address"

namespace gsl_wrapper {
    using namespace std;

    gsl_integration::gsl_integration()= default;

    gsl_integration::gsl_integration(gsl_function *fkt, const vector<double>& range, double abs_error_threshold, double rel_error_threshold, const string& method, int print) {
        set_parameters(fkt,abs_error_threshold,rel_error_threshold);
        comp(method,range,print);
        free();
    }

    gsl_integration::gsl_integration(gsl_function *fkt, const vector<double>& range, double abs_error_threshold, double rel_error_threshold, const string& method,size_t size, int key, int print) {
        set_parameters(fkt,abs_error_threshold,rel_error_threshold, size, key);
        comp(method,range,print);
        free();
    }

    void gsl_integration::set_parameters(gsl_function *fkt, double abs_error_threshold_in, double rel_error_threshold_in, size_t size_in, int key_in) {
        size = size_in;
        key = key_in;

        rel_error_threshold = rel_error_threshold_in;
        abs_error_threshold = abs_error_threshold_in;
        fkt_pp = fkt;

        comp_Int = gsl_integration_workspace_alloc (this->size);
    }

    void gsl_integration::comp(const string& method,const vector<double>& range, int print) {
        if(method=="qng"){
            comp_qng(range, print);
        }else if(method=="qag"){
            comp_qag(range, print);
        }else if(method=="qags"){
            comp_qags(range, print);
        }else if(method=="qagi"){
            comp_qagi(range, print);
        }else if(method=="qagiu"){
            comp_qagiu(range, print);
        }else if(method=="qagil"){
            comp_qagil(range, print);
        }else if(method=="odeiv2"){
            comp_odeiv2(range, print);
        }else if(method=="trap"){
            comp_trap(range, print);
        }else{
            GSL_ERROR_VAL ("gsl_integration::comp: Invalid method!", GSL_FAILURE,);
        }
    }

    void gsl_integration::comp_qng(vector<double> range, int print) {
        if(print){
            printf("|--: gsl_integration::comp_qng: with size %zu \n",size);
        }

        size_t neval;
        gsl_integration_qng(fkt_pp, range[0], range[1], abs_error_threshold, rel_error_threshold, &F, &F_error, &neval);

        F_rel_error=F_error/F;

        if(print){
            printf("|--: gsl_integration::comp_qng: F=%.8E (+-%.4E(=%.4E*F)) with %zu intervals\n",F,F_error,F_rel_error,neval);
        }
    }

    void gsl_integration::comp_qag(vector<double> range, int print) {
        if(print){
            printf("|--: gsl_integration::comp_qag: with size %zu and key=%d \n",size,key);
        }

        gsl_integration_qag(fkt_pp, range[0], range[1], abs_error_threshold, rel_error_threshold, size, key, comp_Int, &F, &F_error);

        F_rel_error=F_error/F;

        if(print){
            printf("|--: gsl_integration::comp_qag: F=%.8E (+-%.4E(=%.4E*F)) in %zu intervals\n",F,F_error,F_rel_error,comp_Int->size);
        }
    }

    void gsl_integration::comp_qags(vector<double> range, int print) {
        if(print){
            printf("|--: gsl_integration::comp_qags: with size %zu \n",size);
        }

        gsl_integration_qags(fkt_pp, range[0], range[1], abs_error_threshold, rel_error_threshold, size, comp_Int, &F, &F_error);

        F_rel_error=F_error/F;

        if(print){
            printf("|--: gsl_integration::comp_qags: F=%.8E (+-%.4E(=%.4E*F)) in %zu intervals\n",F,F_error,F_rel_error,comp_Int->size);
        }
    }

    void gsl_integration::comp_qagi(const vector<double>&, int print) {
        if(print){
            printf("|--: gsl_integration::comp_qagi: with size %zu \n",size);
        }

        gsl_integration_qagi(fkt_pp, abs_error_threshold, rel_error_threshold, size, comp_Int, &F, &F_error);

        F_rel_error=F_error/F;

        if(print){
            printf("|--: gsl_integration::comp_qagi: F=%.8E (+-%.4E(=%.4E*F)) in %zu intervals\n",F,F_error,F_rel_error,comp_Int->size);
        }
    }

    void gsl_integration::comp_qagiu(vector<double> range, int print) {
        if(print){
            printf("|--: gsl_integration::comp_qagiu: with size %zu \n",size);
        }

        gsl_integration_qagiu(fkt_pp, range[0], abs_error_threshold, rel_error_threshold, size, comp_Int, &F, &F_error);

        F_rel_error=F_error/F;

        if(print){
            printf("|--: gsl_integration::comp_qagiu: F=%.8E (+-%.4E(=%.4E*F)) in %zu intervals\n",F,F_error,F_rel_error,comp_Int->size);
        }
    }

    void gsl_integration::comp_qagil(vector<double> range, int print) {
        if(print){
            printf("|--: gsl_integration::comp_qagil: with size %zu \n",size);
        }

        gsl_integration_qagil(fkt_pp, range[1], abs_error_threshold, rel_error_threshold, size, comp_Int, &F, &F_error);

        F_rel_error=F_error/F;

        if(print){
            printf("|--: gsl_integration::comp_qagil: F=%.8E (+-%.4E(=%.4E*F)) in %zu intervals\n",F,F_error,F_rel_error,comp_Int->size);
        }
    }



    void gsl_integration::comp_odeiv2(vector<double> range, int print) {
        if(print){
            printf("|--: gsl_integration::comp_odeiv2: with size %zu \n",size);
        }

        function<double (double, const double *, double *, int)> dQdx_eq = [&](double x , const double * , double *, int i)->double{
            double df_i_dr;

            if(i==0) {
                df_i_dr = GSL_FN_EVAL(fkt_pp,x) ;
            } else {
                GSL_ERROR_VAL ("gsl_integration::comp_odeiv2:dQdx_eq", GSL_FAILURE,GSL_FAILURE);
            }
            return df_i_dr;
        };

        gsl_odeiv2_ode_system Q_eq(1);
        Q_eq.df= &dQdx_eq;
        struct gsl_odeiv2_ode_struct comp_Q_ode = {&Q_eq};

        const double comp_h0 = range[0];   // [km], Initial step size
        const double comp_hmax = (range[1]-range[0])*1E-3;   // [km], Maximal step size
        const double comp_xMax = range[1];    // [km], Maximum value for radial integration
        double comp_odeiv2_control_params[4] = {rel_error_threshold,abs_error_threshold,1,1};  // [error_rel, error_abs, error_ay, error_day]
        double comp_error_scale[1]           = {1};       // Error scale for *comp_l0_control

        auto comp_nMax = (unsigned long)size;   // Maximum number of Iteration steps of gsl_odeiv2_evolve_apply

        gsl_odeiv2_solver comp_solver(1,0,print);
        comp_solver.set_parameters(comp_h0,
                                      comp_hmax,
                                      comp_xMax,
                                      comp_nMax,
                                      (double *)comp_odeiv2_control_params,
                                      (double *)comp_error_scale
        );

        switch (key){
            case 0: // Explicit embedded Runge-Kutta (2, 3) method.
                comp_solver.set_system(&comp_Q_ode,gsl_odeiv2_step_rk2);
                break;
            case 1: // Explicit 4th order (classical) Runge-Kutta. Error estimation is carried out by the step doubling method.
                comp_solver.set_system(&comp_Q_ode,gsl_odeiv2_step_rk4);
                break;
            case 2: // Explicit embedded Runge-Kutta-Fehlberg (4, 5) method. This method is a good general-purpose integrator.
                comp_solver.set_system(&comp_Q_ode,gsl_odeiv2_step_rkf45);
                break;
            case 3: // Explicit embedded Runge-Kutta Cash-Karp (4, 5) method.
                comp_solver.set_system(&comp_Q_ode,gsl_odeiv2_step_rkck);
                break;
            default: // Explicit embedded Runge-Kutta Prince-Dormand (8, 9) method. [DEFAULT]
                comp_solver.set_system(&comp_Q_ode,gsl_odeiv2_step_rk8pd);
        }

        double Q[1] = {0};
        double dQ[1] = {GSL_FN_EVAL(fkt_pp,range[0])};
        double comp_x = range[0];
        double comp_h = (range[1]-range[0])*1E-3;

        if(print){
            printf("|--: gsl_integration::comp_odeiv2: dQ(range[0])=%.8E, dQ(range[1])=%.8E \n",GSL_FN_EVAL(fkt_pp,range[0]),GSL_FN_EVAL(fkt_pp,range[1]));
        }

        unsigned int comp_n=0;

        //region Main comp_l0_solver.comp loop
        comp_solver.set_IC(comp_n,
                              comp_x,
                              comp_h,
                              Q,
                              dQ,
                                nullptr);// Inital conditions for integration
        comp_solver.comp();// Integration
        comp_solver.get_nrh(comp_n,comp_x,comp_h);
        //endregion

        F=Q[0];
        F_error = comp_solver.ferr_max[0];
        F_rel_error = F_error/F;

        if(print){
            printf("|--: gsl_integration::comp_odeiv2: F=%.8E (+-%.4E(=%.4E*F)) in %d intervals\n",F,F_error,F_rel_error,comp_n);
        }
    }

    void gsl_integration::comp_trap(vector<double> range, int print) {
        unsigned long length = range.size()-1;

        if(print){
            printf("|--: gsl_integration::comp_trap: %zu intervals with subdivision %d \n",length,key);
        }

        double a,b,c,d,tij,tij2,f0,f1,f2,F0;
        F0=0;
        F=0;

        for(int i=0; i<length-1; i++){
            a = range[i];
            d = (range[i+1]-range[i])/key;

            for(int j=1;j<=key;j++){
                b = a+d;
                c = a+d/2;

                f0 = GSL_FN_EVAL(fkt_pp,a);
                f1 = GSL_FN_EVAL(fkt_pp,c);
                f2 = GSL_FN_EVAL(fkt_pp,b);

                tij = 0.25*(f0+2*f1+f2)*d;
                tij2 = 0.5*(f0+f2)*d;

                F   += tij;
                F0  += tij2;

                a = b;
            }
        }

        F_rel_error = F0/F-1.;

        if(print){
            printf("|--: gsl_integration::comp_trap: F=%.8E (%.4E*F)\n",F,F_rel_error);
        }
    }

    void gsl_integration::free() const {
        gsl_integration_workspace_free (comp_Int);
    }



}
#pragma clang diagnostic pop