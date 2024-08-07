//
// Created by M. J. Steil on 2017.01.19.
//

#include "../gsl_wrapper.hpp"

namespace gsl_wrapper {

    // gsl_odeiv2 System
    int gsl_odeiv2_ode(double r, const double *f, double *df, void *p) {
        auto * params = (struct gsl_odeiv2_ode_struct *)p;
        gsl_odeiv2_ode_system *system = (params->ode);

        function<double (double, const double *, double *, int)> system_df=*(system->df);

        for(int i=0;i<system->dim;i++){
            df[i] = system_df(r,f,df,i);
        }

        return GSL_SUCCESS;
    }

    gsl_odeiv2_ode_system::gsl_odeiv2_ode_system(int dim) {
        this->dim=dim;
    }


    // gsl_odeiv2_solver
    gsl_odeiv2_solver::gsl_odeiv2_solver(int dim,int mode, int print) {
        this->dim=(size_t)dim;
        this->mode=mode;
        this->print=print;

        ferr_max = vector<double>(this->dim);

        solver_status=1;
    }

    void gsl_odeiv2_solver::set_parameters(double hMin, double hMax, double rMax, unsigned long nMax, double *odeiv2_control_params ,double *error_scale) {
        if(solver_status>0){
            this->hMin=hMin;
            this->hMax=hMax;
            this->rMax=rMax;
            this->nMax=nMax;
            this->odeiv2_control_params=odeiv2_control_params;
            this->error_scale=error_scale;
            solver_status=2;
        }else{
            GSL_ERROR_VAL ("gsl_odeiv2_solver::set_parameters: solve_status>0 required.", GSL_FAILURE,);
        }
    }

    void gsl_odeiv2_solver::set_system(struct gsl_odeiv2_ode_struct *ode_struct, const gsl_odeiv2_step_type *type) {
        if(solver_status>1){
            /* Compatible step types are
                gsl_odeiv2_step_rk2 - Explicit embedded Runge-Kutta (2, 3) method.
                gsl_odeiv2_step_rk4 - Explicit 4th order (classical) Runge-Kutta. Error estimation is carried out by the step doubling method.
                gsl_odeiv2_step_rkf45 - Explicit embedded Runge-Kutta-Fehlberg (4, 5) method. This method is a good general-purpose integrator.
                gsl_odeiv2_step_rkck - Explicit embedded Runge-Kutta Cash-Karp (4, 5) method.
                gsl_odeiv2_step_rk8pd - Explicit embedded Runge-Kutta Prince-Dormand (8, 9) method. [DEFAULT]
             */
            odeiv2_type = type;
            odeiv2_step = gsl_odeiv2_step_alloc(odeiv2_type, dim);
            odeiv2_control = gsl_odeiv2_control_scaled_new(odeiv2_control_params[0], // error_rel
                                                           odeiv2_control_params[1], // error_abs
                                                           odeiv2_control_params[2], // error_ay
                                                           odeiv2_control_params[3], // error_day
                                                           error_scale,
                                                           dim); // Adaptive Step-size Control
            odeiv2_evolve = gsl_odeiv2_evolve_alloc(dim);
            odeiv2_system = {gsl_odeiv2_ode, nullptr, dim,ode_struct};
            solver_status=3;

            if(print){
                printf("|--: gsl_odeiv2_solver::set_system: System setup with %s stepper.\n",gsl_odeiv2_step_name(odeiv2_step));
            }
        }else{
            GSL_ERROR_VAL ("gsl_odeiv2_solver::set_system: solve_status>1 required.", GSL_FAILURE,);
        }
    }

    void gsl_odeiv2_solver::set_IC(unsigned int &n_in,double &r_in, double &h_in , double *f_in, double *df_in, gsl_odeiv2_data *data_in){
        if(solver_status>2){
            n=n_in;
            r=r_in;
            h=h_in;

            f=f_in;
            df=df_in;
            data=data_in;
            odeiv2_evolve->dydt_in = df;

            solver_comp_status = GSL_SUCCESS;
            solver_status=4;
        }else{
            GSL_ERROR_VAL ("gsl_odeiv2_solver::set_IC: solve_status>2 required.", GSL_FAILURE,);
        }
    }

    void gsl_odeiv2_solver::comp_step(){
        if(solver_status>3){
            if (r<rMax) {
                // Step
                if(h>hMax){
                    h=hMax;
                    solver_comp_status = gsl_odeiv2_evolve_apply(odeiv2_evolve, odeiv2_control, odeiv2_step, &odeiv2_system, &r,rMax, &h, f);
                    df=odeiv2_evolve->dydt_out;
                }else{
                    solver_comp_status = gsl_odeiv2_evolve_apply(odeiv2_evolve, odeiv2_control, odeiv2_step, &odeiv2_system, &r,rMax, &h, f);
                    df=odeiv2_evolve->dydt_out;
                }
                n++;

                // Adjust ferr_max = maximal error during integration
                for(int j=0;j<dim;j++){
                    if(odeiv2_evolve->yerr[j]>ferr_max[j]){
                        ferr_max[j] =odeiv2_evolve->yerr[j];
                    }
                }

                // Check for errors
                if(n>nMax){
                    solver_comp_status = GSL_CONTINUE;
                }
                if (solver_comp_status != GSL_SUCCESS) {
                    n--;
                    stringstream error_msg;
                    error_msg << "comp_step: Stopped at Step: n=" << n << " - r=" << r << " - f[0]=" << f[0];
                    GSL_ERROR_VAL (error_msg.str().c_str(), solver_comp_status,);
                }
            }else{
                solver_comp_status = GSL_CONTINUE;
            }
            //endregion
        }else{
            GSL_ERROR_VAL ("gsl_odeiv2_solver::comp_step: solve_status>3 required.", GSL_FAILURE,);
        }
    };

    void gsl_odeiv2_solver::comp(){
        if(solver_status>3){

            if(mode){rhfi = vector<double>((unsigned long)(data->data_dim));};

            //region Main gsl_odeiv2_evolve_apply Loop
            while (solver_comp_status==GSL_SUCCESS) {
                comp_step();
                if(mode){
                    rhfi[0]=r; rhfi[1]=h;
                    for(int i=2;i<data->data_dim;i=i+2){
                        rhfi[i]=f[i/2-1];
                        rhfi[i+1]=(odeiv2_evolve->dydt_out)[i/2-1];
                    }
                    (*data).put(n-1, rhfi);
                };
            }
            solver_status=5;
            if(print){
                printf("|--: gsl_odeiv2_solver::comp: odeiv2_evolve converged at step %.d(+1) nodes (%lu failed (%d%%), %lu computed), r_n=%.4E\n",
                       n-1,odeiv2_evolve->failed_steps,int(double(odeiv2_evolve->failed_steps)/double(odeiv2_evolve->count)*1E2),odeiv2_evolve->count,r);
                for(int i=0;i<dim;i++){
                    printf("|--: gsl_odeiv2_solver::comp: f%d_n=%.4E, df%d_n=%.4E, ferr_max=%d %.4E\n",i,f[i],i,df[i],i,ferr_max[i]);
                }
            }
        }else{
            GSL_ERROR_VAL ("gsl_odeiv2_solver::comp: solve_status>3 required.", GSL_FAILURE,);
        }
    };

    void gsl_odeiv2_solver::get_nrh(unsigned int &n_in, double &r_in, double &h_in) const {
        n_in=n;
        r_in=r;
        h_in=h;
    };
}