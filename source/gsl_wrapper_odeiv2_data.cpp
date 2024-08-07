//
// Created by M. J. Steil on 2017.01.19.
//

#include "../gsl_wrapper.hpp"

namespace gsl_wrapper {
    using namespace std;

    // gsl_odeiv2 Data container

    gsl_odeiv2_data::gsl_odeiv2_data(int dim, unsigned long nMax) {
        mode = 1;

        data_nMax = nMax;
        data_dim = dim+2;
        data_vec_size = (data_nMax+2);

        data={vector<double>(data_vec_size), // r_i
              vector<double>(data_vec_size)};// h_i

        for(int i=2;i<data_dim;i++){
            data.emplace_back(data_vec_size); // f_i
        }
    }

    gsl_odeiv2_data::gsl_odeiv2_data() {
        mode = 0;

        data_nMax = 0;
        data_dim = 0;
        data_vec_size = 0;

        data={};
    }

    void gsl_odeiv2_data::put(unsigned long  n, vector<double> rhfi) {
        this->n=n;
        if(rhfi.size()==data_dim){
            for(int i=0;i<data_dim;i++){
                data[i][n]=rhfi[i];
            }
        }else{
            GSL_ERROR_VAL ("gsl_odeiv2_data::put: dimension mismatch rhfi.size()==data_dim.", GSL_FAILURE,);
        }
    }

    void gsl_odeiv2_data::insert(int data_i,vector<double> *in, int pos) {
        data[data_i].insert(data[data_i].begin()+pos,in->begin(),in->begin());
    }

    void gsl_odeiv2_data::scale(double s, int l) {
        if(l<0||l>data_dim){
            GSL_ERROR_VAL ("gsl_odeiv2_data::scale: dimension mismatch l<0||l>data_dim.", GSL_FAILURE,);
        }else{
            transform(data[l].begin(), data[l].end(),data[l].begin(), bind1st(multiplies<double>(),s));
        }
    }

    void gsl_odeiv2_data::offset(double a, int l) {
        if(l<0||l>data_dim){
            GSL_ERROR_VAL ("gsl_odeiv2_data::offset: dimension mismatch l<0||l>dim*2+2.", GSL_FAILURE,);
        }else{
            transform(data[l].begin(), data[l].end(), data[l].begin(), bind1st(plus<double >(),a));
        }
    }


    void gsl_odeiv2_data::add(int l, int m) {
        if (l < 0 || l > data_dim || m < 0 || m > data_dim) {
            GSL_ERROR_VAL ("gsl_odeiv2_data::add: Invalid vector identifier: l<0||l>data_dim||m<0||m>data_dim.", GSL_FAILURE,);
        } else {
            transform(data[l].begin(), data[l].end(), data[m].begin(), data[l].begin(),plus<double>());
        }
    }

    void gsl_odeiv2_data::add_fkt(int l, int m, function<double(double)> fkt){
        if (l < 0 || l > data_dim || m < 0 || m > data_dim) {
            GSL_ERROR_VAL ("gsl_odeiv2_data::add_fkt: Invalid vector identifier: l<0||l>data_dim||m<0||m>data_dim.", GSL_FAILURE,);
        } else {
            vector<double> fkt_of_m(data_nMax);
            for(int i=0; i<data_nMax; i++){
                fkt_of_m[i]=fkt(data[m][i]);
            }

            transform(data[l].begin(), data[l].end(), fkt_of_m.begin(), data[l].begin(),plus<double>());
        }
    }

    void gsl_odeiv2_data::copy(int data_i,vector<double> *target) {
        target->assign(data[data_i].begin(),data[data_i].begin()+n+1);
    }

}