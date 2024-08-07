//
// Created by M. J. Steil on 2017.08.15.
//

#include "../gsl_wrapper.hpp"

namespace gsl_wrapper{

    gsl_function_pp_fkt::gsl_function_pp_fkt(std::function<double(double)> func) : gsl_function_struct(), _func(std::move(func)){
        function=&gsl_function_pp_fkt::invoke;
        params=this;
    }

    double gsl_function_pp_fkt::invoke(double x, void *params) {
        return static_cast<gsl_function_pp_fkt*>(params)->_func(x);
    }

    gsl_function *gsl_function_pp_fkt::gsl_fkt() {
        return static_cast<gsl_function*>(this);
    }

}