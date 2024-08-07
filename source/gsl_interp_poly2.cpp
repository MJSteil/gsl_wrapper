//
// Quadratic spline
// Created by M. J. Steil on 2017.01.19.
//

#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>

namespace gsl_wrapper {
    static int
    poly2_init(void *,
                      const double [],
                      const double [],
                      size_t ) {
        return GSL_SUCCESS;
    }

    static
    int
    poly2_eval(const void *,
                      const double x_array[], const double y_array[], size_t size,
                      double x,
                      gsl_interp_accel *a,
                      double *y) {
        size_t index;

        if (a != nullptr) {
            index = gsl_interp_accel_find(a, x_array, size, x);
        } else {
            index = gsl_interp_bsearch(x_array, x, 0, size - 1);
        }

        /* evaluate */
        double y0, y1, y2, x0, x1, x2;

        if (index < 1) {
            x0 = x_array[0];
            x1 = x_array[1];
            x2 = x_array[2];

            y0 = y_array[0];
            y1 = y_array[1];
            y2 = y_array[2];
        } else if (index > size - 2) {
            x0 = x_array[size - 3];
            x1 = x_array[size - 2];
            x2 = x_array[size - 1];

            y0 = y_array[size - 3];
            y1 = y_array[size - 2];
            y2 = y_array[size - 1];
        } else {
            x0 = x_array[index - 1];
            x1 = x_array[index];
            x2 = x_array[index + 1];

            y0 = y_array[index - 1];
            y1 = y_array[index];
            y2 = y_array[index + 1];
        }

        if (x2 - x1 > 0) {
            *y = -((y2 * (x0 - x) * (x1 - x)) / ((x0 - x2) * (-x1 + x2)))
                 - (y1 * (x0 - x) * (x2 - x)) / ((x0 - x1) * (x1 - x2))
                 + (y0 * (x1 - x) * (x2 - x)) / ((x0 - x1) * (x0 - x2));
            return GSL_SUCCESS;
        } else {
            *y = 0.0;
            return GSL_EINVAL;
        }

    }


    static
    int
    poly2_eval_deriv(const void *,
                            const double x_array[], const double y_array[], size_t size,
                            double x,
                            gsl_interp_accel *a,
                            double *dydx) {
        size_t index;

        if (a != nullptr) {
            index = gsl_interp_accel_find(a, x_array, size, x);
        } else {
            index = gsl_interp_bsearch(x_array, x, 0, size - 1);
        }

        /* evaluate */
        double y0, y1, y2, x0, x1, x2;

        if (index < 1) {
            x0 = x_array[0];
            x1 = x_array[1];
            x2 = x_array[2];

            y0 = y_array[0];
            y1 = y_array[1];
            y2 = y_array[2];
        } else if (index > size - 2) {
            x0 = x_array[size - 3];
            x1 = x_array[size - 2];
            x2 = x_array[size - 1];

            y0 = y_array[size - 3];
            y1 = y_array[size - 2];
            y2 = y_array[size - 1];
        } else {
            x0 = x_array[index - 1];
            x1 = x_array[index];
            x2 = x_array[index + 1];

            y0 = y_array[index - 1];
            y1 = y_array[index];
            y2 = y_array[index + 1];
        }

        if (x2 - x1 > 0) {
            *dydx = (y2 * (x0 + x1 - 2 * x)) / ((x0 - x2) * (-x1 + x2))
                    - (y1 * (x0 + x2 - 2 * x)) / ((x0 - x1) * (x1 - x2))
                    + (y0 * (x1 + x2 - 2 * x)) / ((x0 - x1) * (x0 - x2));
            return GSL_SUCCESS;
        } else {
            *dydx = 0.0;
            return GSL_EINVAL;
        }
    }
    static
    int
    poly2_eval_deriv2(const void *,
                             const double x_array[], const double y_array[], size_t size,
                             double x,
                             gsl_interp_accel *a,
                             double *y_pp) {
        size_t index;

        if (a != nullptr) {
            index = gsl_interp_accel_find(a, x_array, size, x);
        } else {
            index = gsl_interp_bsearch(x_array, x, 0, size - 1);
        }

        /* evaluate */
        double y0, y1, y2, x0, x1, x2;

        if (index < 1) {
            x0 = x_array[0];
            x1 = x_array[1];
            x2 = x_array[2];

            y0 = y_array[0];
            y1 = y_array[1];
            y2 = y_array[2];
        } else if (index > size - 2) {
            x0 = x_array[size - 3];
            x1 = x_array[size - 2];
            x2 = x_array[size - 1];

            y0 = y_array[size - 3];
            y1 = y_array[size - 2];
            y2 = y_array[size - 1];
        } else {
            x0 = x_array[index - 1];
            x1 = x_array[index];
            x2 = x_array[index + 1];

            y0 = y_array[index - 1];
            y1 = y_array[index];
            y2 = y_array[index + 1];
        }

        if (x2 - x1 > 0) {
            *y_pp = (2 * y0) / ((x0 - x1) * (x0 - x2))
                    - (2 * y1) / ((x0 - x1) * (x1 - x2))
                    - (2 * y2) / ((x0 - x2) * (-x1 + x2));
            return GSL_SUCCESS;
        } else {
            *y_pp = 0.0;
            return GSL_EINVAL;
        }
    }

    static const gsl_interp_type poly2_type =
            {
                    "poly2",
                    2,
                    nullptr, /* alloc, not applicable */
                    &poly2_init,
                    &poly2_eval,
                    &poly2_eval_deriv,
                    &poly2_eval_deriv2,
                    nullptr, /* integ, not applicable */
                    nullptr, /* free, not applicable */
            };

    const gsl_interp_type *gsl_interp_poly2 = &poly2_type;
}