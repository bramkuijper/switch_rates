#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cassert>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include "bramauxiliary.h"


using namespace std;

string filename("iter_freqdep");
string filename_new(create_filename(filename));
ofstream DataFile(filename_new.c_str());  // output file 

struct rparams
{
    // VARS
};

double bound(double val)
{
    val = val < 0 ? 0.0001 : val > 1.0 ? 0.9999 : val;

    return(val);
}


// recursions of all the patch frequencies
int psys_recur(void *params, gsl_vector *f)
{
    //VARFUNC
    
    
    // PATCHRECUR
    
    return GSL_SUCCESS;
}


void reproductive_values(void *params, gsl_vector *f)
{
    // VARFUNC


    // REPVALS
}

void selgrads(void *params, gsl_vector *f)
{
    // VARFUNC

    // SELGRADS

}


void write_params(void *params)
{
    // VARFUNC


    // WRITEPARAMS
}


void write_data(void *params, int time)
{
    // VARFUNC


    if (time < 0)
    {
        // HEADERWRT
    }

    // WRITEDATA
}


int main (int argc, char **argv)
{
    int max_iter = atoi(argv[1]);

    // initialize the vectors that contain the variables
    // functions solve for
    //
    // vector for patch frequencies
    gsl_vector *x_p = gsl_vector_alloc(4);
    
    // vector for reproductive values
    gsl_vector *x_v = gsl_vector_alloc(4);

    // vector for selection gradients
    gsl_vector *x_selgrad = gsl_vector_alloc(4);

    // initialize the struct with parameters
    struct rparams paramstruct; 
  
    // initialize command line argument things
    // see generate_cpp.py 
    // ARGINIT
   

    // ranges for cycling: sometimes numerical iterations
    // won't resolve as solutions slightly cycle around
    // the final value. To see whether cyling behaviour
    // is repetitive over time, however, we create these
    // vectors that stores a series of past values of 
    // the values for the switching rates and compares
    // those to values that are currently found
    gsl_vector *zs1_range = gsl_vector_alloc(10);
    gsl_vector *zs2_range = gsl_vector_alloc(10);
    gsl_vector *p1_range = gsl_vector_alloc(10);
    gsl_vector *p2_range = gsl_vector_alloc(10);

    for (int ik = 0; ik < 10; ++ik)
    {
        // initialize the vector
        gsl_vector_set(zs1_range, ik, 0);        
        gsl_vector_set(zs2_range, ik, 0);        
        gsl_vector_set(p1_range, ik, 0);        
        gsl_vector_set(p2_range, ik, 0);        
    }

    // write the initial data set
    write_data(&paramstruct,-1);

    // iterate
    int iter;
    for (iter = 0; iter < max_iter ; ++iter)
    {
        // patch frequencies
        psys_recur(&paramstruct, x_p);

        paramstruct.f1a = gsl_vector_get(x_p,0);
        paramstruct.f2a = gsl_vector_get(x_p,1);
        paramstruct.f1m = gsl_vector_get(x_p,2);
        paramstruct.f2m = gsl_vector_get(x_p,3);

        assert(isnan(paramstruct.f1a) == 0);
        assert(isnan(paramstruct.f2a) == 0);
        assert(isnan(paramstruct.f1m) == 0);
        assert(isnan(paramstruct.f2m) == 0);

        // reproductive values
        reproductive_values(&paramstruct, x_v);
        
        paramstruct.v1a = gsl_vector_get(x_v, 0);
        paramstruct.v2a = gsl_vector_get(x_v, 1);
        paramstruct.v1m = gsl_vector_get(x_v, 2);
        paramstruct.v2m = gsl_vector_get(x_v, 3);
        
        assert(isnan(paramstruct.v1a) == 0);
        assert(isnan(paramstruct.v2a) == 0);
        assert(isnan(paramstruct.v1m) == 0);
        assert(isnan(paramstruct.v2m) == 0);

        // selection gradients
        selgrads(&paramstruct, x_selgrad);

        bool condition_zs1 = fabs(paramstruct.zs1 - gsl_vector_get(x_selgrad, 0)) < 1e-10; 
        bool condition_zs2 = fabs(paramstruct.zs2 - gsl_vector_get(x_selgrad, 0)) < 1e-10; 
        bool condition_p1 = (fabs(paramstruct.p1 - gsl_vector_get(x_selgrad, 0)) < 1e-10 || paramstruct.p1 >= 0.999) || paramstruct.p1 <= 0.001;
        bool condition_p2 = (fabs(paramstruct.p2 - gsl_vector_get(x_selgrad, 0)) < 1e-10 || paramstruct.p2 >= 0.999) || paramstruct.p2 <= 0.001;

        if (condition_zs1 && condition_zs2 && condition_p1 && condition_p2)
        {
            paramstruct.zs1 = gsl_vector_get(x_selgrad, 0);
            paramstruct.zs2 = gsl_vector_get(x_selgrad, 1);
            paramstruct.p1 = gsl_vector_get(x_selgrad, 2);
            paramstruct.p2 = gsl_vector_get(x_selgrad, 3);

            break;
        }
        paramstruct.zs1 = gsl_vector_get(x_selgrad, 0);
        paramstruct.zs2 = gsl_vector_get(x_selgrad, 1);
        paramstruct.p1 = gsl_vector_get(x_selgrad, 2);
        paramstruct.p2 = gsl_vector_get(x_selgrad, 3);
        
        assert(isnan(paramstruct.zs1) == 0);
        assert(isnan(paramstruct.zs2) == 0);
        assert(isnan(paramstruct.p1) == 0);
        assert(isnan(paramstruct.p2) == 0);
        
        if (iter > 50000)
        {
            bool found_in_range1 = false;
            bool found_in_range2 = false;
            bool found_in_range3 = false;
            bool found_in_range4 = false;

            bool done = false;

            for (int ik = 0; ik < 10; ++ik)
            {
                if (!found_in_range1 && fabs(gsl_vector_get(zs1_range, ik)-paramstruct.zs1) < 1e-10)
                {
                    found_in_range1 = true;
                }

                if (!found_in_range2 && fabs(gsl_vector_get(zs2_range, ik)-paramstruct.zs2) < 1e-10)
                {
                    found_in_range2 = true;
                }
                
                if (!found_in_range3 && fabs(gsl_vector_get(p1_range, ik)-paramstruct.p1) < 1e-10)
                {
                    found_in_range3 = true;
                }
                if (!found_in_range4 && fabs(gsl_vector_get(p2_range, ik)-paramstruct.p2) < 1e-10)
                {
                    found_in_range4 = true;
                }
            }

            if (done)
            {
                break;
            }
        }


        for (int ik = 9; ik > 0; --ik)
        {
            gsl_vector_set(zs1_range, ik, gsl_vector_get(zs1_range, ik - 1));
            gsl_vector_set(zs2_range, ik, gsl_vector_get(zs2_range, ik - 1));
            gsl_vector_set(p1_range, ik, gsl_vector_get(p1_range, ik - 1));
            gsl_vector_set(p2_range, ik, gsl_vector_get(p2_range, ik - 1));
        }

        gsl_vector_set(zs1_range, 0, paramstruct.zs1);
        gsl_vector_set(zs2_range, 0, paramstruct.zs2);
        gsl_vector_set(p1_range, 0, paramstruct.p1);
        gsl_vector_set(p2_range, 0, paramstruct.p2);

        if (iter % 100 == 0)
        {
            write_data(&paramstruct,iter);
        }
    }

    write_data(&paramstruct,iter);
    write_params(&paramstruct);

    gsl_vector_free(x_p);
    gsl_vector_free(x_v);
    gsl_vector_free(x_selgrad);
}
