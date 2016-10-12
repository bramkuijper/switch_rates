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

double bound0(double val)
{
    if (val < 0)
    {
        val = 0.001;
    }

    return(val);
}

double bound01(double val)
{
    val = val < 0 ? 0.0001 : val > 1.0 ? 0.9999 : val;

    return(val);
}


// equations of all the patch frequencies
int psys_equations(void *params, gsl_vector *f)
{
    //VARFUNC
    
    
    // PATCHRECUR
    
    return GSL_SUCCESS;
}

// equations of the relatedness coefficients
int rel_equations(void *params, gsl_vector *f)
{
    //VARFUNC
    
    
    // RELVALS
    
    return GSL_SUCCESS;
}


// equations of the reproductive values
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
    long unsigned int max_iter = atoi(argv[1]);

    // initialize the vectors that contain the variables
    // functions solve for
    //
    // vector for patch frequencies
    gsl_vector *x_p = gsl_vector_alloc(6);
    
    // vector for reproductive values
    gsl_vector *x_v = gsl_vector_alloc(8);
    
    // vector for relatedness
    gsl_vector *x_r = gsl_vector_alloc(6);

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
    long unsigned int iter;
    for (iter = 0; iter < max_iter ; ++iter)
    {
        // patch frequencies
        psys_equations(&paramstruct, x_p);

        paramstruct.f1aa = gsl_vector_get(x_p,0);
        paramstruct.f1am = gsl_vector_get(x_p,1);
        paramstruct.f1mm = gsl_vector_get(x_p,2);
        paramstruct.f2aa = gsl_vector_get(x_p,3);
        paramstruct.f2am = gsl_vector_get(x_p,4);
        paramstruct.f2mm = gsl_vector_get(x_p,5);

        assert(isnan(paramstruct.f1aa) == 0);
        assert(isnan(paramstruct.f1am) == 0);
        assert(isnan(paramstruct.f1mm) == 0);
        assert(isnan(paramstruct.f2aa) == 0);
        assert(isnan(paramstruct.f2am) == 0);
        assert(isnan(paramstruct.f2mm) == 0);

        // reproductive values
        reproductive_values(&paramstruct, x_v);
        
        paramstruct.v1aa = gsl_vector_get(x_v, 0);
        paramstruct.v1am = gsl_vector_get(x_v, 1);
        paramstruct.v1ma = gsl_vector_get(x_v, 2);
        paramstruct.v1mm = gsl_vector_get(x_v, 3);
        paramstruct.v2aa = gsl_vector_get(x_v, 4);
        paramstruct.v2am = gsl_vector_get(x_v, 5);
        paramstruct.v2ma = gsl_vector_get(x_v, 6);
        paramstruct.v2mm = gsl_vector_get(x_v, 7);

        assert(isnan(paramstruct.v1aa) == 0);
        assert(isnan(paramstruct.v1am) == 0);
        assert(isnan(paramstruct.v1ma) == 0);
        assert(isnan(paramstruct.v1mm) == 0);
        assert(isnan(paramstruct.v2aa) == 0);
        assert(isnan(paramstruct.v2am) == 0);
        assert(isnan(paramstruct.v2ma) == 0);
        assert(isnan(paramstruct.v2mm) == 0);

        // relatedness
        rel_equations(&paramstruct, x_r);
        
        paramstruct.raa1 = gsl_vector_get(x_r, 0);
        paramstruct.ram1 = gsl_vector_get(x_r, 1);
        paramstruct.rmm1 = gsl_vector_get(x_r, 2);
        paramstruct.raa2 = gsl_vector_get(x_r, 3);
        paramstruct.ram2 = gsl_vector_get(x_r, 4);
        paramstruct.rmm2 = gsl_vector_get(x_r, 5);
        
        assert(isnan(paramstruct.raa1) == 0);
        assert(isnan(paramstruct.ram1) == 0);
        assert(isnan(paramstruct.rmm1) == 0);
        assert(isnan(paramstruct.raa2) == 0);
        assert(isnan(paramstruct.ram2) == 0);
        assert(isnan(paramstruct.rmm2) == 0);

        // selection gradients
        selgrads(&paramstruct, x_selgrad);

        bool condition_zs1 = fabs(paramstruct.zs1 - gsl_vector_get(x_selgrad, 2)) < 1e-10; 
        bool condition_zs2 = fabs(paramstruct.zs2 - gsl_vector_get(x_selgrad, 3)) < 1e-10; 
        bool condition_p1 = (fabs(paramstruct.p1 - gsl_vector_get(x_selgrad, 0)) < 1e-10 || paramstruct.p1 >= 0.999) || paramstruct.p1 <= 0.001;
        bool condition_p2 = (fabs(paramstruct.p2 - gsl_vector_get(x_selgrad, 1)) < 1e-10 || paramstruct.p2 >= 0.999) || paramstruct.p2 <= 0.001;

        if (condition_zs1 && condition_zs2 && condition_p1 && condition_p2)
        {
            paramstruct.p1 = gsl_vector_get(x_selgrad, 0);
            paramstruct.p2 = gsl_vector_get(x_selgrad, 1);
            paramstruct.zs1 = gsl_vector_get(x_selgrad, 2);
            paramstruct.zs2 = gsl_vector_get(x_selgrad, 3);

            break;
        }
        paramstruct.p1 = gsl_vector_get(x_selgrad, 0);
        paramstruct.p2 = gsl_vector_get(x_selgrad, 1);
        paramstruct.zs1 = gsl_vector_get(x_selgrad, 2);
        paramstruct.zs2 = gsl_vector_get(x_selgrad, 3);
        
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
