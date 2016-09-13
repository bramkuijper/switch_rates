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
    // 
double mum;
double mua;
double d;
double s1;
double s2;
double k;
double p1;
double p2;
double zs1;
double zs2;
double f1a;
double f2a;
double f1m;
double f2m;
double v1a;
double v1m;
double v2a;
double v2m;

};

double bound(double val)
{
    val = val < 0 ? 0.0001 : val > 1.0 ? 0.9999 : val;

    return(val);
}


// recursions of all the patch frequencies
int psys_recur(void *params, gsl_vector *f)
{
    //
double mum = ((struct rparams *) params)->mum;
double mua = ((struct rparams *) params)->mua;
double d = ((struct rparams *) params)->d;
double s1 = ((struct rparams *) params)->s1;
double s2 = ((struct rparams *) params)->s2;
double k = ((struct rparams *) params)->k;
double p1 = ((struct rparams *) params)->p1;
double p2 = ((struct rparams *) params)->p2;
double zs1 = ((struct rparams *) params)->zs1;
double zs2 = ((struct rparams *) params)->zs2;
double f1a = ((struct rparams *) params)->f1a;
double f2a = ((struct rparams *) params)->f2a;
double f1m = ((struct rparams *) params)->f1m;
double f2m = ((struct rparams *) params)->f2m;
double v1a = ((struct rparams *) params)->v1a;
double v1m = ((struct rparams *) params)->v1m;
double v2a = ((struct rparams *) params)->v2a;
double v2m = ((struct rparams *) params)->v2m;

    
    
    // 
double f1atplus1, f1mtplus1, f2atplus1, f2mtplus1;
for (int iter=0; iter<1e08; ++iter) {
f1atplus1 = bound(f1a + 0.01*(-(f1a*mua*(1 - (1 - d)*p1 - d*(f1a*p1 + f2m*p1 + f1m*p2 + f2a*p2))*(1 + k*pow(-1 + zs1,2))) + f1m*mum*((1 - d)*p2 + d*(f1a*p1 + f2m*p1 + f1m*p2 + f2a*p2))*(1 + k*pow(-1 + zs1,2)) - f1a*s1*zs1 + f2m*s2*zs2));

f1mtplus1 = bound(f1m + 0.01*(f1a*mua*(1 - (1 - d)*p1 - d*(f1a*p1 + f2m*p1 + f1m*p2 + f2a*p2))*(1 + k*pow(-1 + zs1,2)) - f1m*mum*((1 - d)*p2 + d*(f1a*p1 + f2m*p1 + f1m*p2 + f2a*p2))*(1 + k*pow(-1 + zs1,2)) - f1m*s1*zs1 + f2a*s2*zs2));

f2atplus1 = bound(f2a + 0.01*(f1m*s1*zs1 + f2m*mum*((1 - d)*(1 - p1) + d*(f1a*(1 - p1) + f2m*(1 - p1) + f1m*(1 - p2) + f2a*(1 - p2)))*(1 + k*pow(-1 + zs2,2)) - f2a*mua*(1 - d*(f1a*(1 - p1) + f2m*(1 - p1) + f1m*(1 - p2) + f2a*(1 - p2)) - (1 - d)*(1 - p2))*(1 + k*pow(-1 + zs2,2)) - f2a*s2*zs2));

f2mtplus1 = bound(f2m + 0.01*(f1a*s1*zs1 - f2m*mum*((1 - d)*(1 - p1) + d*(f1a*(1 - p1) + f2m*(1 - p1) + f1m*(1 - p2) + f2a*(1 - p2)))*(1 + k*pow(-1 + zs2,2)) + f2a*mua*(1 - d*(f1a*(1 - p1) + f2m*(1 - p1) + f1m*(1 - p2) + f2a*(1 - p2)) - (1 - d)*(1 - p2))*(1 + k*pow(-1 + zs2,2)) - f2m*s2*zs2));

if (fabs(f1atplus1 - f1a) < 1e-07 && fabs(f2atplus1 - f2a) < 1e-07 && fabs(f1mtplus1 - f1m) < 1e-07 && fabs(f2mtplus1 - f2m) < 1e-07) {
break;
}
f1a = f1atplus1;
f2a = f2atplus1;
f1m = f1mtplus1;
f2m = f2mtplus1;
}
gsl_vector_set(f, 0, f1a);
gsl_vector_set(f, 1, f2a);
gsl_vector_set(f, 2, f1m);
gsl_vector_set(f, 3, f2m);

    
    return GSL_SUCCESS;
}


void reproductive_values(void *params, gsl_vector *f)
{
    // 
double mum = ((struct rparams *) params)->mum;
double mua = ((struct rparams *) params)->mua;
double d = ((struct rparams *) params)->d;
double s1 = ((struct rparams *) params)->s1;
double s2 = ((struct rparams *) params)->s2;
double k = ((struct rparams *) params)->k;
double p1 = ((struct rparams *) params)->p1;
double p2 = ((struct rparams *) params)->p2;
double zs1 = ((struct rparams *) params)->zs1;
double zs2 = ((struct rparams *) params)->zs2;
double f1a = ((struct rparams *) params)->f1a;
double f2a = ((struct rparams *) params)->f2a;
double f1m = ((struct rparams *) params)->f1m;
double f2m = ((struct rparams *) params)->f2m;
double v1a = ((struct rparams *) params)->v1a;
double v1m = ((struct rparams *) params)->v1m;
double v2a = ((struct rparams *) params)->v2a;
double v2m = ((struct rparams *) params)->v2m;



    // 
double v1atplus1, v2atplus1, v1mtplus1, v2mtplus1;
for (int iter=0; iter<1e08; ++iter) {
v1atplus1 = v1a + 0.01*(-(d*mua*v1a*(1 + k*pow(-1 + zs1,2))) + (1 - d)*mua*(1 - p1)*(-v1a + v1m)*(1 + k*pow(-1 + zs1,2)) + s1*(-v1a + v2m)*zs1 + d*(f1a*mua*(p1*v1a + (1 - p1)*v1m)*(1 + k*pow(-1 + zs1,2)) + f1m*mum*(p1*v1a + (1 - p1)*v1m)*(1 + k*pow(-1 + zs1,2)) + f2a*mua*((1 - p1)*v2a + p1*v2m)*(1 + k*pow(-1 + zs2,2)) + f2m*mum*((1 - p1)*v2a + p1*v2m)*(1 + k*pow(-1 + zs2,2))));

v2atplus1 = v2a + 0.01*(d*(f1a*mua*(p2*v1a + (1 - p2)*v1m)*(1 + k*pow(-1 + zs1,2)) + f1m*mum*(p2*v1a + (1 - p2)*v1m)*(1 + k*pow(-1 + zs1,2)) + f2a*mua*((1 - p2)*v2a + p2*v2m)*(1 + k*pow(-1 + zs2,2)) + f2m*mum*((1 - p2)*v2a + p2*v2m)*(1 + k*pow(-1 + zs2,2))) - d*mua*v2a*(1 + k*pow(-1 + zs2,2)) + (1 - d)*mua*p2*(-v2a + v2m)*(1 + k*pow(-1 + zs2,2)) + s2*(v1m - v2a)*zs2);

v1mtplus1 = v1m + 0.01*((1 - d)*mum*p2*(v1a - v1m)*(1 + k*pow(-1 + zs1,2)) - d*mum*v1m*(1 + k*pow(-1 + zs1,2)) + s1*(-v1m + v2a)*zs1 + d*(f1a*mua*(p2*v1a + (1 - p2)*v1m)*(1 + k*pow(-1 + zs1,2)) + f1m*mum*(p2*v1a + (1 - p2)*v1m)*(1 + k*pow(-1 + zs1,2)) + f2a*mua*((1 - p2)*v2a + p2*v2m)*(1 + k*pow(-1 + zs2,2)) + f2m*mum*((1 - p2)*v2a + p2*v2m)*(1 + k*pow(-1 + zs2,2))));

v2mtplus1 = v2m + 0.01*(d*(f1a*mua*(p1*v1a + (1 - p1)*v1m)*(1 + k*pow(-1 + zs1,2)) + f1m*mum*(p1*v1a + (1 - p1)*v1m)*(1 + k*pow(-1 + zs1,2)) + f2a*mua*((1 - p1)*v2a + p1*v2m)*(1 + k*pow(-1 + zs2,2)) + f2m*mum*((1 - p1)*v2a + p1*v2m)*(1 + k*pow(-1 + zs2,2))) + (1 - d)*mum*(1 - p1)*(v2a - v2m)*(1 + k*pow(-1 + zs2,2)) - d*mum*v2m*(1 + k*pow(-1 + zs2,2)) + s2*(v1a - v2m)*zs2);

if (fabs(v1atplus1 - v1a) < 1e-07 && fabs(v2atplus1 - v2a) < 1e-07 && fabs(v1mtplus1-v1m) < 1e-07 && fabs(v2mtplus1-v2m)) {
break;
}
v1a=v1atplus1;
v2a=v2atplus1;
v1m=v1mtplus1;
v2m=v2mtplus1;
}
gsl_vector_set(f, 0, v1a);
gsl_vector_set(f, 1, v2a);
gsl_vector_set(f, 2, v1m);
gsl_vector_set(f, 3, v2m);

}

void selgrads(void *params, gsl_vector *f)
{
    // 
double mum = ((struct rparams *) params)->mum;
double mua = ((struct rparams *) params)->mua;
double d = ((struct rparams *) params)->d;
double s1 = ((struct rparams *) params)->s1;
double s2 = ((struct rparams *) params)->s2;
double k = ((struct rparams *) params)->k;
double p1 = ((struct rparams *) params)->p1;
double p2 = ((struct rparams *) params)->p2;
double zs1 = ((struct rparams *) params)->zs1;
double zs2 = ((struct rparams *) params)->zs2;
double f1a = ((struct rparams *) params)->f1a;
double f2a = ((struct rparams *) params)->f2a;
double f1m = ((struct rparams *) params)->f1m;
double f2m = ((struct rparams *) params)->f2m;
double v1a = ((struct rparams *) params)->v1a;
double v1m = ((struct rparams *) params)->v1m;
double v2a = ((struct rparams *) params)->v2a;
double v2m = ((struct rparams *) params)->v2m;


    // 
double zs1tplus1 = zs1 + 0.01*(f1a*(s1*(-v1a + v2m) - 2*d*k*mua*v1a*(-1 + zs1) - 2*(-1 + d)*k*mua*(-1 + p1)*(v1a - v1m)*(-1 + zs1)) + f1m*(s1*(-v1m + v2a) - 2*(-1 + d)*k*mum*p2*(v1a - v1m)*(-1 + zs1) - 2*d*k*mum*v1m*(-1 + zs1)));

double zs2tplus1 = zs2 + 0.01*(f1a*(s1*(-v1a + v2m) - 2*d*k*mua*v1a*(-1 + zs1) - 2*(-1 + d)*k*mua*(-1 + p1)*(v1a - v1m)*(-1 + zs1)) + f1m*(s1*(-v1m + v2a) - 2*(-1 + d)*k*mum*p2*(v1a - v1m)*(-1 + zs1) - 2*d*k*mum*v1m*(-1 + zs1)));

double p1tplus1 = bound(p1 + 0.01*(f1a*(-((1 - d)*mua*(-v1a + v1m)*(1 + k*pow(-1 + zs1,2))) + d*(f1a*mua*(v1a - v1m)*(1 + k*pow(-1 + zs1,2)) + f1m*mum*(v1a - v1m)*(1 + k*pow(-1 + zs1,2)) + f2a*mua*(-v2a + v2m)*(1 + k*pow(-1 + zs2,2)) + f2m*mum*(-v2a + v2m)*(1 + k*pow(-1 + zs2,2)))) + f2m*(d*(f1a*mua*(v1a - v1m)*(1 + k*pow(-1 + zs1,2)) + f1m*mum*(v1a - v1m)*(1 + k*pow(-1 + zs1,2)) + f2a*mua*(-v2a + v2m)*(1 + k*pow(-1 + zs2,2)) + f2m*mum*(-v2a + v2m)*(1 + k*pow(-1 + zs2,2))) - (1 - d)*mum*(v2a - v2m)*(1 + k*pow(-1 + zs2,2)))));

double p2tplus1 = bound(p2 + 0.01*(f1m*((1 - d)*mum*(v1a - v1m)*(1 + k*pow(-1 + zs1,2)) + d*(f1a*mua*(v1a - v1m)*(1 + k*pow(-1 + zs1,2)) + f1m*mum*(v1a - v1m)*(1 + k*pow(-1 + zs1,2)) + f2a*mua*(-v2a + v2m)*(1 + k*pow(-1 + zs2,2)) + f2m*mum*(-v2a + v2m)*(1 + k*pow(-1 + zs2,2)))) + f2a*(d*(f1a*mua*(v1a - v1m)*(1 + k*pow(-1 + zs1,2)) + f1m*mum*(v1a - v1m)*(1 + k*pow(-1 + zs1,2)) + f2a*mua*(-v2a + v2m)*(1 + k*pow(-1 + zs2,2)) + f2m*mum*(-v2a + v2m)*(1 + k*pow(-1 + zs2,2))) + (1 - d)*mua*(-v2a + v2m)*(1 + k*pow(-1 + zs2,2)))));

gsl_vector_set(f, 0, zs1tplus1);
gsl_vector_set(f, 1, zs2tplus1);
gsl_vector_set(f, 2, p1tplus1);
gsl_vector_set(f, 3, p2tplus1);


}


void write_params(void *params)
{
    // 
double mum = ((struct rparams *) params)->mum;
double mua = ((struct rparams *) params)->mua;
double d = ((struct rparams *) params)->d;
double s1 = ((struct rparams *) params)->s1;
double s2 = ((struct rparams *) params)->s2;
double k = ((struct rparams *) params)->k;
double p1 = ((struct rparams *) params)->p1;
double p2 = ((struct rparams *) params)->p2;
double zs1 = ((struct rparams *) params)->zs1;
double zs2 = ((struct rparams *) params)->zs2;
double f1a = ((struct rparams *) params)->f1a;
double f2a = ((struct rparams *) params)->f2a;
double f1m = ((struct rparams *) params)->f1m;
double f2m = ((struct rparams *) params)->f2m;
double v1a = ((struct rparams *) params)->v1a;
double v1m = ((struct rparams *) params)->v1m;
double v2a = ((struct rparams *) params)->v2a;
double v2m = ((struct rparams *) params)->v2m;



    // 
DataFile << endl << endl  << "mum;" << mum << endl
 << "mua;" << mua << endl
 << "d;" << d << endl
 << "s1;" << s1 << endl
 << "s2;" << s2 << endl
 << "k;" << k << endl
 << endl;
}


void write_data(void *params, int time)
{
    // 
double mum = ((struct rparams *) params)->mum;
double mua = ((struct rparams *) params)->mua;
double d = ((struct rparams *) params)->d;
double s1 = ((struct rparams *) params)->s1;
double s2 = ((struct rparams *) params)->s2;
double k = ((struct rparams *) params)->k;
double p1 = ((struct rparams *) params)->p1;
double p2 = ((struct rparams *) params)->p2;
double zs1 = ((struct rparams *) params)->zs1;
double zs2 = ((struct rparams *) params)->zs2;
double f1a = ((struct rparams *) params)->f1a;
double f2a = ((struct rparams *) params)->f2a;
double f1m = ((struct rparams *) params)->f1m;
double f2m = ((struct rparams *) params)->f2m;
double v1a = ((struct rparams *) params)->v1a;
double v1m = ((struct rparams *) params)->v1m;
double v2a = ((struct rparams *) params)->v2a;
double v2m = ((struct rparams *) params)->v2m;



    if (time < 0)
    {
        // 
DataFile << "time;f1m;v1a;f2a;zs1;zs2;p1;v1m;p2;f1a;v2a;v2m;f2m;" << endl;
    }

    // 
DataFile << time << ";" << f1m << ";" << 
v1a << ";" << 
f2a << ";" << 
zs1 << ";" << 
zs2 << ";" << 
p1 << ";" << 
v1m << ";" << 
p2 << ";" << 
f1a << ";" << 
v2a << ";" << 
v2m << ";" << 
f2m << ";" << 
 endl;
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
    // 
		paramstruct.mum = atof(argv[2]);
		paramstruct.mua = atof(argv[3]);
		paramstruct.d = atof(argv[4]);
		paramstruct.s1 = atof(argv[5]);
		paramstruct.s2 = atof(argv[6]);
		paramstruct.k = atof(argv[7]);
		paramstruct.p1 = atof(argv[8]);
		paramstruct.p2 = atof(argv[9]);
		paramstruct.zs1 = atof(argv[10]);
		paramstruct.zs2 = atof(argv[11]);
		paramstruct.f1a = atof(argv[12]);
		paramstruct.f2a = atof(argv[13]);
		paramstruct.f1m = atof(argv[14]);
		paramstruct.f2m = atof(argv[15]);
		paramstruct.v1a = atof(argv[16]);
		paramstruct.v1m = atof(argv[17]);
		paramstruct.v2a = atof(argv[18]);
		paramstruct.v2m = atof(argv[19]);

   

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

        bool condition_zs1 = (fabs(paramstruct.zs1 - gsl_vector_get(x_selgrad, 0)) < 1e-7 || paramstruct.zs1 >= 0.999) || paramstruct.zs1 <= 0.001;
        bool condition_zs2 = (fabs(paramstruct.zs2 - gsl_vector_get(x_selgrad, 0)) < 1e-7 || paramstruct.zs2 >= 0.999) || paramstruct.zs2 <= 0.001;
        bool condition_p1 = (fabs(paramstruct.p1 - gsl_vector_get(x_selgrad, 0)) < 1e-7 || paramstruct.p1 >= 0.999) || paramstruct.p1 <= 0.001;
        bool condition_p2 = (fabs(paramstruct.p2 - gsl_vector_get(x_selgrad, 0)) < 1e-7 || paramstruct.p2 >= 0.999) || paramstruct.p2 <= 0.001;

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

        if (iter % 1 == 0)
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

