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
double q;
double eul;
double p1;
double p2;
double zs1;
double zs2;
double f1aa;
double f1am;
double f1mm;
double f2aa;
double f2am;
double f2mm;
double v1aa;
double v1am;
double v1ma;
double v1mm;
double v2aa;
double v2am;
double v2ma;
double v2mm;
double raa1;
double ram1;
double rmm1;
double raa2;
double ram2;
double rmm2;

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
    //
double mum = ((struct rparams *) params)->mum;
double mua = ((struct rparams *) params)->mua;
double d = ((struct rparams *) params)->d;
double s1 = ((struct rparams *) params)->s1;
double s2 = ((struct rparams *) params)->s2;
double k = ((struct rparams *) params)->k;
double q = ((struct rparams *) params)->q;
double eul = ((struct rparams *) params)->eul;
double p1 = ((struct rparams *) params)->p1;
double p2 = ((struct rparams *) params)->p2;
double zs1 = ((struct rparams *) params)->zs1;
double zs2 = ((struct rparams *) params)->zs2;
double f1aa = ((struct rparams *) params)->f1aa;
double f1am = ((struct rparams *) params)->f1am;
double f1mm = ((struct rparams *) params)->f1mm;
double f2aa = ((struct rparams *) params)->f2aa;
double f2am = ((struct rparams *) params)->f2am;
double f2mm = ((struct rparams *) params)->f2mm;
double v1aa = ((struct rparams *) params)->v1aa;
double v1am = ((struct rparams *) params)->v1am;
double v1ma = ((struct rparams *) params)->v1ma;
double v1mm = ((struct rparams *) params)->v1mm;
double v2aa = ((struct rparams *) params)->v2aa;
double v2am = ((struct rparams *) params)->v2am;
double v2ma = ((struct rparams *) params)->v2ma;
double v2mm = ((struct rparams *) params)->v2mm;
double raa1 = ((struct rparams *) params)->raa1;
double ram1 = ((struct rparams *) params)->ram1;
double rmm1 = ((struct rparams *) params)->rmm1;
double raa2 = ((struct rparams *) params)->raa2;
double ram2 = ((struct rparams *) params)->ram2;
double rmm2 = ((struct rparams *) params)->rmm2;

    
    
    // 
double ftplus11aa;

double ftplus11am;

double ftplus11mm;

double ftplus12aa;

double ftplus12am;

double ftplus12mm;

for(int iter = 0; iter < 1e08; ++iter) {

ftplus11aa = bound01(f1aa + eul*(-(pow(exp(1),q*zs1)*f1aa*s1) + pow(exp(1),q*zs2)*f2mm*s2 - 2*f1aa*mua*(1 - (1 - d)*p1 - (d*(2*f1aa*p1 + 2*f2mm*p1 + 2*f1mm*p2 + 2*f2aa*p2 + f1am*(p1 + p2) + f2am*(p1 + p2)))/2.)*(1 + k*pow(zs1,2)) + f1am*mum*(((1 - d)*(p1 + p2))/2. + (d*(2*f1aa*p1 + 2*f2mm*p1 + 2*f1mm*p2 + 2*f2aa*p2 + f1am*(p1 + p2) + f2am*(p1 + p2)))/2.)*(1 + k*pow(zs1,2))));

ftplus11am = bound01(f1am + eul*(-(pow(exp(1),q*zs1)*f1am*s1) + pow(exp(1),q*zs2)*f2am*s2 + 2*f1aa*mua*(1 - (1 - d)*p1 - (d*(2*f1aa*p1 + 2*f2mm*p1 + 2*f1mm*p2 + 2*f2aa*p2 + f1am*(p1 + p2) + f2am*(p1 + p2)))/2.)*(1 + k*pow(zs1,2)) - f1am*mua*(1 - ((1 - d)*(p1 + p2))/2. - (d*(2*f1aa*p1 + 2*f2mm*p1 + 2*f1mm*p2 + 2*f2aa*p2 + f1am*(p1 + p2) + f2am*(p1 + p2)))/2.)*(1 + k*pow(zs1,2)) + 2*f1mm*mum*((1 - d)*p2 + (d*(2*f1aa*p1 + 2*f2mm*p1 + 2*f1mm*p2 + 2*f2aa*p2 + f1am*(p1 + p2) + f2am*(p1 + p2)))/2.)*(1 + k*pow(zs1,2)) - f1am*mum*(((1 - d)*(p1 + p2))/2. + (d*(2*f1aa*p1 + 2*f2mm*p1 + 2*f1mm*p2 + 2*f2aa*p2 + f1am*(p1 + p2) + f2am*(p1 + p2)))/2.)*(1 + k*pow(zs1,2))));

ftplus11mm = bound01(f1mm + eul*(-(pow(exp(1),q*zs1)*f1mm*s1) + pow(exp(1),q*zs2)*f2aa*s2 + f1am*mua*(1 - ((1 - d)*(p1 + p2))/2. - (d*(2*f1aa*p1 + 2*f2mm*p1 + 2*f1mm*p2 + 2*f2aa*p2 + f1am*(p1 + p2) + f2am*(p1 + p2)))/2.)*(1 + k*pow(zs1,2)) - 2*f1mm*mum*((1 - d)*p2 + (d*(2*f1aa*p1 + 2*f2mm*p1 + 2*f1mm*p2 + 2*f2aa*p2 + f1am*(p1 + p2) + f2am*(p1 + p2)))/2.)*(1 + k*pow(zs1,2))));

ftplus12aa = bound01(f2aa + eul*(pow(exp(1),q*zs1)*f1mm*s1 - pow(exp(1),q*zs2)*f2aa*s2 - 2*f2aa*mua*(1 - (d*(f1aa*(2 - 2*p1) + f2mm*(2 - 2*p1) + f1mm*(2 - 2*p2) + f2aa*(2 - 2*p2) + f1am*(2 - p1 - p2) + f2am*(2 - p1 - p2)))/2. - ((1 - d)*(2 - 2*p2))/2.)*(1 + k*pow(zs2,2)) + f2am*mum*((d*(f1aa*(2 - 2*p1) + f2mm*(2 - 2*p1) + f1mm*(2 - 2*p2) + f2aa*(2 - 2*p2) + f1am*(2 - p1 - p2) + f2am*(2 - p1 - p2)))/2. + ((1 - d)*(2 - p1 - p2))/2.)*(1 + k*pow(zs2,2))));

ftplus12am = bound01(f2am + eul*(pow(exp(1),q*zs1)*f1am*s1 - pow(exp(1),q*zs2)*f2am*s2 + 2*f2mm*mum*(((1 - d)*(2 - 2*p1))/2. + (d*(f1aa*(2 - 2*p1) + f2mm*(2 - 2*p1) + f1mm*(2 - 2*p2) + f2aa*(2 - 2*p2) + f1am*(2 - p1 - p2) + f2am*(2 - p1 - p2)))/2.)*(1 + k*pow(zs2,2)) + 2*f2aa*mua*(1 - (d*(f1aa*(2 - 2*p1) + f2mm*(2 - 2*p1) + f1mm*(2 - 2*p2) + f2aa*(2 - 2*p2) + f1am*(2 - p1 - p2) + f2am*(2 - p1 - p2)))/2. - ((1 - d)*(2 - 2*p2))/2.)*(1 + k*pow(zs2,2)) - f2am*mua*(1 - (d*(f1aa*(2 - 2*p1) + f2mm*(2 - 2*p1) + f1mm*(2 - 2*p2) + f2aa*(2 - 2*p2) + f1am*(2 - p1 - p2) + f2am*(2 - p1 - p2)))/2. - ((1 - d)*(2 - p1 - p2))/2.)*(1 + k*pow(zs2,2)) - f2am*mum*((d*(f1aa*(2 - 2*p1) + f2mm*(2 - 2*p1) + f1mm*(2 - 2*p2) + f2aa*(2 - 2*p2) + f1am*(2 - p1 - p2) + f2am*(2 - p1 - p2)))/2. + ((1 - d)*(2 - p1 - p2))/2.)*(1 + k*pow(zs2,2))));

ftplus12mm = bound01(1 - ftplus11aa - ftplus11am - ftplus11mm - ftplus12aa - ftplus12am);

if (

fabs(ftplus11aa - f1aa)<=1e-07&&

fabs(ftplus11am - f1am)<=1e-07&&

fabs(ftplus11mm - f1mm)<=1e-07&&

fabs(ftplus12aa - f2aa)<=1e-07&&

fabs(ftplus12am - f2am)<=1e-07&&

fabs(ftplus12mm - f2mm)<=1e-07)

{break;}

f1aa = ftplus11aa;

f1am = ftplus11am;

f1mm = ftplus11mm;

f2aa = ftplus12aa;

f2am = ftplus12am;

f2mm = ftplus12mm;

}

gsl_vector_set(f,0,f1aa);
gsl_vector_set(f,1,f1am);
gsl_vector_set(f,2,f1mm);
gsl_vector_set(f,3,f2aa);
gsl_vector_set(f,4,f2am);
gsl_vector_set(f,5,f2mm);

    
    return GSL_SUCCESS;
}

// equations of the relatedness coefficients
int rel_equations(void *params, gsl_vector *f)
{
    //
double mum = ((struct rparams *) params)->mum;
double mua = ((struct rparams *) params)->mua;
double d = ((struct rparams *) params)->d;
double s1 = ((struct rparams *) params)->s1;
double s2 = ((struct rparams *) params)->s2;
double k = ((struct rparams *) params)->k;
double q = ((struct rparams *) params)->q;
double eul = ((struct rparams *) params)->eul;
double p1 = ((struct rparams *) params)->p1;
double p2 = ((struct rparams *) params)->p2;
double zs1 = ((struct rparams *) params)->zs1;
double zs2 = ((struct rparams *) params)->zs2;
double f1aa = ((struct rparams *) params)->f1aa;
double f1am = ((struct rparams *) params)->f1am;
double f1mm = ((struct rparams *) params)->f1mm;
double f2aa = ((struct rparams *) params)->f2aa;
double f2am = ((struct rparams *) params)->f2am;
double f2mm = ((struct rparams *) params)->f2mm;
double v1aa = ((struct rparams *) params)->v1aa;
double v1am = ((struct rparams *) params)->v1am;
double v1ma = ((struct rparams *) params)->v1ma;
double v1mm = ((struct rparams *) params)->v1mm;
double v2aa = ((struct rparams *) params)->v2aa;
double v2am = ((struct rparams *) params)->v2am;
double v2ma = ((struct rparams *) params)->v2ma;
double v2mm = ((struct rparams *) params)->v2mm;
double raa1 = ((struct rparams *) params)->raa1;
double ram1 = ((struct rparams *) params)->ram1;
double rmm1 = ((struct rparams *) params)->rmm1;
double raa2 = ((struct rparams *) params)->raa2;
double ram2 = ((struct rparams *) params)->ram2;
double rmm2 = ((struct rparams *) params)->rmm2;

    
    
    // 
double raatplus11;

double ramtplus11;

double rmmtplus11;

double raatplus12;

double ramtplus12;

double rmmtplus12;

for(int iter = 0; iter < 1e08; ++iter) {

raatplus11 = bound01((pow(exp(1),q*zs2)*f2mm*rmm2*s2 + (1 - d)*f1aa*mua*p1*(1 + raa1)*(1 + k*pow(zs1,2)) + (1 - d)*f1am*mum*(p1/2. + (p2*ram1)/2.)*(1 + k*pow(zs1,2)))/(pow(exp(1),q*zs2)*f2mm*s2 + 2*f1aa*mua*((1 - d)*p1 + (d*(2*f1aa*p1 + 2*f2mm*p1 + 2*f1mm*p2 + 2*f2aa*p2 + f1am*(p1 + p2) + f2am*(p1 + p2)))/2.)*(1 + k*pow(zs1,2)) + f1am*mum*(((1 - d)*(p1 + p2))/2. + (d*(2*f1aa*p1 + 2*f2mm*p1 + 2*f1mm*p2 + 2*f2aa*p2 + f1am*(p1 + p2) + f2am*(p1 + p2)))/2.)*(1 + k*pow(zs1,2))));

ramtplus11 = bound01((pow(exp(1),q*zs2)*f2am*ram2*s2 + (1 - d)*f1aa*mua*(1 - p1)*(1 + raa1)*(1 + k*pow(zs1,2)) + (1 - d)*f1am*mua*(p2/2. + (p1*ram1)/2.)*(1 + k*pow(zs1,2)) + (1 - d)*f1am*mum*((1 - p1)/2. + ((1 - p2)*ram1)/2.)*(1 + k*pow(zs1,2)) + (1 - d)*f1mm*mum*p2*(1 + rmm1)*(1 + k*pow(zs1,2)))/(pow(exp(1),q*zs2)*f2am*s2 + 2*f1aa*mua*(1 - (1 - d)*p1 - (d*(2*f1aa*p1 + 2*f2mm*p1 + 2*f1mm*p2 + 2*f2aa*p2 + f1am*(p1 + p2) + f2am*(p1 + p2)))/2.)*(1 + k*pow(zs1,2)) + f1am*mum*(1 - ((1 - d)*(p1 + p2))/2. - (d*(2*f1aa*p1 + 2*f2mm*p1 + 2*f1mm*p2 + 2*f2aa*p2 + f1am*(p1 + p2) + f2am*(p1 + p2)))/2.)*(1 + k*pow(zs1,2)) + 2*f1mm*mum*((1 - d)*p2 + (d*(2*f1aa*p1 + 2*f2mm*p1 + 2*f1mm*p2 + 2*f2aa*p2 + f1am*(p1 + p2) + f2am*(p1 + p2)))/2.)*(1 + k*pow(zs1,2)) + f1am*mua*(((1 - d)*(p1 + p2))/2. + (d*(2*f1aa*p1 + 2*f2mm*p1 + 2*f1mm*p2 + 2*f2aa*p2 + f1am*(p1 + p2) + f2am*(p1 + p2)))/2.)*(1 + k*pow(zs1,2))));

rmmtplus11 = bound01((pow(exp(1),q*zs2)*f2aa*raa2*s2 + ((1 - d)*f1am*mua*(1 - p2 + (1 - p1)*ram1)*(1 + k*pow(zs1,2)))/2. + (1 - d)*f1mm*mum*(1 - p2)*(1 + rmm1)*(1 + k*pow(zs1,2)))/(pow(exp(1),q*zs2)*f2aa*s2 + 2*f1mm*mum*(1 - (1 - d)*p2 - (d*(2*f1aa*p1 + 2*f2mm*p1 + 2*f1mm*p2 + 2*f2aa*p2 + f1am*(p1 + p2) + f2am*(p1 + p2)))/2.)*(1 + k*pow(zs1,2)) + f1am*mua*(1 - ((1 - d)*(p1 + p2))/2. - (d*(2*f1aa*p1 + 2*f2mm*p1 + 2*f1mm*p2 + 2*f2aa*p2 + f1am*(p1 + p2) + f2am*(p1 + p2)))/2.)*(1 + k*pow(zs1,2))));

raatplus12 = bound01((pow(exp(1),q*zs1)*f1mm*rmm1*s1 + (1 - d)*f2aa*mua*(1 - p2)*(1 + raa2)*(1 + k*pow(zs2,2)) + (1 - d)*f2am*mum*((1 - p2)/2. + ((1 - p1)*ram2)/2.)*(1 + k*pow(zs2,2)))/(pow(exp(1),q*zs1)*f1mm*s1 + 2*f2aa*mua*((d*(f1aa*(2 - 2*p1) + f2mm*(2 - 2*p1) + f1mm*(2 - 2*p2) + f2aa*(2 - 2*p2) + f1am*(2 - p1 - p2) + f2am*(2 - p1 - p2)))/2. + ((1 - d)*(2 - 2*p2))/2.)*(1 + k*pow(zs2,2)) + f2am*mum*((d*(f1aa*(2 - 2*p1) + f2mm*(2 - 2*p1) + f1mm*(2 - 2*p2) + f2aa*(2 - 2*p2) + f1am*(2 - p1 - p2) + f2am*(2 - p1 - p2)))/2. + ((1 - d)*(2 - p1 - p2))/2.)*(1 + k*pow(zs2,2))));

ramtplus12 = bound01((pow(exp(1),q*zs1)*f1am*ram1*s1 + (1 - d)*f2aa*mua*p2*(1 + raa2)*(1 + k*pow(zs2,2)) + (1 - d)*f2am*mum*(p2/2. + (p1*ram2)/2.)*(1 + k*pow(zs2,2)) + (1 - d)*f2am*mua*((1 - p1)/2. + ((1 - p2)*ram2)/2.)*(1 + k*pow(zs2,2)) + (1 - d)*f2mm*mum*(1 - p1)*(1 + rmm2)*(1 + k*pow(zs2,2)))/(pow(exp(1),q*zs1)*f1am*s1 + 2*f2mm*mum*(((1 - d)*(2 - 2*p1))/2. + (d*(f1aa*(2 - 2*p1) + f2mm*(2 - 2*p1) + f1mm*(2 - 2*p2) + f2aa*(2 - 2*p2) + f1am*(2 - p1 - p2) + f2am*(2 - p1 - p2)))/2.)*(1 + k*pow(zs2,2)) + 2*f2aa*mua*(1 - (d*(f1aa*(2 - 2*p1) + f2mm*(2 - 2*p1) + f1mm*(2 - 2*p2) + f2aa*(2 - 2*p2) + f1am*(2 - p1 - p2) + f2am*(2 - p1 - p2)))/2. - ((1 - d)*(2 - 2*p2))/2.)*(1 + k*pow(zs2,2)) + f2am*mum*(1 - (d*(f1aa*(2 - 2*p1) + f2mm*(2 - 2*p1) + f1mm*(2 - 2*p2) + f2aa*(2 - 2*p2) + f1am*(2 - p1 - p2) + f2am*(2 - p1 - p2)))/2. - ((1 - d)*(2 - p1 - p2))/2.)*(1 + k*pow(zs2,2)) + f2am*mua*((d*(f1aa*(2 - 2*p1) + f2mm*(2 - 2*p1) + f1mm*(2 - 2*p2) + f2aa*(2 - 2*p2) + f1am*(2 - p1 - p2) + f2am*(2 - p1 - p2)))/2. + ((1 - d)*(2 - p1 - p2))/2.)*(1 + k*pow(zs2,2))));

rmmtplus12 = bound01((pow(exp(1),q*zs1)*f1aa*raa1*s1 + ((1 - d)*f2am*mua*(p1 + p2*ram2)*(1 + k*pow(zs2,2)))/2. + (1 - d)*f2mm*mum*p1*(1 + rmm2)*(1 + k*pow(zs2,2)))/(pow(exp(1),q*zs1)*f1aa*s1 + 2*f2mm*mum*(1 - ((1 - d)*(2 - 2*p1))/2. - (d*(f1aa*(2 - 2*p1) + f2mm*(2 - 2*p1) + f1mm*(2 - 2*p2) + f2aa*(2 - 2*p2) + f1am*(2 - p1 - p2) + f2am*(2 - p1 - p2)))/2.)*(1 + k*pow(zs2,2)) + f2am*mua*(1 - (d*(f1aa*(2 - 2*p1) + f2mm*(2 - 2*p1) + f1mm*(2 - 2*p2) + f2aa*(2 - 2*p2) + f1am*(2 - p1 - p2) + f2am*(2 - p1 - p2)))/2. - ((1 - d)*(2 - p1 - p2))/2.)*(1 + k*pow(zs2,2))));

if (

fabs(raatplus11 - raa1)<=1e-07&&

fabs(ramtplus11 - ram1)<=1e-07&&

fabs(rmmtplus11 - rmm1)<=1e-07&&

fabs(raatplus12 - raa2)<=1e-07&&

fabs(ramtplus12 - ram2)<=1e-07&&

fabs(rmmtplus12 - rmm2)<=1e-07)

{break;}

raa1 = raatplus11;

ram1 = ramtplus11;

rmm1 = rmmtplus11;

raa2 = raatplus12;

ram2 = ramtplus12;

rmm2 = rmmtplus12;

}

gsl_vector_set(f,0,raa1);
gsl_vector_set(f,1,ram1);
gsl_vector_set(f,2,rmm1);
gsl_vector_set(f,3,raa2);
gsl_vector_set(f,4,ram2);
gsl_vector_set(f,5,rmm2);

    return GSL_SUCCESS;
}


// equations of the reproductive values
void reproductive_values(void *params, gsl_vector *f)
{
    // 
double mum = ((struct rparams *) params)->mum;
double mua = ((struct rparams *) params)->mua;
double d = ((struct rparams *) params)->d;
double s1 = ((struct rparams *) params)->s1;
double s2 = ((struct rparams *) params)->s2;
double k = ((struct rparams *) params)->k;
double q = ((struct rparams *) params)->q;
double eul = ((struct rparams *) params)->eul;
double p1 = ((struct rparams *) params)->p1;
double p2 = ((struct rparams *) params)->p2;
double zs1 = ((struct rparams *) params)->zs1;
double zs2 = ((struct rparams *) params)->zs2;
double f1aa = ((struct rparams *) params)->f1aa;
double f1am = ((struct rparams *) params)->f1am;
double f1mm = ((struct rparams *) params)->f1mm;
double f2aa = ((struct rparams *) params)->f2aa;
double f2am = ((struct rparams *) params)->f2am;
double f2mm = ((struct rparams *) params)->f2mm;
double v1aa = ((struct rparams *) params)->v1aa;
double v1am = ((struct rparams *) params)->v1am;
double v1ma = ((struct rparams *) params)->v1ma;
double v1mm = ((struct rparams *) params)->v1mm;
double v2aa = ((struct rparams *) params)->v2aa;
double v2am = ((struct rparams *) params)->v2am;
double v2ma = ((struct rparams *) params)->v2ma;
double v2mm = ((struct rparams *) params)->v2mm;
double raa1 = ((struct rparams *) params)->raa1;
double ram1 = ((struct rparams *) params)->ram1;
double rmm1 = ((struct rparams *) params)->rmm1;
double raa2 = ((struct rparams *) params)->raa2;
double ram2 = ((struct rparams *) params)->ram2;
double rmm2 = ((struct rparams *) params)->rmm2;



    // 
double vtplus11aa;

double vtplus11am;

double vtplus11ma;

double vtplus11mm;

double vtplus12aa;

double vtplus12am;

double vtplus12ma;

double vtplus12mm;

for(int iter = 0; iter < 1e08; ++iter) {

vtplus11aa = bound0(v1aa + eul*(pow(exp(1),q*zs1)*s1*(-v1aa + v2mm) - (1 + (-1 + d)/2.)*mua*v1aa*(1 + k*pow(zs1,2)) + (mua*((1 - d)*(1 - p1) + d*(f1aa*(2 - 2*p1) + f2mm*(2 - 2*p1) + f1mm*(2 - 2*p2) + f2aa*(2 - 2*p2) + f1am*(2 - p1 - p2) + f2am*(2 - p1 - p2)))*(-v1aa + v1am)*(1 + k*pow(zs1,2)))/2. + ((1 - d)*mua*(1 - p1)*(-v1aa + v1ma)*(1 + k*pow(zs1,2)))/2. + ((1 - d)*mua*(p1*v1aa + (1 - p1)*(-v1aa + v1am + v1ma))*(1 + k*pow(zs1,2)))/2. + (d*(2*f1aa*mua*(p1*v1aa + (1 - p1)*v1ma)*(1 + k*pow(zs1,2)) + 2*f1mm*mum*(p1*v1am + (1 - p1)*v1mm)*(1 + k*pow(zs1,2)) + f1am*(mum*(p1*v1aa + (1 - p1)*v1ma)*(1 + k*pow(zs1,2)) + mua*(p1*v1am + (1 - p1)*v1mm)*(1 + k*pow(zs1,2))) + 2*f2aa*mua*((1 - p1)*v2aa + p1*v2ma)*(1 + k*pow(zs2,2)) + 2*f2mm*mum*((1 - p1)*v2am + p1*v2mm)*(1 + k*pow(zs2,2)) + f2am*(mum*((1 - p1)*v2aa + p1*v2ma)*(1 + k*pow(zs2,2)) + mua*((1 - p1)*v2am + p1*v2mm)*(1 + k*pow(zs2,2)))))/2.));

vtplus11am = bound0(v1am + eul*(pow(exp(1),q*zs1)*s1*(-v1am + v2ma) + mum*(((1 - d)*p2)/2. + (d*(2*f1aa*p1 + 2*f2mm*p1 + 2*f1mm*p2 + 2*f2aa*p2 + f1am*(p1 + p2) + f2am*(p1 + p2)))/2.)*(v1aa - v1am)*(1 + k*pow(zs1,2)) - (1 + (-1 + d)/2.)*mua*v1am*(1 + k*pow(zs1,2)) + ((1 - d)*mum*(p1*(2*v1aa - v1am) + (1 - p1)*v1ma)*(1 + k*pow(zs1,2)))/2. + ((1 - d)*mua*(1 - p1)*(-v1am + v1mm)*(1 + k*pow(zs1,2)))/2. + (d*(2*f1aa*mua*(p1*v1aa + (1 - p1)*v1ma)*(1 + k*pow(zs1,2)) + 2*f1mm*mum*(p1*v1am + (1 - p1)*v1mm)*(1 + k*pow(zs1,2)) + f1am*(mum*(p1*v1aa + (1 - p1)*v1ma)*(1 + k*pow(zs1,2)) + mua*(p1*v1am + (1 - p1)*v1mm)*(1 + k*pow(zs1,2))) + 2*f2aa*mua*((1 - p1)*v2aa + p1*v2ma)*(1 + k*pow(zs2,2)) + 2*f2mm*mum*((1 - p1)*v2am + p1*v2mm)*(1 + k*pow(zs2,2)) + f2am*(mum*((1 - p1)*v2aa + p1*v2ma)*(1 + k*pow(zs2,2)) + mua*((1 - p1)*v2am + p1*v2mm)*(1 + k*pow(zs2,2)))))/2.));

vtplus11ma = bound0(v1ma + eul*(pow(exp(1),q*zs1)*s1*(-v1ma + v2am) + ((1 - d)*mum*p2*(v1aa - v1ma)*(1 + k*pow(zs1,2)))/2. - (1 + (-1 + d)/2.)*mum*v1ma*(1 + k*pow(zs1,2)) + mua*(((1 - d)*(1 - p1))/2. + (d*(f1aa*(2 - 2*p1) + f2mm*(2 - 2*p1) + f1mm*(2 - 2*p2) + f2aa*(2 - 2*p2) + f1am*(2 - p1 - p2) + f2am*(2 - p1 - p2)))/2.)*(-v1ma + v1mm)*(1 + k*pow(zs1,2)) + ((1 - d)*mua*(p2*v1am + (1 - p2)*(-v1ma + 2*v1mm))*(1 + k*pow(zs1,2)))/2. + (d*(2*f1aa*mua*(p2*v1aa + (1 - p2)*v1ma)*(1 + k*pow(zs1,2)) + 2*f1mm*mum*(p2*v1am + (1 - p2)*v1mm)*(1 + k*pow(zs1,2)) + f1am*(mum*(p2*v1aa + (1 - p2)*v1ma)*(1 + k*pow(zs1,2)) + mua*(p2*v1am + (1 - p2)*v1mm)*(1 + k*pow(zs1,2))) + 2*f2aa*mua*((1 - p2)*v2aa + p2*v2ma)*(1 + k*pow(zs2,2)) + 2*f2mm*mum*((1 - p2)*v2am + p2*v2mm)*(1 + k*pow(zs2,2)) + f2am*(mum*((1 - p2)*v2aa + p2*v2ma)*(1 + k*pow(zs2,2)) + mua*((1 - p2)*v2am + p2*v2mm)*(1 + k*pow(zs2,2)))))/2.));

vtplus11mm = bound0(v1mm + eul*(pow(exp(1),q*zs1)*s1*(-v1mm + v2aa) + ((1 - d)*mum*p2*(v1am - v1mm)*(1 + k*pow(zs1,2)))/2. + (mum*((1 - d)*p2 + d*(2*f1aa*p1 + 2*f2mm*p1 + 2*f1mm*p2 + 2*f2aa*p2 + f1am*(p1 + p2) + f2am*(p1 + p2)))*(v1ma - v1mm)*(1 + k*pow(zs1,2)))/2. - (1 + (-1 + d)/2.)*mum*v1mm*(1 + k*pow(zs1,2)) + ((1 - d)*mum*(p2*(v1am + v1ma - v1mm) + (1 - p2)*v1mm)*(1 + k*pow(zs1,2)))/2. + (d*(2*f1aa*mua*(p2*v1aa + (1 - p2)*v1ma)*(1 + k*pow(zs1,2)) + 2*f1mm*mum*(p2*v1am + (1 - p2)*v1mm)*(1 + k*pow(zs1,2)) + f1am*(mum*(p2*v1aa + (1 - p2)*v1ma)*(1 + k*pow(zs1,2)) + mua*(p2*v1am + (1 - p2)*v1mm)*(1 + k*pow(zs1,2))) + 2*f2aa*mua*((1 - p2)*v2aa + p2*v2ma)*(1 + k*pow(zs2,2)) + 2*f2mm*mum*((1 - p2)*v2am + p2*v2mm)*(1 + k*pow(zs2,2)) + f2am*(mum*((1 - p2)*v2aa + p2*v2ma)*(1 + k*pow(zs2,2)) + mua*((1 - p2)*v2am + p2*v2mm)*(1 + k*pow(zs2,2)))))/2.));

vtplus12aa = bound0(v2aa + eul*(pow(exp(1),q*zs2)*s2*(v1mm - v2aa) - (1 + (-1 + d)/2.)*mua*v2aa*(1 + k*pow(zs2,2)) + (mua*((1 - d)*p2 + d*(2*f1aa*p1 + 2*f2mm*p1 + 2*f1mm*p2 + 2*f2aa*p2 + f1am*(p1 + p2) + f2am*(p1 + p2)))*(-v2aa + v2am)*(1 + k*pow(zs2,2)))/2. + ((1 - d)*mua*p2*(-v2aa + v2ma)*(1 + k*pow(zs2,2)))/2. + ((1 - d)*mua*((1 - p2)*v2aa + p2*(-v2aa + v2am + v2ma))*(1 + k*pow(zs2,2)))/2. + (d*(2*f1aa*mua*(p2*v1aa + (1 - p2)*v1ma)*(1 + k*pow(zs1,2)) + 2*f1mm*mum*(p2*v1am + (1 - p2)*v1mm)*(1 + k*pow(zs1,2)) + f1am*(mum*(p2*v1aa + (1 - p2)*v1ma)*(1 + k*pow(zs1,2)) + mua*(p2*v1am + (1 - p2)*v1mm)*(1 + k*pow(zs1,2))) + 2*f2aa*mua*((1 - p2)*v2aa + p2*v2ma)*(1 + k*pow(zs2,2)) + 2*f2mm*mum*((1 - p2)*v2am + p2*v2mm)*(1 + k*pow(zs2,2)) + f2am*(mum*((1 - p2)*v2aa + p2*v2ma)*(1 + k*pow(zs2,2)) + mua*((1 - p2)*v2am + p2*v2mm)*(1 + k*pow(zs2,2)))))/2.));

vtplus12am = bound0(v2am + eul*(pow(exp(1),q*zs2)*s2*(v1ma - v2am) + mum*(((1 - d)*(1 - p1))/2. + (d*(f1aa*(2 - 2*p1) + f2mm*(2 - 2*p1) + f1mm*(2 - 2*p2) + f2aa*(2 - 2*p2) + f1am*(2 - p1 - p2) + f2am*(2 - p1 - p2)))/2.)*(v2aa - v2am)*(1 + k*pow(zs2,2)) - (1 + (-1 + d)/2.)*mua*v2am*(1 + k*pow(zs2,2)) + ((1 - d)*mum*((1 - p2)*(2*v2aa - v2am) + p2*v2ma)*(1 + k*pow(zs2,2)))/2. + ((1 - d)*mua*p2*(-v2am + v2mm)*(1 + k*pow(zs2,2)))/2. + (d*(2*f1aa*mua*(p2*v1aa + (1 - p2)*v1ma)*(1 + k*pow(zs1,2)) + 2*f1mm*mum*(p2*v1am + (1 - p2)*v1mm)*(1 + k*pow(zs1,2)) + f1am*(mum*(p2*v1aa + (1 - p2)*v1ma)*(1 + k*pow(zs1,2)) + mua*(p2*v1am + (1 - p2)*v1mm)*(1 + k*pow(zs1,2))) + 2*f2aa*mua*((1 - p2)*v2aa + p2*v2ma)*(1 + k*pow(zs2,2)) + 2*f2mm*mum*((1 - p2)*v2am + p2*v2mm)*(1 + k*pow(zs2,2)) + f2am*(mum*((1 - p2)*v2aa + p2*v2ma)*(1 + k*pow(zs2,2)) + mua*((1 - p2)*v2am + p2*v2mm)*(1 + k*pow(zs2,2)))))/2.));

vtplus12ma = bound0(v2ma + eul*(pow(exp(1),q*zs2)*s2*(v1am - v2ma) + ((1 - d)*mum*(1 - p1)*(v2aa - v2ma)*(1 + k*pow(zs2,2)))/2. - (1 + (-1 + d)/2.)*mum*v2ma*(1 + k*pow(zs2,2)) + mua*(((1 - d)*p2)/2. + (d*(2*f1aa*p1 + 2*f2mm*p1 + 2*f1mm*p2 + 2*f2aa*p2 + f1am*(p1 + p2) + f2am*(p1 + p2)))/2.)*(-v2ma + v2mm)*(1 + k*pow(zs2,2)) + ((1 - d)*mua*((1 - p1)*v2am + p1*(-v2ma + 2*v2mm))*(1 + k*pow(zs2,2)))/2. + (d*(2*f1aa*mua*(p1*v1aa + (1 - p1)*v1ma)*(1 + k*pow(zs1,2)) + 2*f1mm*mum*(p1*v1am + (1 - p1)*v1mm)*(1 + k*pow(zs1,2)) + f1am*(mum*(p1*v1aa + (1 - p1)*v1ma)*(1 + k*pow(zs1,2)) + mua*(p1*v1am + (1 - p1)*v1mm)*(1 + k*pow(zs1,2))) + 2*f2aa*mua*((1 - p1)*v2aa + p1*v2ma)*(1 + k*pow(zs2,2)) + 2*f2mm*mum*((1 - p1)*v2am + p1*v2mm)*(1 + k*pow(zs2,2)) + f2am*(mum*((1 - p1)*v2aa + p1*v2ma)*(1 + k*pow(zs2,2)) + mua*((1 - p1)*v2am + p1*v2mm)*(1 + k*pow(zs2,2)))))/2.));

vtplus12mm = bound0((2 - 2*f1aa*vtplus11aa - f1am*vtplus11am - f1am*vtplus11ma - 2*f1mm*vtplus11mm - 2*f2aa*vtplus12aa - f2am*vtplus12am - f2am*vtplus12ma)/(2.*f2mm));

if (

fabs(vtplus11aa - v1aa)<=1e-07&&

fabs(vtplus11am - v1am)<=1e-07&&

fabs(vtplus11ma - v1ma)<=1e-07&&

fabs(vtplus11mm - v1mm)<=1e-07&&

fabs(vtplus12aa - v2aa)<=1e-07&&

fabs(vtplus12am - v2am)<=1e-07&&

fabs(vtplus12ma - v2ma)<=1e-07&&

fabs(vtplus12mm - v2mm)<=1e-07)

{break;}

v1aa = vtplus11aa;

v1am = vtplus11am;

v1ma = vtplus11ma;

v1mm = vtplus11mm;

v2aa = vtplus12aa;

v2am = vtplus12am;

v2ma = vtplus12ma;

v2mm = vtplus12mm;

}

gsl_vector_set(f,0,v1aa);
gsl_vector_set(f,1,v1am);
gsl_vector_set(f,2,v1ma);
gsl_vector_set(f,3,v1mm);
gsl_vector_set(f,4,v2aa);
gsl_vector_set(f,5,v2am);
gsl_vector_set(f,6,v2ma);
gsl_vector_set(f,7,v2mm);

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
double q = ((struct rparams *) params)->q;
double eul = ((struct rparams *) params)->eul;
double p1 = ((struct rparams *) params)->p1;
double p2 = ((struct rparams *) params)->p2;
double zs1 = ((struct rparams *) params)->zs1;
double zs2 = ((struct rparams *) params)->zs2;
double f1aa = ((struct rparams *) params)->f1aa;
double f1am = ((struct rparams *) params)->f1am;
double f1mm = ((struct rparams *) params)->f1mm;
double f2aa = ((struct rparams *) params)->f2aa;
double f2am = ((struct rparams *) params)->f2am;
double f2mm = ((struct rparams *) params)->f2mm;
double v1aa = ((struct rparams *) params)->v1aa;
double v1am = ((struct rparams *) params)->v1am;
double v1ma = ((struct rparams *) params)->v1ma;
double v1mm = ((struct rparams *) params)->v1mm;
double v2aa = ((struct rparams *) params)->v2aa;
double v2am = ((struct rparams *) params)->v2am;
double v2ma = ((struct rparams *) params)->v2ma;
double v2mm = ((struct rparams *) params)->v2mm;
double raa1 = ((struct rparams *) params)->raa1;
double ram1 = ((struct rparams *) params)->ram1;
double rmm1 = ((struct rparams *) params)->rmm1;
double raa2 = ((struct rparams *) params)->raa2;
double ram2 = ((struct rparams *) params)->ram2;
double rmm2 = ((struct rparams *) params)->rmm2;


    // 
double p1tplus1 = bound01(p1 + eul * (-((1 - d)*f1am*mua*ram1*(-v1ma + v1mm)*(1 + k*pow(zs1,2)))/2. - ((1 - d)*f2am*mum*ram2*(v2aa - v2am)*(1 + k*pow(zs2,2)))/2. + 2*f1aa*(-((1 - d)*mua*raa1*(-v1aa + v1am)*(1 + k*pow(zs1,2)))/2. + ((1 - d)*mua*(2*v1aa - v1am - v1ma)*(1 + k*pow(zs1,2)))/2. - ((1 - d)*mua*(-v1aa + v1ma)*(1 + k*pow(zs1,2)))/2. + (d*(2*f1aa*mua*(v1aa - v1ma)*(1 + k*pow(zs1,2)) + 2*f1mm*mum*(v1am - v1mm)*(1 + k*pow(zs1,2)) + f1am*(mum*(v1aa - v1ma)*(1 + k*pow(zs1,2)) + mua*(v1am - v1mm)*(1 + k*pow(zs1,2))) + 2*f2aa*mua*(-v2aa + v2ma)*(1 + k*pow(zs2,2)) + 2*f2mm*mum*(-v2am + v2mm)*(1 + k*pow(zs2,2)) + f2am*(mum*(-v2aa + v2ma)*(1 + k*pow(zs2,2)) + mua*(-v2am + v2mm)*(1 + k*pow(zs2,2)))))/2.) + f1am*(((1 - d)*mum*(2*v1aa - v1am - v1ma)*(1 + k*pow(zs1,2)))/2. - ((1 - d)*mua*(-v1am + v1mm)*(1 + k*pow(zs1,2)))/2. + (d*(2*f1aa*mua*(v1aa - v1ma)*(1 + k*pow(zs1,2)) + 2*f1mm*mum*(v1am - v1mm)*(1 + k*pow(zs1,2)) + f1am*(mum*(v1aa - v1ma)*(1 + k*pow(zs1,2)) + mua*(v1am - v1mm)*(1 + k*pow(zs1,2))) + 2*f2aa*mua*(-v2aa + v2ma)*(1 + k*pow(zs2,2)) + 2*f2mm*mum*(-v2am + v2mm)*(1 + k*pow(zs2,2)) + f2am*(mum*(-v2aa + v2ma)*(1 + k*pow(zs2,2)) + mua*(-v2am + v2mm)*(1 + k*pow(zs2,2)))))/2.) + f2am*(-((1 - d)*mum*(v2aa - v2ma)*(1 + k*pow(zs2,2)))/2. + ((1 - d)*mua*(-v2am - v2ma + 2*v2mm)*(1 + k*pow(zs2,2)))/2. + (d*(2*f1aa*mua*(v1aa - v1ma)*(1 + k*pow(zs1,2)) + 2*f1mm*mum*(v1am - v1mm)*(1 + k*pow(zs1,2)) + f1am*(mum*(v1aa - v1ma)*(1 + k*pow(zs1,2)) + mua*(v1am - v1mm)*(1 + k*pow(zs1,2))) + 2*f2aa*mua*(-v2aa + v2ma)*(1 + k*pow(zs2,2)) + 2*f2mm*mum*(-v2am + v2mm)*(1 + k*pow(zs2,2)) + f2am*(mum*(-v2aa + v2ma)*(1 + k*pow(zs2,2)) + mua*(-v2am + v2mm)*(1 + k*pow(zs2,2)))))/2.) + 2*f2mm*(-((1 - d)*mum*(v2am - v2mm)*(1 + k*pow(zs2,2)))/2. - ((1 - d)*mum*rmm2*(v2ma - v2mm)*(1 + k*pow(zs2,2)))/2. + ((1 - d)*mum*(-v2am - v2ma + 2*v2mm)*(1 + k*pow(zs2,2)))/2. + (d*(2*f1aa*mua*(v1aa - v1ma)*(1 + k*pow(zs1,2)) + 2*f1mm*mum*(v1am - v1mm)*(1 + k*pow(zs1,2)) + f1am*(mum*(v1aa - v1ma)*(1 + k*pow(zs1,2)) + mua*(v1am - v1mm)*(1 + k*pow(zs1,2))) + 2*f2aa*mua*(-v2aa + v2ma)*(1 + k*pow(zs2,2)) + 2*f2mm*mum*(-v2am + v2mm)*(1 + k*pow(zs2,2)) + f2am*(mum*(-v2aa + v2ma)*(1 + k*pow(zs2,2)) + mua*(-v2am + v2mm)*(1 + k*pow(zs2,2)))))/2.)));

double p2tplus1 = bound01(p2 + eul * (((1 - d)*f1am*mum*ram1*(v1aa - v1am)*(1 + k*pow(zs1,2)))/2. + ((1 - d)*f2am*mua*ram2*(-v2ma + v2mm)*(1 + k*pow(zs2,2)))/2. + f1am*(((1 - d)*mum*(v1aa - v1ma)*(1 + k*pow(zs1,2)))/2. + ((1 - d)*mua*(v1am + v1ma - 2*v1mm)*(1 + k*pow(zs1,2)))/2. + (d*(2*f1aa*mua*(v1aa - v1ma)*(1 + k*pow(zs1,2)) + 2*f1mm*mum*(v1am - v1mm)*(1 + k*pow(zs1,2)) + f1am*(mum*(v1aa - v1ma)*(1 + k*pow(zs1,2)) + mua*(v1am - v1mm)*(1 + k*pow(zs1,2))) + 2*f2aa*mua*(-v2aa + v2ma)*(1 + k*pow(zs2,2)) + 2*f2mm*mum*(-v2am + v2mm)*(1 + k*pow(zs2,2)) + f2am*(mum*(-v2aa + v2ma)*(1 + k*pow(zs2,2)) + mua*(-v2am + v2mm)*(1 + k*pow(zs2,2)))))/2.) + 2*f1mm*(((1 - d)*mum*(v1am + v1ma - 2*v1mm)*(1 + k*pow(zs1,2)))/2. + ((1 - d)*mum*(v1am - v1mm)*(1 + k*pow(zs1,2)))/2. + ((1 - d)*mum*rmm1*(v1ma - v1mm)*(1 + k*pow(zs1,2)))/2. + (d*(2*f1aa*mua*(v1aa - v1ma)*(1 + k*pow(zs1,2)) + 2*f1mm*mum*(v1am - v1mm)*(1 + k*pow(zs1,2)) + f1am*(mum*(v1aa - v1ma)*(1 + k*pow(zs1,2)) + mua*(v1am - v1mm)*(1 + k*pow(zs1,2))) + 2*f2aa*mua*(-v2aa + v2ma)*(1 + k*pow(zs2,2)) + 2*f2mm*mum*(-v2am + v2mm)*(1 + k*pow(zs2,2)) + f2am*(mum*(-v2aa + v2ma)*(1 + k*pow(zs2,2)) + mua*(-v2am + v2mm)*(1 + k*pow(zs2,2)))))/2.) + 2*f2aa*(((1 - d)*mua*raa2*(-v2aa + v2am)*(1 + k*pow(zs2,2)))/2. + ((1 - d)*mua*(-v2aa + v2ma)*(1 + k*pow(zs2,2)))/2. + ((1 - d)*mua*(-2*v2aa + v2am + v2ma)*(1 + k*pow(zs2,2)))/2. + (d*(2*f1aa*mua*(v1aa - v1ma)*(1 + k*pow(zs1,2)) + 2*f1mm*mum*(v1am - v1mm)*(1 + k*pow(zs1,2)) + f1am*(mum*(v1aa - v1ma)*(1 + k*pow(zs1,2)) + mua*(v1am - v1mm)*(1 + k*pow(zs1,2))) + 2*f2aa*mua*(-v2aa + v2ma)*(1 + k*pow(zs2,2)) + 2*f2mm*mum*(-v2am + v2mm)*(1 + k*pow(zs2,2)) + f2am*(mum*(-v2aa + v2ma)*(1 + k*pow(zs2,2)) + mua*(-v2am + v2mm)*(1 + k*pow(zs2,2)))))/2.) + f2am*(((1 - d)*mum*(-2*v2aa + v2am + v2ma)*(1 + k*pow(zs2,2)))/2. + ((1 - d)*mua*(-v2am + v2mm)*(1 + k*pow(zs2,2)))/2. + (d*(2*f1aa*mua*(v1aa - v1ma)*(1 + k*pow(zs1,2)) + 2*f1mm*mum*(v1am - v1mm)*(1 + k*pow(zs1,2)) + f1am*(mum*(v1aa - v1ma)*(1 + k*pow(zs1,2)) + mua*(v1am - v1mm)*(1 + k*pow(zs1,2))) + 2*f2aa*mua*(-v2aa + v2ma)*(1 + k*pow(zs2,2)) + 2*f2mm*mum*(-v2am + v2mm)*(1 + k*pow(zs2,2)) + f2am*(mum*(-v2aa + v2ma)*(1 + k*pow(zs2,2)) + mua*(-v2am + v2mm)*(1 + k*pow(zs2,2)))))/2.)));

double zs1tplus1 = zs1 + eul * (f1am*((pow(exp(1),q*zs1)*q*s1*(-v1ma + v2am))/2. + (1 - d)*k*mum*p2*(v1aa - v1ma)*zs1 - 2*(1 + (-1 + d)/2.)*k*mum*v1ma*zs1) + f1am*((pow(exp(1),q*zs1)*q*ram1*s1*(-v1am + v2ma))/2. + 2*k*mum*(((1 - d)*p2)/2. + (d*(2*f1aa*p1 + 2*f2mm*p1 + 2*f1mm*p2 + 2*f2aa*p2 + f1am*(p1 + p2) + f2am*(p1 + p2)))/2.)*ram1*(v1aa - v1am)*zs1 + (1 - d)*k*mum*ram1*(p1*(2*v1aa - v1am) + (1 - p1)*v1ma)*zs1) + 2*f1aa*((pow(exp(1),q*zs1)*q*(1 + raa1)*s1*(-v1aa + v2mm))/2. - 2*(1 + (-1 + d)/2.)*k*mua*v1aa*zs1 + k*mua*((1 - d)*(1 - p1) + d*(f1aa*(2 - 2*p1) + f2mm*(2 - 2*p1) + f1mm*(2 - 2*p2) + f2aa*(2 - 2*p2) + f1am*(2 - p1 - p2) + f2am*(2 - p1 - p2)))*raa1*(-v1aa + v1am)*zs1 + (1 - d)*k*mua*(1 - p1)*(-v1aa + v1ma)*zs1 + (1 - d)*k*mua*raa1*(p1*v1aa + (1 - p1)*(-v1aa + v1am + v1ma))*zs1) + f1am*((pow(exp(1),q*zs1)*q*s1*(-v1am + v2ma))/2. - 2*(1 + (-1 + d)/2.)*k*mua*v1am*zs1 + (1 - d)*k*mua*(1 - p1)*(-v1am + v1mm)*zs1) + 2*f1mm*((pow(exp(1),q*zs1)*q*(1 + rmm1)*s1*(-v1mm + v2aa))/2. + (1 - d)*k*mum*p2*(v1am - v1mm)*zs1 + k*mum*((1 - d)*p2 + d*(2*f1aa*p1 + 2*f2mm*p1 + 2*f1mm*p2 + 2*f2aa*p2 + f1am*(p1 + p2) + f2am*(p1 + p2)))*rmm1*(v1ma - v1mm)*zs1 - 2*(1 + (-1 + d)/2.)*k*mum*v1mm*zs1 + (1 - d)*k*mum*rmm1*(p2*(v1am + v1ma - v1mm) + (1 - p2)*v1mm)*zs1) + f1am*((pow(exp(1),q*zs1)*q*ram1*s1*(-v1ma + v2am))/2. + 2*k*mua*(((1 - d)*(1 - p1))/2. + (d*(f1aa*(2 - 2*p1) + f2mm*(2 - 2*p1) + f1mm*(2 - 2*p2) + f2aa*(2 - 2*p2) + f1am*(2 - p1 - p2) + f2am*(2 - p1 - p2)))/2.)*ram1*(-v1ma + v1mm)*zs1 + (1 - d)*k*mua*ram1*(p2*v1am + (1 - p2)*(-v1ma + 2*v1mm))*zs1));

double zs2tplus1 = zs2 + eul * (f2am*((pow(exp(1),q*zs2)*q*s2*(v1am - v2ma))/2. + (1 - d)*k*mum*(1 - p1)*(v2aa - v2ma)*zs2 - 2*(1 + (-1 + d)/2.)*k*mum*v2ma*zs2) + f2am*((pow(exp(1),q*zs2)*q*ram2*s2*(v1ma - v2am))/2. + 2*k*mum*(((1 - d)*(1 - p1))/2. + (d*(f1aa*(2 - 2*p1) + f2mm*(2 - 2*p1) + f1mm*(2 - 2*p2) + f2aa*(2 - 2*p2) + f1am*(2 - p1 - p2) + f2am*(2 - p1 - p2)))/2.)*ram2*(v2aa - v2am)*zs2 + (1 - d)*k*mum*ram2*((1 - p2)*(2*v2aa - v2am) + p2*v2ma)*zs2) + 2*f2aa*((pow(exp(1),q*zs2)*q*(1 + raa2)*s2*(v1mm - v2aa))/2. - 2*(1 + (-1 + d)/2.)*k*mua*v2aa*zs2 + k*mua*((1 - d)*p2 + d*(2*f1aa*p1 + 2*f2mm*p1 + 2*f1mm*p2 + 2*f2aa*p2 + f1am*(p1 + p2) + f2am*(p1 + p2)))*raa2*(-v2aa + v2am)*zs2 + (1 - d)*k*mua*p2*(-v2aa + v2ma)*zs2 + (1 - d)*k*mua*raa2*((1 - p2)*v2aa + p2*(-v2aa + v2am + v2ma))*zs2) + f2am*((pow(exp(1),q*zs2)*q*s2*(v1ma - v2am))/2. - 2*(1 + (-1 + d)/2.)*k*mua*v2am*zs2 + (1 - d)*k*mua*p2*(-v2am + v2mm)*zs2) + 2*f2mm*((pow(exp(1),q*zs2)*q*(1 + rmm2)*s2*(v1aa - v2mm))/2. + (1 - d)*k*mum*(1 - p1)*(v2am - v2mm)*zs2 + k*mum*((1 - d)*(1 - p1) + d*(f1aa*(2 - 2*p1) + f2mm*(2 - 2*p1) + f1mm*(2 - 2*p2) + f2aa*(2 - 2*p2) + f1am*(2 - p1 - p2) + f2am*(2 - p1 - p2)))*rmm2*(v2ma - v2mm)*zs2 - 2*(1 + (-1 + d)/2.)*k*mum*v2mm*zs2 + (1 - d)*k*mum*rmm2*((1 - p1)*(v2am + v2ma - v2mm) + p1*v2mm)*zs2) + f2am*((pow(exp(1),q*zs2)*q*ram2*s2*(v1am - v2ma))/2. + 2*k*mua*(((1 - d)*p2)/2. + (d*(2*f1aa*p1 + 2*f2mm*p1 + 2*f1mm*p2 + 2*f2aa*p2 + f1am*(p1 + p2) + f2am*(p1 + p2)))/2.)*ram2*(-v2ma + v2mm)*zs2 + (1 - d)*k*mua*ram2*((1 - p1)*v2am + p1*(-v2ma + 2*v2mm))*zs2));

gsl_vector_set(f, 0, p1tplus1);
gsl_vector_set(f, 1, p2tplus1);
gsl_vector_set(f, 2, zs1tplus1);
gsl_vector_set(f, 3, zs2tplus1);


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
double q = ((struct rparams *) params)->q;
double eul = ((struct rparams *) params)->eul;
double p1 = ((struct rparams *) params)->p1;
double p2 = ((struct rparams *) params)->p2;
double zs1 = ((struct rparams *) params)->zs1;
double zs2 = ((struct rparams *) params)->zs2;
double f1aa = ((struct rparams *) params)->f1aa;
double f1am = ((struct rparams *) params)->f1am;
double f1mm = ((struct rparams *) params)->f1mm;
double f2aa = ((struct rparams *) params)->f2aa;
double f2am = ((struct rparams *) params)->f2am;
double f2mm = ((struct rparams *) params)->f2mm;
double v1aa = ((struct rparams *) params)->v1aa;
double v1am = ((struct rparams *) params)->v1am;
double v1ma = ((struct rparams *) params)->v1ma;
double v1mm = ((struct rparams *) params)->v1mm;
double v2aa = ((struct rparams *) params)->v2aa;
double v2am = ((struct rparams *) params)->v2am;
double v2ma = ((struct rparams *) params)->v2ma;
double v2mm = ((struct rparams *) params)->v2mm;
double raa1 = ((struct rparams *) params)->raa1;
double ram1 = ((struct rparams *) params)->ram1;
double rmm1 = ((struct rparams *) params)->rmm1;
double raa2 = ((struct rparams *) params)->raa2;
double ram2 = ((struct rparams *) params)->ram2;
double rmm2 = ((struct rparams *) params)->rmm2;



    // 
DataFile << endl << endl  << "mum;" << mum << endl
 << "mua;" << mua << endl
 << "d;" << d << endl
 << "s1;" << s1 << endl
 << "s2;" << s2 << endl
 << "k;" << k << endl
 << "q;" << q << endl
 << "eul;" << eul << endl
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
double q = ((struct rparams *) params)->q;
double eul = ((struct rparams *) params)->eul;
double p1 = ((struct rparams *) params)->p1;
double p2 = ((struct rparams *) params)->p2;
double zs1 = ((struct rparams *) params)->zs1;
double zs2 = ((struct rparams *) params)->zs2;
double f1aa = ((struct rparams *) params)->f1aa;
double f1am = ((struct rparams *) params)->f1am;
double f1mm = ((struct rparams *) params)->f1mm;
double f2aa = ((struct rparams *) params)->f2aa;
double f2am = ((struct rparams *) params)->f2am;
double f2mm = ((struct rparams *) params)->f2mm;
double v1aa = ((struct rparams *) params)->v1aa;
double v1am = ((struct rparams *) params)->v1am;
double v1ma = ((struct rparams *) params)->v1ma;
double v1mm = ((struct rparams *) params)->v1mm;
double v2aa = ((struct rparams *) params)->v2aa;
double v2am = ((struct rparams *) params)->v2am;
double v2ma = ((struct rparams *) params)->v2ma;
double v2mm = ((struct rparams *) params)->v2mm;
double raa1 = ((struct rparams *) params)->raa1;
double ram1 = ((struct rparams *) params)->ram1;
double rmm1 = ((struct rparams *) params)->rmm1;
double raa2 = ((struct rparams *) params)->raa2;
double ram2 = ((struct rparams *) params)->ram2;
double rmm2 = ((struct rparams *) params)->rmm2;



    if (time < 0)
    {
        // 
DataFile << "time;f1mm;zs2;ram2;v2ma;raa1;v2mm;v1aa;v2am;rmm1;v2aa;zs1;p1;f1am;f2mm;v1mm;f1aa;raa2;v1am;f2am;ram1;f2aa;p2;v1ma;rmm2;" << endl;
    }

    // 
DataFile << time << ";" << f1mm << ";" << 
zs2 << ";" << 
ram2 << ";" << 
v2ma << ";" << 
raa1 << ";" << 
v2mm << ";" << 
v1aa << ";" << 
v2am << ";" << 
rmm1 << ";" << 
v2aa << ";" << 
zs1 << ";" << 
p1 << ";" << 
f1am << ";" << 
f2mm << ";" << 
v1mm << ";" << 
f1aa << ";" << 
raa2 << ";" << 
v1am << ";" << 
f2am << ";" << 
ram1 << ";" << 
f2aa << ";" << 
p2 << ";" << 
v1ma << ";" << 
rmm2 << ";" << 
 endl;
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
    // 
		paramstruct.mum = atof(argv[2]);
		paramstruct.mua = atof(argv[3]);
		paramstruct.d = atof(argv[4]);
		paramstruct.s1 = atof(argv[5]);
		paramstruct.s2 = atof(argv[6]);
		paramstruct.k = atof(argv[7]);
		paramstruct.q = atof(argv[8]);
		paramstruct.eul = atof(argv[9]);
		paramstruct.p1 = atof(argv[10]);
		paramstruct.p2 = atof(argv[11]);
		paramstruct.zs1 = atof(argv[12]);
		paramstruct.zs2 = atof(argv[13]);
		paramstruct.f1aa = atof(argv[14]);
		paramstruct.f1am = atof(argv[15]);
		paramstruct.f1mm = atof(argv[16]);
		paramstruct.f2aa = atof(argv[17]);
		paramstruct.f2am = atof(argv[18]);
		paramstruct.f2mm = atof(argv[19]);
		paramstruct.v1aa = atof(argv[20]);
		paramstruct.v1am = atof(argv[21]);
		paramstruct.v1ma = atof(argv[22]);
		paramstruct.v1mm = atof(argv[23]);
		paramstruct.v2aa = atof(argv[24]);
		paramstruct.v2am = atof(argv[25]);
		paramstruct.v2ma = atof(argv[26]);
		paramstruct.v2mm = atof(argv[27]);
		paramstruct.raa1 = atof(argv[28]);
		paramstruct.ram1 = atof(argv[29]);
		paramstruct.rmm1 = atof(argv[30]);
		paramstruct.raa2 = atof(argv[31]);
		paramstruct.ram2 = atof(argv[32]);
		paramstruct.rmm2 = atof(argv[33]);

   

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

