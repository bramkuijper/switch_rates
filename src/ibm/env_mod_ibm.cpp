// Gillespie simulation 
// for the evolution of nongenetic effects and enviromental modulation
// (c) Bram Kuijper 2016

#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cassert>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "bramauxiliary.h"

//#define NDEBUG

using namespace std;


// random number generator 
// see http://www.gnu.org/software/gsl/manual/html_node/Random-Number-Generation.html#Random-Number-Generation 
gsl_rng_type const * T; // gnu scientific library rng type
gsl_rng *r; // gnu scientific rng 

const size_t Npatch = 4000; 
const size_t Npp = 2;
unsigned long int NumGen = 1e10; // maximum number of generations
int sample = 20;

double p0 = 0.5;

int seed = -1; // the random seed
time_t total_time;  // track the total time the simulation is running
unsigned long int generation = 0; // track the current generation 

// output once every 10^6 generations
unsigned long int skip = 1e6;

double d = 0.1; // dispersal probability
double q = 0.1; // efficacy of environmental modulation
double k = 0.1; // survival cost of environmental modulation
double mum = 2.0; // mortality rate when maladapted
double mua = 2.0; // mortality rate when adapted


// mutation rates
double mu = 0.01;
double sdmu = 0.01;

struct Individual {

    double p1[2];
    double p2[2];
    double zs1[2];
    double zs2[2];

    bool phen;
};

struct Patch {
    Individual locals[Npp];

    bool environment;
};


Patch MetaPop[Npatch];

// distribution of patches
// numbers of patches
// in environmental state ei={1,2}
// having 0, 1, 2 adapted individuals
size_t Npatches[2][3];

// next to counts store lists with the ids
// of the corresponding patches
size_t patch_ids[2][3][Npatch];

// give the outputfile a unique name
string filename("sim_envmod");
string filename_new(create_filename(filename));
ofstream DataFile(filename_new.c_str());  // output file 

// initialize the command line arguments to vary 
// the parameters
void init_arguments(int argc, char *argv[])
{
    mu = atof(argv[1]);
    sdmu = atof(argv[2]);
    d = atof(argv[3]);
    q = atof(argv[4]);
    k = atof(argv[5]);
    p0 = atof(argv[6]);
}

// initialization function, runs at start of sim
// gives all the individuals some starting
// value
void init_pop()
{ 
    // start the time
    total_time = time(NULL);

    // obtain a seed from current nanosecond count
	seed = get_nanoseconds();

    // set up the random number generators
    // (from the gnu gsl library)
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, seed);

    // initialize patch counters and set them to 0
    for (int n_envt = 0; n_envt < 2; ++n_envt)
    {
        for (int n_adapted = 0; n_adapted <= Npp; ++n_adapted)
        {
            Npatches[n_envt][n_adapted] = 0;
        }
    }

    size_t n_adapted;

    // initialize all the patches
    for (size_t i = 0; i < Npatch; ++i)
    {
        // keep track of the number of adapted
        // individuals in this patch
        n_adapted = 0;

        // randomly assign environment
        MetaPop[i].environment = gsl_rng_uniform(r) < 0.5;

        //initialize breeders
        for (size_t j = 0; j < Npp; ++j)
        {
            // randomly assign initial phenotypes
            MetaPop[i].locals[j].phen = gsl_rng_uniform(r) < 0.5;

            if (MetaPop[i].locals[j].phen == MetaPop[i].environment)
            {
                ++n_adapted;
            }

            // initialize allelic values
            for (int allele = 0; allele < 2; ++allele)
            {
                MetaPop[i].locals[j].p1[allele] = p0;
                MetaPop[i].locals[j].p2[allele] = p0;
                MetaPop[i].locals[j].zs1[allele] = 0;
                MetaPop[i].locals[j].zs2[allele] = 0;
            }
        }

        // increase the counter of patches 
        // in environmental state ei, 
        // containing n_adapted individuals
        //
        size_t Npatches_of_this_type = 
            Npatches[MetaPop[i].environment][n_adapted];

        // append this patch to the corresponding stack of patch id's
        patch_ids[MetaPop[i].environment][n_adapted][Npatches_of_this_type] = i;

        // update counter
        ++Npatches[MetaPop[i].environment][n_adapted];
    }
}

// mutation for p1,p2 (bounded between 0 and 1)
double mut_pi(double p)
{
    assert(p >= 0 && p <= 1);
    p+= gsl_rng_uniform(r) < mu ? gsl_ran_gaussian(r,sdmu_h) : 0;
    p = p < 0 ? 0 : p > 1 ? 1 : p;

    return(p);
}

// mutation of the environmental modulation traits
double mut_zsi(double zsi)
{
    zsi+= gsl_rng_uniform(r) < mu ? gsl_ran_gaussian(r,sdmu_h) : 0;
    return(zsi);
}

// make a kid out of parent asexually
void create_kid(Individual &mother, Individual &Kid)
{
    // if mother.phen == true, then mother bears
    // z2 phenotype, hence we have to use p2
    // otherwise she bears a z1 phenotype and we
    // have to use p1;
    double expr_p = mother.phen ? 
        .5 * (mother.p2[0] + mother.p2[1])
        :
        .5 * (mother.p1[0] + mother.p1[1]);

    // offspring attains z1 phenotype (i.e., 0 or false) 
    // with probability expr_p. Otherwise it will be z2
    Kid.hawk = gsl_rng_uniform(r) < expr_p ? 0 : 1;

    // inherit alleles
    for (int allele_i = 0; allele_i < 2; ++allele_i)
    {
        Kid.p1[allele_i] = mut_pi(mother.p1[allele_i]);
        Kid.p2[allele_i] = mut_pi(mother.p2[allele_i]);

        Kid.zs1[allele_i] = mut_zsi(mother.zs1[allele_i]);
        Kid.zs2[allele_i] = mut_zsi(mother.zs2[allele_i]);
    }
}

// do the stats on the data
void write_data()
{
    // get variance and means
    double p1, p2, zs1, zs2;
    double varp1, varp2, varzs1, varzs2;
    double meanp1, meanp2, meanzs1, meanzs2;
    double ssp1, ssp2, sszs1, sszs2;

    meanp1 = 0;
    meanp2 = 0;
    meanzs1 = 0;
    meanzs2 = 0;
    ssp1 = 0;
    ssp2 = 0;
    sszs1 = 0;
    sszs2 = 0;

    // loop through patches and individuals
    // and get stats on genotypes and patch freqs.
    for (size_t i = 0; i < Npatch; ++i)
    {
        for (size_t j = 0; j < Npp; ++j)
        {
            p1 = 0.5 * (MetaPop[i].locals[j].p1[0] + MetaPop[i].locals[j].p1[1]);
            p2 = 0.5 * (MetaPop[i].locals[j].p2[0] + MetaPop[i].locals[j].p2[1]);
            zs1 = 0.5 * (MetaPop[i].locals[j].zs1[0] + MetaPop[i].locals[j].zs1[1]);
            zs2 = 0.5 * (MetaPop[i].locals[j].zs2[0] + MetaPop[i].locals[j].zs2[1]);

            meanp1+=p1;
            meanp2+=p2;
            meanzs1+=zs1;
            meanzs2+=zs2;

            ssp1+=p1*p1;
            ssp2+=p2*p2;
            sszs1+=zs1*zs1;
            sszs2+=zs2*zs2;
        }    
    }

    DataFile << generation << ";"; 

    // write down all the patch frequency combinations
    for (size_t n_envt = 0; n_envt < 2; ++n_envt)
    {
        for (size_t n_adapted = 0; n_adapted <= Npp; n_adapted++)
        {
            DataFile << setprecision(5) 
                << (double) Npatches[n_envt][n_adapted] / Npatch << ";";
        }
    }

    meanp1/=Npp * Npatch;
    meanp2/=Npp * Npatch;
    meanzs1/=Npp * Npatch;
    meanzs2/=Npp * Npatch;

    varp1 = ssp1 / (Npatch * Npp) - meanp1 *meanp1;
    varp2 = ssp2 / (Npatch * Npp) - meanp2 *meanp2;
    varzs1 = sszs1 / (Npatch * Npp) - meanzs1 *meanzs1;
    varzs2 = sszs2 / (Npatch * Npp) - meanzs2 *meanzs2;

    DataFile 
        << setprecision(5) << meanp1 << ";"
        << setprecision(5) << meanp2 << ";"
        << setprecision(5) << meanzs1 << ";"
        << setprecision(5) << meanzs2 << ";"
        << setprecision(5) << varp1 << ";"
        << setprecision(5) << varp2 << ";"
        << setprecision(5) << varzs1 << ";" 
        << setprecision(5) << varzs2 << ";" << endl;
}

// headers of the datafile
void write_data_headers()
{
    DataFile << "generation;";

    for (size_t n_envt = 0; n_envt < 2; ++n_envt)
    {
        for (size_t n_adapted = 0; n_adapted <= Npp; n_adapted++)
        {
            DataFile << "f" << (n_envt + 1) << n_adapted << ";";
        }
    }
    
    DataFile << "meanp1;meanp2;meanzs1;meanzs2;varp1;varp2;varzs1;varzs2;" << endl;
}

// parameters at the end of the sim
void write_parameters()
{
    total_time = time(NULL) - total_time;
    DataFile << endl << endl 
                << "d;" << d << endl
                << "k;" << k << endl
                << "q;" << q << endl
                << "p0;" << p0 << endl
                << "mu;" << mu << endl
                << "sdmu;" << sdmu << endl
                << "NumGen;" << NumGen << endl
                << "seed;" << seed << endl
                << "Npatches;" << Npatch << endl
                << "runtime;" << total_time << endl;
}


void mortality(size_t n_adapted, bool environment, bool mortality_adapted)
{
    assert(n_adapted == 0 && !mortality_adapted);
    assert(n_adapted == 2 && mortality_adapted);
    
    // sample a random patch that fulfills the conditions
    size_t random_patch = 
        gsl_rng_uniform_int(r, 
            Npatches[environment][n_adapted]
        );

    // get the patch_id 
    size_t random_patch_id = 
        patch_ids[environment][n_adapted][random_patch];

    // error checking
    assert(MetaPop[random_patch_id].environment == environment);

#ifndef NDEBUG

    size_t n_adapted_check = 0;

    for (size_t breeder_i = 0; breeder_i < Npp; ++breeder_i)
    {
        n_adapted_check += MetaPop[random_patch_id].locals[breeder_i].phen == environment; 
    }

    assert(n_adapted_check == n_adapted);
#endif

    // array with the individuals that may die
    size_t mortality_candidate_ids[Npp];
    size_t n_mortality_candidates = 0;

    // at the same time calculate which local breeders
    // might die
    for (size_t local_i = 0; local_i < Npp; ++local_i)
    {
        // store the mortality candidates
        //
        // mortality adapted
        if (
                (
                 mortality_adapted && 
                    MetaPop[random_patch_id].locals[local_i].phen == environment
                )
                ||
                (!mortality_adapted && 
                    MetaPop[random_patch_id].locals[local_i].phen != environment
                )

            )
        {
            mortality_candidate_ids[n_mortality_candidates++] = local_i;
        }
    }


    // now perform mortality and sample recruit

    // keep track of new state variables
    // to update patch type statistics and counters (see below)
    size_t n_adapted_new = n_adapted;

    // then mortality 
    assert(n_mortality_candidates > 0 && n_mortality_candidates <= Npp);

    // sample a random candidate for mortality
    size_t candidate_id = 
        mortality_candidate_ids[gsl_rng_uniform_int(r, n_mortality_candidates)];

    assert(candidate_id >= 0 && candidate_id < Npp);

    assert(
            (mortality_adapted &&
            MetaPop[random_patch_id].locals[candidate_id].phen == environment
            )
            ||
            (!mortality_adapted &&
            !MetaPop[random_patch_id].locals[candidate_id].phen != environment
            )
    );

    if (mortality_adapted)
    {
        --n_adapted_new;
    }


    // make new Kid
    Individual Kid;

    size_t random_remote_patch_id = 0;

    // sample local kid
    if (gsl_rng_uniform(r) < 1.0 - d)
    {
        create_kid(
                MetaPop[random_patch_id].locals[gsl_rng_uniform_int(r, Npp)],
                Kid
                );
    }
    else
    {
        // birth from a remote parent
        random_remote_patch_id = gsl_rng_uniform_int(r, Npatch);

        create_kid(
                MetaPop[random_remote_patch_id].locals[gsl_rng_uniform_int(r, Npp)], 
                Kid
                );
    }

    // mortality and replacement
    MetaPop[random_patch_id].locals[candidate_id] = Kid;

    // if newborn is hawk, increment count
    if (MetaPop[random_patch_id].locals[candidate_id].phen == environment)    
    {
        ++n_adapted_new;
    }

    // update statistics
    //
    // delete patch id in the corresponding stack
    // by replacing it with the last patch id in the stack
    // and then deleting the last element (by reducing the counter by 1)
    patch_ids[environment][n_adapted][random_patch] =
        patch_ids[environment][n_adapted][
            --Npatches[environment][n_adapted] 
        ];

    // add patch id to the correct stack
    patch_ids[environment][n_adapted_new][
        Npatches[environment][n_adapted_new]++
    ] = random_patch_id;

    // we're done.
    // error checking
#ifndef NDEBUG

    size_t n_adapted_check = 0;

    for (size_t breeder_i = 0; breeder_i < Npp; ++breeder_i)
    {
        n_adapted_check += MetaPop[random_patch_id].locals[breeder_i].phen == environment; 
    }

    assert(n_adapted_check == n_adapted_new);
#endif
}

// main function body
int main(int argc, char **argv)
{
    init_arguments(argc, argv);
    init_pop();
    write_parameters();
    write_data_headers();

    // there are this number of n probabilities:
    // 1. mortality hawk on a 1 hawk 1 dove patch
    // 2. mortality hawk on a 2 hawk patch
    // 3. mortality dove on a 1 hawk 1 dove patch
    // 4. mortality dove on a 2 dove patch
    //
    size_t nprobs = 4;

    double probs[nprobs]; // vector with probabilities

    double prob; // actual event probability sampled from 
                // the cumulative prob distribution

    // run the markov process
    for (generation = 0; generation < NumGen; ++generation)
    {
        // generate cumulative prob. distribution of possible
        // events
        double sum_probs = 0;

        /// 1. mortality maladapted on a mm, e1 patch
        sum_probs += 2 * mum * Npatches[0][0];
        probs[0] = sum_probs;

        // 2. mortality maladapted on a mm, e2 patch
        sum_probs += (1.0-(.5*v - .5*c)) * Npatches[2];
        probs[1] = sum_probs;

        // 3. mortality dove on a 1 hawk 1 dove patch
        sum_probs += .5 * (1.0-0) * Npatches[1];
        probs[2] = sum_probs;
    
        // 4. mortality dove on a 2 dove patch
        sum_probs += (1.0-.5*v) * Npatches[0];
        probs[3] = sum_probs;

        // decide which event will happen
        //
        // draw value from cumul distribution
        prob = gsl_rng_uniform(r) *  sum_probs;

        // 1. mortality hawk on a 1 hawk 1 dove patch
        //
        if (prob <= probs[0])
        {
            mortality(1, true);
        } else if (prob <= probs[1])
        {
            mortality(2, true);
        }
        else if (prob <= probs[2])
        {
            mortality(1, false);
        }
        else if (prob <= probs[3])
        {
            mortality(0, false);
        }
    
        if (generation % skip == 0)
        {
            write_data();
        }
    } // end for size_t _generation

    write_data();

    exit(1);
}
