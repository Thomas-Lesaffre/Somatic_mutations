#ifndef DEPRESSION_H
#define DEPRESSION_H

#include <vector>
#include <iostream>
// #include "MersenneTwister.h"
#include "mt.h"
//#include <gmp.h>
//#include <mpfr.h>

using namespace std;

// global variables

#define fichierLecture "parameters.txt" 

// definition of structure "chr" representing a chromosome:
// "M" is the Modifier locus
// "sel" is a vector containing the positions of deleterious alleles along the chromosome

struct mut
{
	double pos;
	double age;
};

struct chr
{
      vector<mut> sel; // Selected loci
      double u; // Mutation rate modifier
      
      double S;
      double rep;
      double ind_age;

      chr(){}; // Default constructor
      ~chr(){}; // Default destructor
};

struct result
{
       double pouet;
};



// Prototypes of functions

void save_pop(chr P[], int G, int N);
void read_pop(char * file_name_mut, char * file_name_car, chr * P);

void ouvrirFichierE();

void ouvrirFichierS();

bool lireFichier(int &Nr, int &Ngr, int &step_mesr, double &ar, double &sigr, double &sr, double &hr, double &cr, double &tr, double &Er, double &Lr, double &init_ur, double &ur, double &step_ur, int &trigger_ur, double &g_tr, double &g_gar, double &g_gbr, double &sgr, int &itr);

result recursion(int Nv, int Ngv, int step_mesv, double av, double sigv, double sv, double hv, double cv, double tv, double Ev, double Lv, double init_uv, double uv, double step_uv, int trigger_uv,  double g_tv, double g_gav, double g_gbv, double sgv, int itv);

double cost(double U, double c, double t);
double gam(double u, double ga, double gb, double t);

double gammln(const double xx);
double poisdev(const double xm);
double binldev(const double pp, const int n);
double gaussdev();

double fitness(chr &c1, chr &c2, double wHe, double wHo);

void rec_vec(chr &res, vector<double> &v1, vector<double> &v2, double Sz);

double fitness_vec(vector<double> &c1, vector<double> &c2, double wHe, double wHo);

void cntl_c_handler(int bidon);


#endif
