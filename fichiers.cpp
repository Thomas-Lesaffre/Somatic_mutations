// Functions to open input and output files,
// read parameter values and write them in output file

#include "header.h"
#include <iostream>
#include <fstream>
#include <sstream>

//#include <gmp.h>
//#include <mpfr.h>

using namespace std;

extern FILE * fichierE;

//Opens input file:

void ouvrirFichierE()
{
	fichierE = fopen(fichierLecture,"r");
	if (!fichierE)
		cout << "The file " << fichierLecture << " doesn't exist!\n\n";
}

// Reads parameter values from input file,
// returns 1 if it reaches the end of input file, else returns 0:

bool lireFichier(int &Nr, int &Ngr, int &step_mesr, double &ar, double &sigr, double &sr, double &hr, double &cr, double &tr, double &Er, double &Lr, double &init_ur, double &ur, double &step_ur, int &trigger_ur, double &g_tr, double &g_gar, double &g_gbr, double &sgr,  int &itr)
{
	int z;
	bool term;
	do {z = fgetc(fichierE);} while (!((z == '*') || (z == EOF)));
		// Lines with parameter sets must begin with *
	if (z == EOF)
	{
		cout << "\nEnd of input file\n";
		term = true;
	}
	else
	{
		fscanf(fichierE," %d",&Nr);
		fscanf(fichierE," %d",&Ngr);
		fscanf(fichierE," %d",&step_mesr);
		fscanf(fichierE," %lf",&ar);
		fscanf(fichierE," %lf",&sigr);
		fscanf(fichierE," %lf",&sr);
		fscanf(fichierE," %lf",&hr);
		fscanf(fichierE," %lf",&cr);
		fscanf(fichierE," %lf",&tr);
		fscanf(fichierE," %lf",&Er);
		fscanf(fichierE," %lf",&Lr);
		fscanf(fichierE," %lf",&init_ur);	
		fscanf(fichierE," %lf",&ur);
		fscanf(fichierE," %lf",&step_ur);
		fscanf(fichierE," %d",&trigger_ur);							
		fscanf(fichierE," %lf",&g_tr);
		fscanf(fichierE," %lf",&g_gar);
		fscanf(fichierE," %lf",&g_gbr);
		fscanf(fichierE," %lf",&sgr);
		fscanf(fichierE," %d",&itr);							

        term = false;
	}
	
	return term;
}
