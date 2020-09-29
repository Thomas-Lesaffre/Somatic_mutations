#include "header.h"
#include "mt.h"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <sstream>
//#include <gmp.h>
//#include <mpfr.h>

using namespace std;

// Random number generator:

MTRand rnd;

// Pointers on input and output files:

FILE * fichierE;

int main()
{

    cout << "Program initialization" << "\n";

	// Parameters:

	int N, Ng, step_mes, it, iteration, trigger_u, iter;
	double a, sig, s, h, c, t, E, L, init_u, u, step_u, gt, gga, ggb, sg;

	result Rslt;
	
	system("./sauvegarde.sh");

	// Opens input and output files:

	bool end;
	ouvrirFichierE();
	end = false;
	
	int no = 1;

	do
	{
        //reads parameter values from input file:

		end = lireFichier(N, Ng, step_mes, a, sig, s, h, c, t, E, L, init_u, u, step_u, trigger_u, gt, gga, ggb, sg, it);
					   	 // end = true if end of input file
        if(!end) cout << "\n___________________________________\n" << "\n";

        if(!end)
            for (iter = 0; iter < it; iter++)
            {

            iteration = iter+1;

            cout << "\n" << "Iteration Number : "<< iter << "\n";

			// Simulation:

 			Rslt = recursion(N, Ng, step_mes, a, sig, s, h, c, t, E, L, init_u, u, step_u, trigger_u, gt, gga, ggb, sg, iteration);
 						
            }

            no++;

	}while(!end);

	// Closes files:
	
    fclose(fichierE);
		
	return 0;
}


