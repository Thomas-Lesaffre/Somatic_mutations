// Functions for survival and fitness computation and for recombination:

#include "header.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include "mt.h"
#include <fstream>
#include <sstream>

//#include <gmp.h>
//#include <mpfr.h>

using namespace std;

extern MTRand rnd;

// Calculate a fitness using vectors and not chrs.

double fitness_vec(vector<double> &c1, vector<double> &c2, double wHe, double wHo)
{

	double w = 1.0;
	
	int s1 = c1.size() - 1;
	int s2 = c2.size() - 1;

	while ((s1 > -1) || (s2 > -1)) // As long as there is at least one deleterious mutation,
	{
		if (s1 == -1) // If there is no del. mut. on chr1, then all mutations on chr2 are in heterozygous state
		{
			w *= pow(wHe, s2 + 1); //  Thus, fitness is a power function : w = w*(1-hs)^(s2_a+1) where s2_a+1 is nb. of mut. on chr2.

			break;
		}

		if (s2 == -1) // If there is no del. mut. on chr2, same process as before but with the other chr.
		{
			w *= pow(wHe, s1 + 1);
			break;
		}

		// heterozygous mutation on c1:

		if (c1[s1] > c2[s2]) // If the last mutation of chr1 is further along the chr. than that of chr2,
		{
			w *= wHe; // Multiply fitness by (1-hs), that is, add an heterozygous mutation
			s1--; // Diminish s1_a by 1, so that next time, the previous mutation in the vector is analysed.
			continue; // Skip the rest of the loop, go back to top.
		}

		// heterozygous mutation on c2:

		if (c1[s1] < c2[s2]) // Same process but with chr2 this time.
		{
			w *= wHe;
			s2--;
			continue;
		}

		// When we get here, it means we have found c1.B[s1] = c2.B[s2] (none of the 'if' were found true), that is, an homozygous mutation.

		w *= wHo; // Multiply fitness by (1-s)

		s1--; // Diminish mutation count by one on chr1
		s2--; // and on chr2.
	
	} // Back to the top.


	return w;
}

// Recombination

void rec_vec(chr &res, vector<double> &v1, vector<double> &v2, double Sz)
{
	int sz1 = v1.size();
	int sz2 = v2.size();
	
	int s1 = sz1 - 1;
	int s2 = sz2 - 1;
	
	res.sel.clear();

	res.rep = 0;
	
	int j;
	vector<double> Co;
	mut tmp;
	
	tmp.age = 1;

	double twoS = 2 * Sz;

	// number of cross-overs:

	int nbCo = int(poisdev(Sz));  // Nb. of crossing-overs sampled from a Poisson law with mean "Sz".

	// positions of cross-overs (between 0 and 2L) are put in the vector Co.sel

	for (j = 0; j < nbCo; j++)
    {
        Co.push_back(twoS * rnd.rand()); // Positions of C-O are are sampled at random in twoS.
        sort(Co.begin(), Co.end()); // Positions are sorted.
    }

	// Counting the number of cross-overs on the right of the Modifier locus:

	int cmpt = 0; // set C-O counter to 0.
	j = Co.size() - 1;
	
	while ((j > -1) && (Co[j] > Sz)) // While there is at least 1 C-O and and the considered one is located on right of the M-locus.
	{
		cmpt++; // raise number of C-Os by 1.
		j--; // lower the number of C-O yet to be considered by one.
	}

	// If the number of C-Os on the right is even:

	if (cmpt % 2 == 0)
	{

		for (j = 1; j <= cmpt; j++)
		{
			if (j % 2 == 1) // First C-O results in offspring getting part of chr1, then 2nd C-O (where 'j' is even, see below)
			{               // gives part of chr2, up until the next C-O...
				// All mutations on the right of cross-over on chromosome 1 are incorporated

				while ((s1 > -1) && (v1[s1] > Co[nbCo - j]))
				{
					tmp.pos = v1[s1];
					res.sel.push_back(tmp);
					s1--;
				}
				while ((s2 > -1) && (v2[s2] > Co[nbCo - j]))
					s2--;

			}
			else
			{
				// all mutations on the right of cross-over on chromosome 2 are incorporated

				while ((s2 > -1) && (v2[s2] > Co[nbCo - j]))
				{
					tmp.pos = v2[s2];
					res.sel.push_back(tmp);
					s2--;
				}
				while ((s1 > -1) && (v1[s1] > Co[nbCo - j]))
					s1--;
			}

		}

	}

	// If the number of cross-overs on the right is odd:

	else
	{
		// For each cross-over (starting on extreme right):

		for (j = 1; j <= cmpt; j++)
		{
			if (j % 2 == 1)
			{
				// all mutations on the right of cross-over on chromosome 2 are incorporated

				while ((s2 > -1) && (v2[s2] > Co[nbCo - j]))
				{
					tmp.pos = v2[s2];
					res.sel.push_back(tmp);
					s2--;
				}
				while ((s1 > -1) && (v1[s1] > Co[nbCo - j]))
					s1--;
			}
			else
			{
				// all mutations on the right of cross-over on chromosome 1 are incorporated

				while ((s1 > -1) && (v1[s1] > Co[nbCo - j]))
				{
					tmp.pos = v1[s1];
					res.sel.push_back(tmp);
					s1--;
				}
				while ((s2 > -1) && (v2[s2] > Co[nbCo - j]))
					s2--;
			}

		}

	}


	// mutations between the Modifier locus and the nearest cross-over on the right (on chromosome 1):

	while ((s1 > -1) && (v1[s1] > Sz))
	{
		tmp.pos = v1[s1];
		res.sel.push_back(tmp);
		s1--;
	}
	while ((s2 > -1) && (v2[s2] > Sz))
		s2--;

	// number of cross-overs on the left of the locus:

	int frst = nbCo - cmpt;

	// for each cross-over (starting with the nearest to the locus on the left):

	for (j = 1; j <= frst; j++)
	{
		if (j % 2 == 1)
		{
			// all mutations on the right of cross-over on chromosome 1 are incorporated

			while ((s1 > -1) && (v1[s1] > Co[frst - j]))
			{
				tmp.pos = v1[s1];
				res.sel.push_back(tmp);
				s1--;
			}
			while ((s2 > -1) && (v2[s2] > Co[frst - j]))
				s2--;
		}
		else
		{
			// all mutations on the right of cross-over on chromosome 2 are incorporated

			while ((s2 > -1) && (v2[s2] > Co[frst - j]))
			{
				tmp.pos=v2[s2];
				res.sel.push_back(tmp);
				s2--;
			}
			while ((s1 > -1) && (v1[s1] > Co[frst - j]))
				s1--;
		}


	}


	// mutations on the left of the left-most cross-over:

	if (frst % 2 == 0)
	{

		while (s1 > -1)
		{
			tmp.pos=v1[s1];
			res.sel.push_back(tmp);
			s1--;
		}

    }
    else
    {

		while (s2 > -1)
		{
			tmp.pos=v2[s2];
			res.sel.push_back(tmp);
			s2--;
		}

    }

	// sorts mutations on offspring chromosome:

	sort(res.sel.begin(), res.sel.end(), [](const mut & mut1, const mut & mut2) { return ( mut1.pos < mut2.pos ); }); 
}

// Fidelity cost function: various types
	// Type 1: Power function.
	// Type 2: Exponential function.

double cost(double U, double c, double t)
{
	double f;
	
	if(t == 1)
	{
		if(U == 0)
		{
			f = 0;
		}
		else
		{
			f = pow(U,c);
		}
	}
	else if(t == 2)
	{
		if(U == 0)
		{
			f = 0;
		}
		else
		{
			f = exp(-c/U);
		}
	}
	
	return f;
}

// Meiotic mut. rate: various types
	// Type 1: Linear function.
	// Type 2: Affine function.
	// Type 3: Exponential function.
	// Type 4: Power function.

double gam(double u, double ga, double gb, double t)
{
	double mu;
	
	if(t == 1)
	{
		mu = ga*u;	
	}
	else if(t == 2)
	{
		mu = ga*u + gb;		
	}
	else if(t == 3)
	{
		mu = exp(u/ga) - 1;
	}
	else if(t == 4)
	{
		mu = pow(mu,ga);
	}
	
	return mu;
}


void save_pop(chr P[], int G, int N)
{
	// Declaring elements
	char mut[256]; // Chars for names
	char car[256];
	
	int i,j;
	
	ofstream fmut, fcar; // Creating outward stream
	
	stringstream strn; // Stringstream to build names
	
	system("rm SAVE_*"); // Remove previous saves.
	
	// Creating save files
	strn << "SAVE_MUT_" << G;
	strn >> mut; 
	
	strn.str(""); // Resetting stringstream
	strn.clear();
	
	strn << "SAVE_CAR_" << G;
	strn >> car;
	
	// Filling files
	fmut.open(mut); // Open streams to files
	fcar.open(car); 
	
	for(i=0; i<N; i++)
	{
		// Mutations
		fmut << "*";
		for(j=0; j<P[i].sel.size(); j++)
		{
			if(j < P[i].sel.size() - 1)
			{
				fmut << P[i].sel[j].pos << " " << P[i].sel[j].age << " ";
			}
			else
			{
				fmut << P[i].sel[j].pos << " " << P[i].sel[j].age << endl;				
			}
		}
		
		// Characteristics
		fcar << "*" << P[i].u << " " << P[i].S << " " << P[i].rep << " " << P[i].ind_age << endl;
	}
	
	fmut.close();
	fcar.close();
}

void read_pop(char * file_name_mut, char * file_name_car, chr * P)
{

	FILE * file_mut;
	FILE * file_car;		
	
	// Reading and copying mutations
	
	cout << "Ouverture..." << endl;
	
	file_mut = fopen(file_name_mut,"r");
	
	cout << "... Ouvert !" << endl;
	
	stringstream tmp;
	int z,i,n;
	char ch;
	double val;
	mut tmp_mut;
	
	i=0;
	n=-1;
	
	cout << "Lecture..." << endl;
	
	do{		
		tmp.clear();
		z = fgetc(file_mut);
				
		if(z != EOF)
		{
			ch=z;
			
			if((ch != ' ')&&(ch != '*')&&(z != EOF))
			{
				tmp << ch;
			}
			else if((ch == ' ')&&(z != EOF))
			{
				tmp >> val;
				tmp.str("");
				
				if(i % 2 == 0)
				{
					tmp_mut.pos = val;
				}
				else
				{
					tmp_mut.age = val;
					
					P[n].sel.push_back(tmp_mut);
				}
				
				i++;
			}
			else if((ch == '*')&&(z != EOF))
			{
				if(i != 0)
				{
					tmp >> val;
					tmp.str("");

					tmp_mut.age = val;	
					
					P[n].sel.push_back(tmp_mut);
				
					i=0;			
				}
				
				n++;
				// cout << "n = " << n << endl;
			}
		}
		else
		{
			tmp >> val;
			tmp.str("");
			tmp_mut.age = val;			
			P[n].sel.push_back(tmp_mut);	
		}
		
				
	}while(z != EOF);
	
	cout << "... Lu !" << endl;
	
	// Reading and copying characteristics
	
	file_car = fopen(file_name_car,"r");
	
	n=-1;
	
	vector<double> vec;
	
	do{
		tmp.clear();
		z = fgetc(file_car);
		
		if(z != EOF)
		{
			ch=z;
			
			if((ch != ' ')&&(ch != '*')&&(z != EOF))
			{
				tmp << ch;
			}
			else if((ch == ' ')&&(z != EOF))
			{
				tmp >> val;
				tmp.str("");
	
				vec.push_back(val);
			}
			else if((ch == '*')&&(z != EOF))
			{
				if(n == -1)
				{
					n++;
				}
				else
				{
					tmp >> val;
					tmp.str("");
					
					vec.push_back(val);

					P[n].u = vec[0];
					P[n].S = vec[1];
					P[n].rep = vec[2];
					P[n].ind_age = vec[3];
					n++;
					vec.clear();
				}
			}
		}
		else
		{
			tmp >> val;
			tmp.str("");
			vec.push_back(val);
			
			P[n].u = vec[0];
			P[n].S = vec[1];
			P[n].rep = vec[2];
			P[n].ind_age = vec[3];
			
			vec.clear();
		}
				
	}while(z != EOF);
}

