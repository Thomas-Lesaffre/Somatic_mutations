#include "header.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <vector>
#include <cmath>
#include <csignal>
#include <algorithm>
#include "mt.h"
#include <cstring>
#include <string>


//#include <gmp.h>
//#include <mpfr.h>

using namespace std;

extern MTRand rnd;

extern FILE * fichierE;
extern FILE * fichierS;

extern FILE * file_mut;
extern FILE * file_car;

// For stopping the program with Ctrl-C:

bool cntl_c_bool = false;
void cntl_c_handler (int bidon)
{
	cntl_c_bool = true;
}

result recursion(int Nv, int Ngv, int step_mesv, double av, double sigv, double sv, double hv, double cv, double tv, double Ev, double Lv, double init_uv, double uv, double step_uv, int trigger_uv, double g_tv, double g_gav, double g_gbv, double sgv, int itv)
{

    // Naming the output file :

	char nomFichier[256];
	
	stringstream nomF;
	
	nomF << "mut_rate_multi_N" << Nv << "_Ng" << Ngv << "_step" << step_mesv << "_a" << av << "_sig" << sigv << "_s" << sv << "_h" << hv << "_c" << cv << "_t" << tv << "_E" << Ev << "_L" << Lv << "_initU" << init_uv << "_umod" << uv << "_du" << step_uv << "_trig" << trigger_uv << "_gt" << g_tv << "_gga" << g_gav << "_ggb" << g_gbv << "_sg" << sgv << "_it" << itv << ".txt";
	nomF >> nomFichier;
		
	ofstream fout(nomFichier);
	
	fout << "mc u del01 del02 del12" << endl;
		
   	// For stopping prog. with Ctrl-C :
   	
	int stp;
	
   	cntl_c_bool = false;
	signal(SIGINT, cntl_c_handler);
	
	// Declaring stuff
	
	int twoN = 2*Nv; 
	double twoL = 2*Lv; 
	
	int i, j, k, nb, g, g_s, m, k1, k2, n_par1, n_par2, n_sec1, n_sec2, k_first, k_last;
	double S, max_fec, max_fec_within, mc, n_mc, n0, n1, n2, r0, r1, r2, wHet, wHom, U, del_01, del_02, del_12, prev_fec, Nmut, mes_U, new_u1, new_u2, age_ind, age_mut, age_sec, tmc;
	
	mut tmp_mut, tmp_mut2;
	
	vector<double> v_fec, v_u, tmp_sel1, tmp_sel2, fec_within, sel_pos_par1_1, sel_pos_par1_2, sel_pos_par2_1, sel_pos_par2_2;
	vector<int> vec_death, v_ind, v_sec;
	
	chr * pop = new chr[twoN];
	chr * par = new chr[twoN];
	
	chr off1, off2;
	
	result Res;
			
	// Fitness effects 
	
	wHet = 1 - sv*hv;
	wHom = 1 - sv;
		
	// Initialization
	
	S = (Ev - 1)/Ev; // Retrieve survival probability from life expectancy specified in the parameters.
	
	if(sgv == 0)
	{
	
		for(i=0; i<Nv; i++)
		{
			nb = 2*i;
		
			// Set reproduction states (S=0,1,2) stochastically according to parameters.
		
			if(rnd.rand() < av)
			{
				if(rnd.rand() < sigv)
				{
					pop[nb].S = 1;
					pop[nb+1].S = 1;
				}
				else
				{
					pop[nb].S = 2;
					pop[nb+1].S = 2;			
				}
			}
			else
			{
				pop[nb].S = 0;
				pop[nb+1].S = 0;
			}
		
			// Set individuals' mutation rate to initial value.
		
			pop[nb].u = init_uv;
			pop[nb+1].u = init_uv;
		
			// Set reproduction counter to zero.
		
			pop[nb].rep = 0;
			pop[nb+1].rep = 0;
		
			// Set individual age (that is, rank of the most recent section), to one.
		
			pop[nb].ind_age = 1;
			pop[nb+1].ind_age = 1;
		}
	
		// Set number of mutations per genome and number of measured individuals to zero.
	
		mc=0;
		n_mc=0;
	}
	else
	{
	
		// Reading back up files
		
		char fm[256];
		char fc[256];
		
		nomF.str("");
		nomF.clear();
		
		nomF << "SAVE_MUT_" << sgv << endl;
		nomF >> fm;
		nomF.str("");
		nomF.clear();
		
		nomF << "SAVE_CAR_" << sgv << endl;
		nomF >> fc;
		nomF.str("");
		nomF.clear();
		
		read_pop(fm,fc,pop);
		
		cout << "reading completed." << endl;
		
		// Set number of mutations per genome and number of measured individuals to zero.
	
		mc=0;
		n_mc=0;		
	}
	
	// Begin recursion, iterating over generations.
	
	for(g=sgv; g<Ngv; g++)
	{	
		if (cntl_c_bool) // Control for stopping: if Ctrl+C is pressed, the program will break the loop next time it reachs this point.
		{
			stp = 1;
			break;
		}
		
		// cout << "Gen: " << g << endl;
		
		// Build the fecundity table, which has three columns. For each line, we have:
		// Cumulated fecundity | the individual number (from 0 to Nv - 1) | the section number within individual (from 1 to .ind_age).
		
		//First, make sure the vectors are clear.
		
		v_fec.clear();
		v_ind.clear();
		v_sec.clear();
		
		for(i=0; i<Nv; i++) // For each individual,
		{						
			nb = 2*i;
			
			// Copy individual into the 'par' vector.
			
			par[nb] = pop[nb];
			par[nb+1] = pop[nb+1];
			
			U = (par[nb].u + par[nb+1].u)/2; // Individual [i] mutation rate.
			
			for(j=1; j<=par[nb].ind_age; j++) // For each section (from 1 to .ind_age), we create a subset vector of the .sel vectors, with the mutations present in said section, and calculate fitness.
			{
				// Clear temporary vectors.
				
				tmp_sel1.clear();
				tmp_sel2.clear();
				
				// Retrieve section genotype on the first chromosome.
				
				if(par[nb].sel.size() > 0) // If there is at least one mutation (this 'if' is there to prevent segmentation fault).
				{
					for(k=0; k<par[nb].sel.size(); k++) // iterate over the mutations present.
					{
						if(par[nb].sel[k].age <= j) // if the k^th mutation is of .age <= to the considered age 'j', copy its position in the temporary vector.
						{
							tmp_sel1.push_back(par[nb].sel[k].pos);
						}
					}
				}
				
				// Do the same for the second chromosome.
				
				if(par[nb+1].sel.size() > 0) 
				{
					for(k=0; k<par[nb+1].sel.size(); k++)
					{
						if(par[nb+1].sel[k].age <= j)
						{
							tmp_sel2.push_back(par[nb+1].sel[k].pos);
						}
					}
				}
				
				// Using the chromosomes and the individual's mutation rate, calculate fitness,
				
				v_fec.push_back(cost(gam(U, g_gav, g_gbv, g_tv),cv,tv)*fitness_vec(tmp_sel1, tmp_sel2, wHet, wHom)); // push_back in the v_fec vector.
					// The 'cost()' function calculates the cost of replication fidelity.
					// The 'gam()' function calculates the meiotic mutation rate using the somatic mutation rate.
					// The 'fitness_vec()' function calculates the fitness effect of deleterious mutations using the temporary vectors constructed above.
				v_ind.push_back(i); // push_back individual number in the v_ind vector.
				v_sec.push_back(j); // push_back section number in the v_sec vector.
			}
		}
		
		max_fec = *std::max_element(v_fec.begin(), v_fec.end()); // Retrieve the maximal fitness from the v_fec vector.
		
		vec_death.clear(); // Clear the vector indicating death.
		
		for(i=0; i<Nv; i++) // For each individual,
		{
			nb = 2*i;
			
			if(rnd.rand() < S) // If the individual survives.
			{
			
				U = (pop[nb].u + pop[nb+1].u)/2; // Determine its mutation rate.
				
				vec_death.push_back(0); // Indicate that it survived.
				
				// Increase individual age by 1.
				
				pop[nb].ind_age += 1;
				pop[nb+1].ind_age += 1;
				
				tmp_mut.age = pop[nb].ind_age;

				// Somatic mutations occur and are sorted.
								
				Nmut=poisdev(U);
				for(m=0; m<Nmut; m++)
				{
					tmp_mut.pos = twoL*rnd.rand();
					pop[nb].sel.push_back(tmp_mut);
				}
				if(Nmut > 0)
				{
					sort(pop[nb].sel.begin(),pop[nb].sel.end(), [](const mut & mut1, const mut & mut2) { return ( mut1.pos < mut2.pos ); });	
				}
							
				Nmut=poisdev(U);
				for(m=0; m<Nmut; m++)
				{
					tmp_mut.pos = twoL*rnd.rand();
					pop[nb+1].sel.push_back(tmp_mut);

				}

				if(Nmut > 0)
				{
					sort(pop[nb+1].sel.begin(),pop[nb+1].sel.end(), [](const mut & mut1, const mut & mut2) { return ( mut1.pos < mut2.pos ); });
				}
			}
			else
			{
			
				// If the individual dies.
				
				vec_death.push_back(1); // Indicate that it died.
				
				// Sample a mother : select a position in the v_fec vector.
				
				do{
					k1=rnd.randInt(v_fec.size() - 1);
				}while(v_fec[k1]/max_fec < rnd.rand());	
				
				// Determine how the offspring is produced.
				
				if(rnd.rand() < av)
				{
													
					if(rnd.rand() < sigv)
					{
						// If selfing occurs within the same section, the father is the mother.
						
						k2=k1;	
						
						off1.S=1;
						off2.S=1;	
					}
					else
					{
												
						// If selfing occurs between different sections within the individual, we need to compute a new fitness vector, to select a father.
						
						fec_within.clear(); // Make sure the vector starts empty.
						
						age_ind = par[2*v_ind[k1]].ind_age; // Retrieve individual age.
						
						k_first = k1 - v_sec[k1] + 1; // Position of the first section within the individual in the general fitness vector.
						k_last = k1 + age_ind - v_sec[k1]; // Position of the last section within the individual in the general fitness vector.

						if(k_first != k_last)
						{						
							fec_within.assign(v_fec.begin() + k_first, v_fec.begin() + k_last); // Extract the subvector of fecundities,
						
							max_fec_within = *std::max_element(fec_within.begin(), fec_within.end());
						
							do{
								k2=rnd.randInt(fec_within.size() - 1);
							}while(fec_within[k2]/max_fec_within < rnd.rand()); // here, we select a section within the individual.
						
							k2 += k_first; // That way, we retrieve the position of the selected section in the general fitness table.
						}
						else
						{
							k2 = k1;
						}
						
						// Set reproduction status accordingly.

						off1.S = 2;
						off2.S = 2;
					}
				}
				else
				{
					// If outcrossing occurs, select a father.
					
					do{
						k2=rnd.randInt(v_fec.size() - 1);
					}while(v_fec[k2]/max_fec < rnd.rand());	
					
					off1.S = 0;
					off2.S = 0;
				}
				
				// Now that we have selected the parents, we need to reconstruct the chromosomes of the sections.
				
				n_par1=v_ind[k1]; // Individual number of parent 1.
				n_sec1=v_sec[k1]; // Section number.

				n_par2=v_ind[k2]; // Individual number of parent 2.
				n_sec2=v_sec[k2]; // Section number.
				
				// Parent1, chromosome 1.
				
				sel_pos_par1_1.clear(); // Make sur the temporary mutations vector is clear.
				
				if(par[2*n_par1].sel.size() != 0) // First chromosome of the mother bears at the least one mutation,
				{
					// Sort mutations with respect to age
					sort(par[2*n_par1].sel.begin(), par[2*n_par1].sel.end(), [](const mut & mut1, const mut & mut2) { return ( mut1.age < mut2.age ); });
					
					age_mut=par[2*n_par1].sel.back().age; // Retrieve the rank of the most recent mutation.

					if(n_sec1 < age_mut) // If the age of the section is < the rank of the most recent mutation,
					{
						k=-1;
						do{
							k++; // Increase k by 1.
							
							if(par[2*n_par1].sel[k].age <= n_sec1) // If the rank of the mutation is < rank of the section, copy to temp. vector.
							{
								sel_pos_par1_1.push_back(par[2*n_par1].sel[k].pos);
							}
							
						}while(par[2*n_par1].sel[k].age <= n_sec1); // Copy to position vector up to the section.			
					}
					else
					{
						// If n_sec1 >= age_mut, copy all the mutation in the temp. vector.
						for(k=0; k<par[2*n_par1].sel.size(); k++)
						{
							sel_pos_par1_1.push_back(par[2*n_par1].sel[k].pos);
						}
					}
					
					sort(sel_pos_par1_1.begin(), sel_pos_par1_1.end()); // Sort the new vector of positions, to prepare it for recombination.					
					sort(par[2*n_par1].sel.begin(), par[2*n_par1].sel.end(), [](const mut & mut1, const mut & mut2) { return ( mut1.pos < mut2.pos ); }); // Sort back to 'pos' sorting.
				}
				
				// Parent1, chromosome 2.
				
				sel_pos_par1_2.clear();
				if(par[2*n_par1 + 1].sel.size() != 0)
				{	
					sort(par[2*n_par1+1].sel.begin(), par[2*n_par1+1].sel.end(), [](const mut & mut1, const mut & mut2) { return ( mut1.age < mut2.age ); });
					age_mut=par[2*n_par1+1].sel.back().age;

					if(n_sec1 < age_mut)
					{
						k=-1;
						do{
							k++;
							if(par[2*n_par1+1].sel[k].age <= n_sec1)
							{
								sel_pos_par1_2.push_back(par[2*n_par1+1].sel[k].pos);
							}
						}while(par[2*n_par1+1].sel[k].age <= n_sec1); // Copy to position vector up to the section.			
					}
					else
					{
						for(k=0; k<par[2*n_par1+1].sel.size(); k++)
						{
							sel_pos_par1_2.push_back(par[2*n_par1+1].sel[k].pos);
						}
					}
					
					sort(sel_pos_par1_2.begin(), sel_pos_par1_2.end()); // Sort the new vector of positions, to prepare it for recombination.					
					sort(par[2*n_par1+1].sel.begin(), par[2*n_par1+1].sel.end(), [](const mut & mut1, const mut & mut2) { return ( mut1.pos < mut2.pos ); }); // Sort back to 'pos' sorting
				}
				
				
				// Parent2, chromosome 1.

				sel_pos_par2_1.clear();
				if(par[2*n_par2].sel.size() != 0)
				{			
					sort(par[2*n_par2].sel.begin(), par[2*n_par2].sel.end(), [](const mut & mut1, const mut & mut2) { return ( mut1.age < mut2.age ); });
					age_mut=par[2*n_par2].sel.back().age;

					if(n_sec2 < age_mut)
					{
						k=-1;
						do{
							k++;
							if(par[2*n_par2].sel[k].age <= n_sec2)
							{
								sel_pos_par2_1.push_back(par[2*n_par2].sel[k].pos);
							}
						}while(par[2*n_par2].sel[k].age <= n_sec2); // Copy to position vector up to the section.
					}
					else
					{
						for(k=0; k<par[2*n_par2].sel.size(); k++)
						{
							sel_pos_par2_1.push_back(par[2*n_par2].sel[k].pos);
						}
					}
					
					sort(sel_pos_par2_1.begin(), sel_pos_par2_1.end()); // Sort the new vector of positions, to prepare it for recombination.
					sort(par[2*n_par2].sel.begin(), par[2*n_par2].sel.end(), [](const mut & mut1, const mut & mut2) { return ( mut1.pos < mut2.pos ); }); // Sort back to 'pos' sorting.			
				}
				
				// Parent2, chromosome 2.
				
				sel_pos_par2_2.clear();
				if(par[2*n_par2+1].sel.size() != 0)
				{
					sort(par[2*n_par2 + 1].sel.begin(), par[2*n_par2 + 1].sel.end(), [](const mut & mut1, const mut & mut2) { return ( mut1.age < mut2.age ); });
					age_mut=par[2*n_par2 + 1].sel.back().age;

					if(n_sec2 < age_mut)
					{
						k=-1;
						do{
							k++;
							if(par[2*n_par2+1].sel[k].age <= n_sec2)
							{
								sel_pos_par2_2.push_back(par[2*n_par2+1].sel[k].pos);
							}
						}while(par[2*n_par2+1].sel[k].age <= n_sec2); // Copy to position vector up to the section.
					}
					else
					{
						for(k=0; k<par[2*n_par2+1].sel.size(); k++)
						{
							sel_pos_par2_2.push_back(par[2*n_par2+1].sel[k].pos);
						}
					}
					
					sort(sel_pos_par2_2.begin(), sel_pos_par2_2.end()); // Sort the new vector of positions, to prepare it for recombination
					sort(par[2*n_par2+1].sel.begin(), par[2*n_par2+1].sel.end(), [](const mut & mut1, const mut & mut2) { return ( mut1.pos < mut2.pos ); }); // Sort back to 'pos' sorting.
				}
				
				// Mutation in Parent 1
				
				// Mutation at the modifier.
				
				new_u1 = par[2*n_par1].u;
				new_u2 = par[2*n_par1 + 1].u;
				
				U = (new_u1 + new_u2)/2; // Retrieve mutation rate.
				
				if(g > trigger_uv) // If mutation at the modifier is triggered,
				{
					if(rnd.rand() < gam(U, g_gav, g_gbv, g_tv)*uv)
					{
						new_u1 = new_u1 + step_uv*gaussdev();

						if(new_u1 < 0)
						{
							new_u1 = 0;
						}
					}
					
					if(rnd.rand() < gam(U, g_gav, g_gbv, g_tv)*uv)
					{
						new_u2 = new_u2 + step_uv*gaussdev();

						if(new_u2 < 0)
						{
							new_u2 = 0;
						}
					}
				}
				
				// Mutations at selected loci.
				
				U = (new_u1 + new_u2)/2; // Retrieve mutation rate.
						
				Nmut=int(poisdev(gam(U, g_gav, g_gbv, g_tv))); // Sample the number of mutations.	
				for(m=0; m<Nmut; m++)
				{
					sel_pos_par1_1.push_back(twoL*rnd.rand());
				}

				if(Nmut > 0)
				{
					sort(sel_pos_par1_1.begin(), sel_pos_par1_1.end());				
				}
								
				Nmut=int(poisdev(gam(U, g_gav, g_gbv, g_tv)));		
				for(m=0; m<Nmut; m++)
				{
					sel_pos_par1_2.push_back(twoL*rnd.rand());
				}

				if(Nmut > 0)
				{
					sort(sel_pos_par1_2.begin(), sel_pos_par1_2.end());
				}
				
				// Recombination
				
				off1.sel.clear();
				
				// Off1
				if(rnd.rand() < 0.5) // Each chromosome has a 50% chance of being transmitted.
				{
					rec_vec(off1,sel_pos_par1_1,sel_pos_par1_2,	Lv);
					off1.u = new_u1;	
				}
				else
				{
					rec_vec(off1,sel_pos_par1_2,sel_pos_par1_1,	Lv);
					off1.u = new_u2;								
				}
				
				// Mutation in Parent 2
				
				// Mutation at the modifier.
				
				new_u1 = par[2*n_par2].u;
				new_u2 = par[2*n_par2 + 1].u;
				
				U = (new_u1 + new_u2)/2; // Retrieve mutation rate.
				
				if(g > trigger_uv) // If mutation at the modifier is triggered,
				{
					if(rnd.rand() < gam(U, g_gav, g_gbv, g_tv)*uv)
					{
						new_u1 = new_u1 + step_uv*gaussdev();

						if(new_u1 < 0)
						{
							new_u1 = 0;
						}
					}
					
					if(rnd.rand() < gam(U, g_gav, g_gbv, g_tv)*uv)
					{
						new_u2 = new_u2 + step_uv*gaussdev();

						if(new_u2 < 0)
						{
							new_u2 = 0;
						}
					}
				}
				
				// Mutations at selected loci.
				
				U = (new_u1 + new_u2)/2;
						
				Nmut=int(poisdev(gam(U, g_gav, g_gbv, g_tv)));		
				for(m=0; m<Nmut; m++)
				{
					sel_pos_par2_1.push_back(twoL*rnd.rand());
				}
				
				if(Nmut > 0)
				{
					sort(sel_pos_par2_1.begin(), sel_pos_par2_1.end());	
				}
				
				Nmut=int(poisdev(gam(U, g_gav, g_gbv, g_tv)));		
				for(m=0; m<Nmut; m++)
				{
					sel_pos_par2_2.push_back(twoL*rnd.rand());
				}		
				
				if(Nmut > 0)
				{
					sort(sel_pos_par2_2.begin(), sel_pos_par2_2.end());
				}
				
				// Recombination
				
				off2.sel.clear();	
				
				// Off2
				if(rnd.rand() < 0.5) // Each chromosome has a 50% chance of being transmitted.
				{
					rec_vec(off2,sel_pos_par2_1,sel_pos_par2_2,	Lv);
					off2.u = new_u1;	
				}
				else
				{
					rec_vec(off2,sel_pos_par2_2,sel_pos_par2_1,	Lv);
					off2.u = new_u2;								
				}

				// The offspring is ready to be incorporated into the population.
				
				off1.ind_age = 1;
				off2.ind_age = 1;
				
				pop[nb] = off1;
				pop[nb+1] = off2;
				
				mc += (pop[nb].sel.size() + pop[nb+1].sel.size())/2;
				n_mc += 1;
				
				// Somatic mutation rate to reach maturity
				
				U = (pop[nb].u + pop[nb+1].u)/2; // Determine its mutation rate.
				
				tmp_mut.age = pop[nb].ind_age;

				// Somatic mutations occur and are sorted.
								
				Nmut=poisdev(U);
				for(m=0; m<Nmut; m++)
				{
					tmp_mut.pos = twoL*rnd.rand();
					pop[nb].sel.push_back(tmp_mut);
				}

				if(Nmut > 0)
				{
					sort(pop[nb].sel.begin(),pop[nb].sel.end(), [](const mut & mut1, const mut & mut2) { return ( mut1.pos < mut2.pos ); });
				}
								
				Nmut=poisdev(U);
				for(m=0; m<Nmut; m++)
				{
					tmp_mut.pos = twoL*rnd.rand();
					pop[nb+1].sel.push_back(tmp_mut);

				}

				if(Nmut > 0)
				{
					sort(pop[nb+1].sel.begin(),pop[nb+1].sel.end(), [](const mut & mut1, const mut & mut2) { return ( mut1.pos < mut2.pos ); });
				}
				
				// Increase reproduction counter for parents.
				
				par[2*n_par1].rep +=1;
				par[2*n_par1+1].rep +=1;				
				
				if(n_par1 != n_par2)
				{
					par[2*n_par2].rep +=1;
					par[2*n_par2+1].rep +=1;
				}
			} // End of 'if' over survival.
		} // End of loop over 'pop'.
		
		// At each timestep we need to check who died, and transfer .rep
		
		for(i=0; i<Nv; i++)
		{
			nb=2*i;
			
			if(vec_death[i] == 0) // If the individual survived, we transfer the .rep counter
			{
				pop[nb].rep = par[nb].rep;
				pop[nb+1].rep = par[nb+1].rep;
			}
			else
			{
				if(par[nb].S == 0)
				{
					n0 += 1;
					r0 += (par[nb].rep + par[nb+1].rep)/2;
				}
				else if(par[nb].S == 1)
				{
					n1 += 1;
					r1 += (par[nb].rep + par[nb+1].rep)/2;
				}
				else if(par[nb].S == 2)
				{
					n2 += 1;
					r2 += (par[nb].rep + par[nb+1].rep)/2;
				}
			}
			
			mes_U += (pop[nb].u + pop[nb+1].u)/2;
		}
		
		mes_U /= Nv;
		
		cout << mes_U << endl;
		
		v_u.push_back(mes_U);
		mes_U = 0;
		
		if((g % step_mesv == 0)&&(g > 0))
		{
			for(i=0; i<v_u.size(); i++)
			{
				mes_U += v_u[i];
			}
			mes_U/=v_u.size();
			v_u.clear();
			
			del_01 = 1 - (r1/n1)/(r0/n0);
			del_02 = 1 - (r2/n2)/(r0/n0);
			del_12 = 1 - (r1/n1)/(r2/n2);
			
			mc /= n_mc;
			
			if(sigv == 1)
			{
				fout << mc << " "  << mes_U << " " << del_01 << " " << 0 << " " << 0 << endl;
			}
			else
			{
				fout << mc << " "  << mes_U << " " << del_01 << " " << del_02 << " " << del_12 << endl;			
			}
			
			r0=0;
			n0=0;
			
			r1=0;
			n1=0;
			
			r2=0;
			n2=0;
			
			mc=0;
			n_mc=0;
			
			mes_U=0;
			
			// For now, we save every measuring step.
			
			save_pop(pop, g, 2*Nv);
			g_s=g;
		}
		
		if(g == Ngv - 1) // If it is the last time we run, rename save file !
		{
			char fm[256];
			char fm_n[256];
			
			char fc[256];
			char fc_n[256];
		
			nomF.str("");
			nomF.clear();
		
			nomF << "SAVE_MUT_" << g_s << endl;
			nomF >> fm;
			nomF.str("");
			nomF.clear();

			nomF << itv << "_E" << Ev << "_SAVE_MUT_" << g_s << endl;
			nomF >> fm_n;
			nomF.str("");
			nomF.clear();
			
			nomF << "SAVE_CAR_" << g_s << endl;
			nomF >> fc;
			nomF.str("");
			nomF.clear();
			
			nomF << itv << "_E" << Ev << "_SAVE_CAR_" << g_s << endl;
			nomF >> fc_n;
			nomF.str("");
			nomF.clear();
			
			rename(fm,fm_n);
			rename(fc,fc_n);
		}
	} // End of loop over time.

    delete[] pop;
	delete[] par;
		
    return(Res);
}

