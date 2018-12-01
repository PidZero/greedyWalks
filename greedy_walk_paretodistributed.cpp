/* (c) Johannes Neidhart, 2012 
 * This little piece of code simulates a greedy adaptive walk in a pareto distributed 
 * Rough Mount Fuki Fitness landscape. For details see for example
 * http://kups.ub.uni-koeln.de/5878/
 */

// load the necessary libraries

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

// GSL will be used to generate random numbers. The library needs
// some variables, which are better defined globally. It would not be good
// start from the beginning, everytime a subroutine is called,
// that's why I chose them to be global.
long seed;	
const gsl_rng_type * T;
gsl_rng * ran;

// dist_NH( vector <double> &F, int d, double c, double a)
// vector <double> &F is used to return the new fitness values
// int d is the new distance to the reference sequence
// double c is the smoothness of the landscape
// double a is the parameter of the Pareto distribution
//
// This function calculates the new neighbourhood after a step is taken
void dist_NH( vector <double> &F, int d, double c, double a){
	for(int i=0; i< d;i++){
		F[i] = -c * (double)(d-1) + gsl_ran_pareto(ran, a, 1);
	}
	
	for(int j=d; j< (int)F.size();j++){
		F[j] = -c * (double)(d+1) + gsl_ran_pareto(ran, a, 1);
	}
}


// double smax(vector <double> F, int &scnd_largest)
// vector <double> F is the momentary fitness landscape
// int &scnd_largest gives back the second largest value in the neighbourhood
// output gives back the largest value in the neighbourhood
//
// This function calculates the largest fitness value in the
// current neighbourhood
double smax(vector <double> F, int &scnd_largest){
	double h=F[0];
	for(int i=1; i<F.size(); i++){
		if(h<F[i]){
			h=F[i];
			scnd_largest=i;
		}
	}
	return(h);
}


// int walk(int L, double c, double a)
// int L is the dimension of the hypercube
// double c is the smoothness of the landscape
// double a is the parameter of the Pareto distribution
// output is the greedy walk length of one walk in such a landscape
//
// This function calculates simulates one walk in a Rough Mount Fuji Landscape
// and returns the number of steps, until a fitness maximum is reached
int walk(int L, double c, double a){
	vector <double> F;
	vector <double> F_sort;
	int d=L;
	double F_a;
	bool nomax;
	bool step;
	int l=0;
	int h;
	double F_h;
	double sm;
	l=0;
	
	F.resize(L);
	dist_NH(F,d,c,a);
	sm = smax(F, h);
	F_a = -c * (L) + gsl_ran_pareto (ran,a,  1);
	nomax=true;
	step=true;
	while(nomax == true){	
		if(step == true){ 				
			if(sm>F_a){
				nomax = true;			
				step=false;
			}else{
				nomax=false;
			}
		}
		
		F_h = sm;					
		if(h<d){
			d--;
		}else{
			d++;
		}				
		dist_NH(F,d,c, a);			
						
		F_a = F_h;			
		l++;				
		step=true;
		sm = smax(F, h);
	}
	
	return(l);		
}
	
	
// This program calcluates the mean greedy adaptive walk length in a 
// Rough Mount Fuji landscape.
int main(){
	// Setup of the RNG
	T = gsl_rng_default;			
	ran = gsl_rng_alloc (T);
	seed = time (NULL) * getpid();
	gsl_rng_set (ran, seed);
	
	int h;		// a helper variable
	int p=1000; // the statistics, how many runs should be simulated
	int L=1000; // the dimension of the hypercube
	double c;	// Sets how smooth the landscape is. The larger, the smoother.
	ofstream f; // Used to write results to file
	double a=5; // Pareto parameter
	f.open("output.dat");
	// This loop does a sweep through different rough landscapes and simulates
	// p walks throughout to calculate the mean walk length. The result 
	// is then written to f.
	for(c=0.01; c <=1; c+=0.01){
		h=0;
		for(int i=0; i<p; i++){
			h+=walk(L,c, a);
			
		}
		f<<c<<"\t"<<(double)h/(double)p<<endl;
	}
	

	f.close();
}
