#include<armadillo>
#include<iostream>
#include<vector>
#include<cstdlib>
#include<math.h>
#include<complex>

using namespace arma;
using namespace std;

#define PI 3.141592653589793238462643





//First, we introduce the desired settings of the system that we want to study

const int q = 3;  //Number of values each spin can take

const int L = 16;  //System size

const double k_B = 1.380649 * pow(10, -23);  //Boltzmann constant

const double T = 0.5; //Temperature in units of J		NEEDS TO BE REPLACE FOR T/J=0.5

const int N = L; //Total number of spins

const double pconnect = 1 - exp(1/T);  //Connection probability

const int NCLUSTERS = 1;  //Number of cluster builds in one MC step

const int NESTEPS = 10000;  //Number of equilibrium MC steps

const int NMSTEPS = 10000;  //Number of measurement MC steps

const int NBINS = 10;  //Number of measurement bins

vector <int> S(N);  //The spin array

vector <int> M(q);  //Number of spins in each state

vector <complex<double> > W(q);  //Order parameter weights





//Now we will define some methods in order to manipulate the lattice

enum dirs { RIGHT, LEFT };

int indx(int x) {

	return x;  //Make an indx on every site

}


int xpos(int i) {

	return i % L;

}


int Nbr(int i, int dir) {

	int x = xpos(i);

	switch (dir) {

	case RIGHT: return indx((x + 1) % L);
	case LEFT: return indx((x - 1 + L) % L);

	}

	return 0;

}


void FlipandBuildFrom(int s) {

	int oldstate(S[s]), newstate((S[s] + 1) % q);

	S[s] = newstate;  //Flip spin

	M[oldstate]--; M[newstate]++;  //Update spin counts

	for (int dir = 0; dir < 2; dir++) {

		int j = Nbr(s, dir);

		if (S[j] == oldstate) {

			if (rand() / (RAND_MAX + 1.) < pconnect) {

				FlipandBuildFrom(j);

			}

		}

	}

}
