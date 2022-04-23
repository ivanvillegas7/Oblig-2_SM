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

const double T = 0.25; //Temperature in units of J		NEEDS TO BE REPLACE FOR T/J=0.5

const int N = L; //Total number of spins

const double pconnect = 1 - exp(1 / T);  //Connection probability

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


cx_vec k(N);


void C(cx_vec& k) {

	ofstream ofile;
	ofile.open("C(r)_T_0.25.txt");		//NEEDS TO BE REPLACE FOR T / J = 0.5
	ofile << "r" << "	" << "Real part" << endl;
	ofile << scientific;

	for (int r = 0; r < N; r++) {

		complex<double>m_0(0., 0.);
		complex<double>m_r(0., 0.);
		complex<double>m_p(0., 0.);

		for (int t = 0; t < NMSTEPS; t++) {


			for (int c = 0; c < NCLUSTERS; c++) {

				FlipandBuildFrom(rand() % N);

			}


			m_0 = m_0 + conj(complex<double>(cos(2 * PI * S[0] / q), sin(2 * PI * S[0] / q)));
			m_r = m_r + complex<double>(cos(2 * PI * S[r] / q), sin(2 * PI * S[r] / q));
			m_p = m_p + m_0 * m_r;

		}

		complex<double>division(NMSTEPS, 0.);

		m_0 = m_0 / division;
		m_r = m_r / division;
		m_p = m_p / division;


		k(r) = m_p - (m_0 * m_r);

		ofile << r << "		" << real(k(r)) << endl;

	}



}

//We finally start to perform the MC method

int main() {


	//Initialize the weights

	for (int s = 0; s < q; s++) {

		W[s] = complex<double>(cos(2 * PI * s / q), sin(2 * PI * s / q));

	}



	//Initialize to the all spins=zero state

	for (int i = 0; i < N; i++) {

		S[i] = 0;

	}



	//Initialize the counters

	for (int s = 1; s < q; s++) {

		M[s] = 0;

	}

	M[0] = N;

	srand((unsigned)time(0));  //Initialize random number generation

	//We now equilibrate the matrix fliping random clusters NESTEPS time

	for (int t = 0; t < NESTEPS; t++) {

		for (int c = 0; c < NCLUSTERS; c++) {

			FlipandBuildFrom(rand() % N);

		}

	}



	//Finally we perform the measure


	//We define some variables that we will fill out with our measurements


	for (int n = 0; n < NBINS; n++) {

		complex<double> m(0., 0.);

		double m1 = 0, m2 = 0, m4 = 0;



		//Then we make NMSTEPS flips and measure

		for (int t = 0; t < NMSTEPS; t++) {

			for (int c = 0; c < NCLUSTERS; c++) {

				FlipandBuildFrom(rand() % N);

			}

			complex<double> tm(0., 0.);

			for (int s = 0; s < q; s++) {

				tm += W[s] * double(M[s]);

			}

			tm /= N;

			double tm1 = abs(tm);

			double tm2 = tm1 * tm1;

			m += tm;
			m1 += tm1;
			m2 += tm2;
			m4 += tm2 * tm2;
		}

		m /= NMSTEPS;
		m1 /= NMSTEPS;
		m2 /= NMSTEPS;
		m4 /= NMSTEPS;

		cout << m << "		" << m1 << "		" << m2 << "		" << m4 << endl;

	}

	C(k);

	return 0;

}
