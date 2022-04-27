#include <armadillo>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <math.h>
#include <complex>

using namespace arma;
using namespace std;

#define PI 3.141592653589793238462643

vector <double> TJ(double start_in, double end_in, int num_in)
{

	std::vector<double> linspaced;

	double start = static_cast<double>(start_in);
	double end = static_cast<double>(end_in);
	double num = static_cast<double>(num_in);

	if (num == 0) { return linspaced; }
	if (num == 1)
	{
		linspaced.push_back(start);
		return linspaced;
	}

	double delta = (end - start) / (num - 1);

	for (int i = 0; i < num - 1; ++i)
	{
		linspaced.push_back(start + delta * i);
	}
	linspaced.push_back(end); // I want to ensure that start and end are exactly the same as the input

	return linspaced;
}

//First, we introduce the desired settings of the system that we want to study

const int q = 3;  //Number of values each spin can take

const int L = 16;  //System size		NEEDS TO BE CHANGED FOR DIFFERENT SIZES L={8, 16, 32}

const int N = L;  //Total number of spins

const int NCLUSTERS = 3;  //Number of cluster builds in one MC step

const int NESTEPS = 10000;  //Number of equilibrium MC steps

const int NMSTEPS = 10000;  //Number of measurement MC steps

const int NBINS = 10;  //Number of measurement bins

const int points = pow(10, 3);  //Number of points

Mat <int> S(N, N);  //The spin matrix

vec M(q);  //Number of spins in each state

cx_vec W(q);  //Order parameter weights

vector <int> pos(N);

vector <double> T = TJ(0, pow(10, 5), points);  // NEEDS TO BE CHANGED FOR EVALUATING THE GAMMA RATIO TO A MAX. T = 2

vector <double> pconnect;


//Now we will define some methods in order to manipulate the lattice

enum dirs { RIGHT, LEFT, DOWN, UP };

int indx(int x) {

	return x;  //Make an indx on every site

}


int xpos(int i) {

	return i % L;

}


int Nbr(int i, int j, int dir) {

	if (dir == 0) {

		pos[0] = i;

		pos[1] = (xpos(j) + 1) % L;

	}


	else if (dir == 1) {

		pos[0] = i;

		pos[1] = (xpos(j) - 1 + L) % L;

	}


	else if (dir == 2) {

		pos[0] = (xpos(i) + 1) % L;

		pos[1] = j;

	}


	else {

		pos[0] = (xpos(i) - 1 + L) % L;

		pos[1] = j;

	}

	return 0;

}


void FlipandBuildFrom(int s, int u, int k) {

	int oldstate(S[s, u]), newstate((S[s, u] + 1) % q);

	S[s, u] = newstate;  //Flip spin

	M[oldstate]--; M[newstate]++;  //Update spin counts

	for (int dir = 0; dir < 4; dir++) {

		Nbr(s, u, dir);

		int i = pos[0];

		int j = pos[1];

		if (S[i, j] == oldstate) {

			if (rand() / (RAND_MAX + 1.) < pconnect[k]) {

				FlipandBuildFrom(i, j, k);

			}

		}

	}

}


cx_vec k(N);


void C(cx_vec& k, int h) {

	ofstream ofile;
	ofile.open("C(r)_2D.txt");
	ofile << "r" << "	" << "Real part" << endl;
	ofile << scientific;

	for (int r = 0; r < N; r++) {

		complex<double>m_0(0., 0.);
		complex<double>m_r(0., 0.);
		complex<double>m_p(0., 0.);

		for (int t = 0; t < NMSTEPS; t++) {


			for (int c = 0; c < NCLUSTERS; c++) {

				FlipandBuildFrom(rand() % N, rand() % N, h);

			}


			m_0 = m_0 + conj(complex<double>(cos(2 * PI * S(0,0) / q), sin(2 * PI * S(0, 0) / q)));
			m_r = m_r + complex<double>(cos(2 * PI * S(r, r) / q), sin(2 * PI * S(r, r) / q));
			m_p = m_p + conj(complex<double>(cos(2 * PI * S(0, 0) / q), sin(2 * PI * S(0, 0) / q))) * complex<double>(cos(2 * PI * S(r, r) / q), sin(2 * PI * S(r, r) / q));

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

	ofstream ofile;
	ofile.open("m_2D.txt");//NEEDS TO BE CHANGED FOR DIFFERENT SIZES L={8, 16, 32} WHEN EVALUATING GAMMA WITH NAME "m_2D_L.txt"
	ofile << "T" << "	" << "Real part of <m>" << "	" << "<|m|^2>" << "	" << "<|m|^4>" << endl;
	ofile << scientific;
	

	for (int i = 0; i < points; i++) {

		pconnect.push_back(1 - exp(-1 / T[i]));

	}

	for (int h = 0; h < points; h++) {


		//Initialize the weights

		for (int s = 0; s < q; s++) {

			W(s) = complex<double>(cos(2 * PI * s / q), sin(2 * PI * s / q));

		}



		//Initialize to the all spins=zero state (SUSTITUUIR POR CREAR DIRECTAMENTE UNA MATRIZ DE CEROS)

		for (int i = 0; i < N; i++) {

			for (int j = 0; j < N; j++) {

				S[i, j] = 0;

			}

		}



		//Initialize the counters

		for (int s = 1; s < q; s++) {

			M[s] = 0;

		}

		M[0] = N * N;



		srand((unsigned)time(0));  //Initialize random number generation



		//We now equilibrate the matrix fliping random clusters NESTEPS time

		for (int t = 0; t < NESTEPS; t++) {

			for (int c = 0; c < NCLUSTERS; c++) {

				FlipandBuildFrom(rand() % N, rand() % N, h);

			}

		}


		//We define some variables that we will fill out with our measurements

		complex<double> m(0., 0.);

		double m1 = 0, m2 = 0, m4 = 0;


		//Finally we perform the measure

		for (int n = 0; n < NBINS; n++) {


			//Then we make NMSTEPS flips and measure

			for (int t = 0; t < NMSTEPS; t++) {

				for (int c = 0; c < NCLUSTERS; c++) {

					FlipandBuildFrom(rand() % N, rand() % N, h);

				}

				complex<double> tm(0., 0.);

				for (int s = 0; s < q; s++) {

					tm += W[s] * double(M[s]);

				}

				tm /= N*N;

				m += tm;
				m1 += abs(tm);
				m2 += abs(tm) * abs(tm);
				m4 += pow(abs(tm), 4);

			}

			m /= NMSTEPS;
			m1 /= NMSTEPS;
			m2 /= NMSTEPS;
			m4 /= NMSTEPS;

			std::cout << m << "		" << m1 << "		" << m2 << "		" << m4 << endl;

		}


		ofile << T[h] << "		" << real(m) << "		" << m2 << "		" << m4 << endl;

		//C(k, h);

	}

	

	return 0;

}
