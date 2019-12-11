#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include "structures.h"

using namespace std;
const double eps = 1e-12;

GlobalData Data;
El_Universal Univ;

void readFile() {
	ifstream file;
	file.open("input.txt");
	file >> Data.H;
	file >> Data.W;
	file >> Data.nH;
	file >> Data.nW;
	file >> Univ.Npc;
	file >> Data.Cw;
	file >> Data.Ro;
	file >> Data.k;
	file >> Data.a;
	file >> Data.T0;
	file >> Data.ambT;
	file >> Data.time;
	file >> Data.dTau;

	Univ.fill();
	for (int i = 0; i < Univ.Npc*Univ.Npc; i++) {
		file >> Univ.Pc[i][0];
		file >> Univ.Wc[i][0];
		file >> Univ.Pc[i][1];
		file >> Univ.Wc[i][1];
	}
	for (int i = 0; i < 2 * Univ.Npc * Univ.Npc; i++) {
		file >> Univ.Pc_BC[i][0];
		file >> Univ.Wc_BC[i][0];
		file >> Univ.Pc_BC[i][1];
		file >> Univ.Wc_BC[i][1];
	}
	file.close();

	Univ.calculate();
	Data.nN = Data.nH * Data.nW;
	Data.nE = (Data.nH - 1) * (Data.nW - 1);
	Data.dX = Data.W / (double(Data.nW) -1);
	Data.dY = Data.H / (double(Data.nH) -1);
}

double* gauss(int n, double** AB, double* X) {
	int i, j, k;
	double m, s;

	// eliminacja wspó³czynników
	for (i = 0; i < n - 1; i++)	{
		for (j = i + 1; j < n; j++)		{
			if (fabs(AB[i][i]) < eps) return NULL;
			m = -AB[j][i] / AB[i][i];
			for (k = i + 1; k <= n; k++)
				AB[j][k] += m * AB[i][k];
		}
	}

	// wyliczanie niewiadomych
	for (i = n - 1; i >= 0; i--){
		s = AB[i][n];
		for (j = n - 1; j >= i + 1; j--)
			s -= AB[i][j] * X[j];
		if (fabs(AB[i][i]) < eps) return NULL;
		X[i] = s / AB[i][i];
	}
	return X;
}

void step(Grid &g, double* t0) {
	double** H_roof = new double* [Data.nN];
	double** X = new double* [Data.nN];
	for (int i = 0; i < Data.nN; i++) {
		H_roof[i] = new double[Data.nN];
		X[i] = new double[Data.nN];
	}
	double* P_roof = new double[Data.nN];

	for (int i = 0; i < Data.nN; i++) {
		double C_new = 0;
		for (int j = 0; j < Data.nN; j++) {
			H_roof[i][j] = g.H[i][j] + (g.C[i][j] / Data.dTau);
			C_new += g.C[i][j] / Data.dTau * t0[j];
		}
		P_roof[i] = g.P[i] + C_new;
	}

	/*cout << "\nMacierz H dach:\n";
	for (int i = 0; i < Data.nN; i++) {
		for (int j = 0; j < Data.nN; j++) {
			cout << setw(10) <<  H_roof[i][j] << "  ";
		}
		cout << "\n";
	}*/

	cout << "\nWektor P dach:\n";
	for (int i = 0; i < Data.nN; i++) {
		cout << setw(10) << P_roof[i] << "  ";
	}
	cout << "\n\n";

	/* METODA ELIMINACJI GAUSSA */
	double** AB = new double* [Data.nN];
	for (int j = 0; j < Data.nN; j++) {
		AB[j] = new double[Data.nN + 1];
	}
	for (int i = 0; i < Data.nN; i++) {
		for (int j = 0; j < Data.nN; j++) {
			AB[i][j] = H_roof[i][j];
		}
	}
	for (int i = 0; i < Data.nN; i++) {
		AB[i][16] = P_roof[i];
	}

	t0 = gauss(Data.nN, AB, t0);

	double t_min = t0[0];
	double t_max = t0[0];
	for (int i = 1; i < Data.nN; i++) {
		if (t_min > t0[i])
			t_min = t0[i];
		if (t_max < t0[i])
			t_max = t0[i];
	}
	cout << "T min: " << t_min << "\nT max: " << t_max << "\n";
}

int main()
{
	readFile();
	Grid g;
	g.fillGrid();
	g.calculate();
	
	cout << "\nMacierz H globalna:\n";
	for (int i = 0; i < Data.nN; i++) {
		for (int j = 0; j < Data.nN; j++) {
			cout << setw(10) << g.H[i][j] << "  ";
		}
		cout << "\n";
	}

	cout << "\nMacierz C globalna:\n";
	for (int i = 0; i < Data.nN; i++) {
		for (int j = 0; j < Data.nN; j++) {
			cout << setw(10) << g.C[i][j] << "  ";
		}
		cout << "\n";
	}

	cout << "\nWektor P globalny:\n";
	for (int i = 0; i < Data.nN; i++) {
		cout << setw(10) << g.P[i] << "  ";
	}
	cout << "\n\n";

	double* t = new double[Data.nN];
	for (int i = 0; i < Data.nN; i++)
		t[i] = Data.T0;

	for (int i = Data.dTau; i <= Data.time; i+=Data.dTau) {
		cout << "\n\nIteracja "<<i/Data.dTau<<". czas: "<<i<<"s\n";
		step(g, t);
	}

	return 0;
}
