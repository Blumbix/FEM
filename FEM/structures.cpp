#include "structures.h"
#include <iostream>
#include <vector>
using namespace std;

extern GlobalData Data;
extern El_Universal Univ;

void El_Universal::fill() {
	int n = Npc * Npc;
	Pc = new double* [n];
	Wc = new double* [n];
	tab_N = new double* [n];
	dN_dKsi = new double* [n];
	dN_dEta = new double* [n];
	for (int i = 0; i < n; i++) {
		Pc[i] = new double[Npc];
		Wc[i] = new double[Npc];
		tab_N[i] = new double[n];
		dN_dKsi[i] = new double[n];
		dN_dEta[i] = new double[n];
	}
}

void El_Universal::calculate() {
	for (int i = 0; i < Npc * Npc; i++) {
		tab_N[0][i] = 0.25 * (1 - Pc[i][0]) * (1 - Pc[i][1]);
		tab_N[1][i] = 0.25 * (1 + Pc[i][0]) * (1 - Pc[i][1]);
		tab_N[2][i] = 0.25 * (1 + Pc[i][0]) * (1 + Pc[i][1]);
		tab_N[3][i] = 0.25 * (1 - Pc[i][0]) * (1 + Pc[i][1]);

		dN_dKsi[0][i] = -0.25 * (1 - Pc[i][1]);
		dN_dKsi[1][i] = 0.25 * (1 - Pc[i][1]);
		dN_dKsi[2][i] = 0.25 * (1 + Pc[i][1]);
		dN_dKsi[3][i] = -0.25 * (1 + Pc[i][1]);

		dN_dEta[0][i] = -0.25 * (1 - Pc[i][0]);
		dN_dEta[1][i] = -0.25 * (1 + Pc[i][0]);
		dN_dEta[2][i] = 0.25 * (1 + Pc[i][0]);
		dN_dEta[3][i] = 0.25 * (1 - Pc[i][0]);
	}
}

void Grid::fillGrid() {
	nodes = new Node[Data.nN];
	elements = new Element[Data.nE];

	int n = 0;
	for (int i = 0; i < Data.nW; i++) {
		for (int j = 0; j < Data.nH; j++) {
			nodes[n].x = i;
			nodes[n].y = j;
			nodes[n].T = 100;
			if (i == 0 || i == Data.nW - 1 || j == 0 || j == Data.nH - 1) {
				nodes[n].BC = true;
			}
			else {
				nodes[n].BC = false;
			}
			n++;
		}
	}
	int x = 0;
	for (int i = 0; i < Data.nE; i++) {
		if (i % (Data.nH - 1) == 0 && i != 0) x++;
		elements[i].ID[0] = i + x;
		elements[i].ID[1] = elements[i].ID[0] + Data.nH;
		elements[i].ID[2] = elements[i].ID[1] + 1;
		elements[i].ID[3] = elements[i].ID[0] + 1;
	}

	for (int i = 0; i < Data.nE; i++) {
		elements[i].calculate(nodes);
	}
}


void Element::calculate(Node *nodes) {
	int n = Univ.Npc * Univ.Npc;
	double*** tempH = new double** [n];
	double*** tempC = new double** [n];
	double** jacob = new double* [Univ.Npc];
	double* dX_dKsi = new double[n];
	double* dY_dEta = new double[n];
	double* dX_dEta = new double[n];
	double* dY_dKsi = new double[n];
	double* detJ = new double[n];
	double** dN_dX = new double* [n];
	double** dN_dY = new double* [n];
	H = new double* [n];
	C = new double* [n];

	for (int i = 0; i < n; i++) {
		H[i] = new double[n];
		C[i] = new double[n];
		dN_dX[i] = new double[n];
		dN_dY[i] = new double[n];
		tempH[i] = new double* [n];
		tempC[i] = new double* [n];
		for (int j = 0; j < n; j++) {
			tempH[i][j] = new double[n];
			tempC[i][j] = new double[n];
		}
	}

	for (int i = 0; i < Univ.Npc; i++) {
		jacob[i] = new double[Univ.Npc];
	}

	for (int i = 0; i < n; i++) {
		dX_dKsi[i] = Univ.dN_dKsi[0][i] * (nodes[ID[0]].x * Data.dX) +
			Univ.dN_dKsi[1][i] * (nodes[ID[1]].x * Data.dX) +
			Univ.dN_dKsi[2][i] * (nodes[ID[2]].x * Data.dX) +
			Univ.dN_dKsi[3][i] * (nodes[ID[3]].x * Data.dX);
		dX_dEta[i] = Univ.dN_dEta[0][i] * (nodes[ID[0]].x * Data.dX) +
			Univ.dN_dEta[1][i] * (nodes[ID[1]].x * Data.dX) +
			Univ.dN_dEta[2][i] * (nodes[ID[2]].x * Data.dX) +
			Univ.dN_dEta[3][i] * (nodes[ID[3]].x * Data.dX);
		dY_dEta[i] = Univ.dN_dEta[0][i] * (nodes[ID[0]].y * Data.dY) +
			Univ.dN_dEta[1][i] * (nodes[ID[1]].y * Data.dY) +
			Univ.dN_dEta[2][i] * (nodes[ID[2]].y * Data.dY) +
			Univ.dN_dEta[3][i] * (nodes[ID[3]].y * Data.dY);
		dY_dKsi[i] = Univ.dN_dKsi[0][i] * (nodes[ID[0]].y * Data.dY) +
			Univ.dN_dKsi[1][i] * (nodes[ID[1]].y * Data.dY) +
			Univ.dN_dKsi[2][i] * (nodes[ID[2]].y * Data.dY) +
			Univ.dN_dKsi[3][i] * (nodes[ID[3]].y * Data.dY);
	}

	for (int i = 0; i < n; i++) {
		detJ[i] = dX_dKsi[i] * dY_dEta[i] - dX_dEta[i] * dY_dKsi[i];
		for (int j = 0; j < n; j++) {
			dN_dX[i][j] = (1 / detJ[i]) * ((dY_dEta[i] * Univ.dN_dKsi[i][j]) + (-dY_dKsi[i] * Univ.dN_dEta[i][j]));
			dN_dY[i][j] = (1 / detJ[i]) * ((-dX_dEta[i] * Univ.dN_dKsi[i][j]) + (dX_dKsi[i] * Univ.dN_dEta[i][j]));
		}
	}

	for (int m = 0; m < n; m++) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				tempH[m][i][j] = Data.k * (dN_dX[i][m] * dN_dX[j][m] + dN_dY[i][m] * dN_dY[j][m]);
				tempC[m][i][j] = Data.Cw * Data.Ro * Univ.tab_N[i][m] * Univ.tab_N[j][m];
			}
		}
	}

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			double sumH = 0, sumC = 0;
			for (int m = 0; m < n; m++) {
				sumH += tempH[m][i][j] * Univ.Wc[m][0] * Univ.Wc[m][1] * detJ[m];
				sumC += tempC[m][i][j] * Univ.Wc[m][0] * Univ.Wc[m][1] * detJ[m];
			}
			H[i][j] = sumH;
			C[i][j] = sumC;
		}
	}
}