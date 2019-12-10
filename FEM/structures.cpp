#include "structures.h"
#include <iostream>
using namespace std;

extern GlobalData Data;
extern El_Universal Univ;

void El_Universal::fill() {
	int n = Npc * Npc;
	Pc = new double* [n];
	Wc = new double* [n];
	Pc_BC = new double* [2*n];
	Wc_BC = new double* [2*n];
	tab_N = new double* [n];
	tab_N_BC = new double* [n];
	dN_dKsi = new double* [n];
	dN_dEta = new double* [n];
	for (int i = 0; i < n; i++) {
		Pc[i] = new double[Npc];
		Wc[i] = new double[Npc];
		tab_N[i] = new double[n];
		tab_N_BC[i] = new double[2 * n];
		dN_dKsi[i] = new double[n];
		dN_dEta[i] = new double[n];
	}
	for (int i = 0; i < 2*n; i++) {
		Pc_BC[i] = new double[Npc];
		Wc_BC[i] = new double[Npc];
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
	for (int i = 0; i < 2 * Npc * Npc; i++) {
		tab_N_BC[0][i] = 0.25 * (1 - Pc_BC[i][0]) * (1 - Pc_BC[i][1]);
		tab_N_BC[1][i] = 0.25 * (1 + Pc_BC[i][0]) * (1 - Pc_BC[i][1]);
		tab_N_BC[2][i] = 0.25 * (1 + Pc_BC[i][0]) * (1 + Pc_BC[i][1]);
		tab_N_BC[3][i] = 0.25 * (1 - Pc_BC[i][0]) * (1 + Pc_BC[i][1]);
	}
}

void Grid::fillGrid() {
	nodes = new Node[Data.nN];
	elements = new Element[Data.nE];

	int n = 0;
	for (int i = 0; i < Data.nW; i++) {
		for (int j = 0; j < Data.nH; j++) {
			nodes[n].x = i * Data.dX;
			nodes[n].y = j * Data.dY;
			nodes[n].T = Data.T0;
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
	// Alokacja danych
	int n = Univ.Npc * Univ.Npc;
	double*** tempH = new double** [n];
	double*** tempH_BC = new double** [2*n];
	double*** tempC = new double** [n];
	double** jacob = new double* [Univ.Npc];
	double* dX_dKsi = new double[n];
	double* dY_dEta = new double[n];
	double* dX_dEta = new double[n];
	double* dY_dKsi = new double[n];
	double* detJ = new double[n];
	double** dN_dX = new double* [n];
	double** dN_dY = new double* [n];
	double*** temp_BC = new double**[n];
	double** temp_P = new double* [n];
	H = new double* [n];
	C = new double* [n];
	P = new double[n];

	for (int i = 0; i < n; i++) {
		H[i] = new double[n];
		C[i] = new double[n];
		dN_dX[i] = new double[n];
		dN_dY[i] = new double[n];
		tempH[i] = new double* [n];
		tempC[i] = new double* [n];
		temp_BC[i] = new double* [n];
		temp_P[i] = new double[n];
		for (int j = 0; j < n; j++) {
			tempH[i][j] = new double[n];
			tempC[i][j] = new double[n];
			temp_BC[i][j] = new double[n];
		}
	}

	// H dla kazdego punktu calkowania
	for (int i = 0; i < 2 * n; i++) {
		tempH_BC[i] = new double* [n];
		for (int j = 0; j < n; j++) {
			tempH_BC[i][j] = new double[n];
		}
	}

	for (int i = 0; i < Univ.Npc; i++) {
		jacob[i] = new double[Univ.Npc];
	}

	for (int i = 0; i < n; i++) {
		dX_dKsi[i] = Univ.dN_dKsi[0][i] * nodes[ID[0]].x +
			Univ.dN_dKsi[1][i] * nodes[ID[1]].x +
			Univ.dN_dKsi[2][i] * nodes[ID[2]].x +
			Univ.dN_dKsi[3][i] * nodes[ID[3]].x;
		dX_dEta[i] = Univ.dN_dEta[0][i] * nodes[ID[0]].x +
			Univ.dN_dEta[1][i] * nodes[ID[1]].x +
			Univ.dN_dEta[2][i] * nodes[ID[2]].x +
			Univ.dN_dEta[3][i] * nodes[ID[3]].x;
		dY_dEta[i] = Univ.dN_dEta[0][i] * nodes[ID[0]].y +
			Univ.dN_dEta[1][i] * nodes[ID[1]].y +
			Univ.dN_dEta[2][i] * nodes[ID[2]].y +
			Univ.dN_dEta[3][i] * nodes[ID[3]].y;
		dY_dKsi[i] = Univ.dN_dKsi[0][i] * nodes[ID[0]].y +
			Univ.dN_dKsi[1][i] * nodes[ID[1]].y +
			Univ.dN_dKsi[2][i] * nodes[ID[2]].y +
			Univ.dN_dKsi[3][i] * nodes[ID[3]].y;
	}

	for (int i = 0; i < n; i++) {
		detJ[i] = dX_dKsi[i] * dY_dEta[i] - dX_dEta[i] * dY_dKsi[i];
		for (int j = 0; j < n; j++) {
			dN_dX[i][j] = (1 / detJ[i]) * ((dY_dEta[i] * Univ.dN_dKsi[i][j]) + (-dY_dKsi[i] * Univ.dN_dEta[i][j]));
			dN_dY[i][j] = (1 / detJ[i]) * ((-dX_dEta[i] * Univ.dN_dKsi[i][j]) + (dX_dKsi[i] * Univ.dN_dEta[i][j]));
		}
	}

	for (int pc = 0; pc < n; pc++) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				tempH[pc][i][j] = Data.k * (dN_dX[i][pc] * dN_dX[j][pc] + dN_dY[i][pc] * dN_dY[j][pc]);
				tempC[pc][i][j] = Data.Cw * Data.Ro * Univ.tab_N[i][pc] * Univ.tab_N[j][pc];
			}
		}
	}

	// Warunek brzegowy
	for (int pc = 0; pc < 2 * n; pc++) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				tempH_BC[pc][i][j] = Data.a * Univ.tab_N_BC[i][pc] * Univ.tab_N_BC[j][pc];
			}
		}
	}

	if (nodes[ID[0]].BC && nodes[ID[1]].BC) {
		for (int i = 0; i < n; i++) {
			temp_P[0][i] = (Univ.tab_N_BC[i][0] * Univ.Wc_BC[0][0] + Univ.tab_N_BC[i][1] * Univ.Wc_BC[1][0]) * ((nodes[ID[1]].x - nodes[ID[0]].x) / 2);
			for (int j = 0; j < n; j++) {
				temp_BC[0][i][j] = (tempH_BC[0][i][j] * Univ.Wc_BC[0][0] + tempH_BC[1][i][j] * Univ.Wc_BC[1][0]) * ((nodes[ID[1]].x - nodes[ID[0]].x) / 2);

			}
		}
	} else {
		for (int i = 0; i < n; i++) {
			temp_P[0][i] = 0;
			for (int j = 0; j < n; j++) {
				temp_BC[0][i][j] = 0;
			}
		}
	}
	if (nodes[ID[1]].BC && nodes[ID[2]].BC) {
		for (int i = 0; i < n; i++) {
			temp_P[1][i] = (Univ.tab_N_BC[i][2] * Univ.Wc_BC[2][1] + Univ.tab_N_BC[i][3] * Univ.Wc_BC[3][1]) * ((nodes[ID[2]].y - nodes[ID[1]].y) / 2);
			for (int j = 0; j < n; j++) {
				temp_BC[1][i][j] = (tempH_BC[2][i][j] * Univ.Wc_BC[2][1] + tempH_BC[3][i][j] * Univ.Wc_BC[3][1]) * ((nodes[ID[2]].y - nodes[ID[1]].y) / 2);
			}
		}
	}
	else {
		for (int i = 0; i < n; i++) {
			temp_P[1][i] = 0;
			for (int j = 0; j < n; j++) {
				temp_BC[1][i][j] = 0;
			}
		}
	}
	if (nodes[ID[2]].BC && nodes[ID[3]].BC) {
		for (int i = 0; i < n; i++) {
			temp_P[2][i] = (Univ.tab_N_BC[i][4] * Univ.Wc_BC[4][0] + Univ.tab_N_BC[i][5] * Univ.Wc_BC[5][0]) * ((nodes[ID[2]].x - nodes[ID[3]].x) / 2);
			for (int j = 0; j < n; j++) {
				temp_BC[2][i][j] = (tempH_BC[4][i][j] * Univ.Wc_BC[4][0] + tempH_BC[5][i][j] * Univ.Wc_BC[5][0]) * ((nodes[ID[2]].x - nodes[ID[3]].x) / 2);
			}
		}
	}
	else {
		for (int i = 0; i < n; i++) {
			temp_P[2][i] = 0;
			for (int j = 0; j < n; j++) {
				temp_BC[2][i][j] = 0;
			}
		}
	}
	if (nodes[ID[3]].BC && nodes[ID[0]].BC) {
		for (int i = 0; i < n; i++) {
			temp_P[3][i] = (Univ.tab_N_BC[i][6] * Univ.Wc_BC[6][1] + Univ.tab_N_BC[i][7] * Univ.Wc_BC[7][1]) * ((nodes[ID[3]].y - nodes[ID[0]].y) / 2);
			for (int j = 0; j < n; j++) {
				temp_BC[3][i][j] = (tempH_BC[6][i][j] * Univ.Wc_BC[6][1] + tempH_BC[7][i][j] * Univ.Wc_BC[7][1]) * ((nodes[ID[3]].y - nodes[ID[0]].y) / 2);
			}
		}
	}
	else {
		for (int i = 0; i < n; i++) {
			temp_P[3][i] = 0;
			for (int j = 0; j < n; j++) {
				temp_BC[3][i][j] = 0;
			}
		}
	}

	// Zliczanie wynikow do glownych zmiennych
	for (int i = 0; i < n; i++) {
		double sumP = 0;
		for (int j = 0; j < n; j++) {
			double sumH = 0, sumC = 0;
			for (int pc = 0; pc < n; pc++) {
				sumH += (tempH[pc][i][j] * Univ.Wc[pc][0] * Univ.Wc[pc][1] * detJ[pc]) + temp_BC[pc][i][j];
				//sumH += (tempH[pc][i][j] * Univ.Wc[pc][0] * Univ.Wc[pc][1] * detJ[pc]);
				sumC += tempC[pc][i][j] * Univ.Wc[pc][0] * Univ.Wc[pc][1] * detJ[pc];
			}
			H[i][j] = sumH;
			C[i][j] = sumC;

			sumP += Data.a * temp_P[j][i] * Data.ambT;
		}
		P[i] = sumP;
	}
}

void Grid::calculate() {
	int n = Data.nN, ID[4];
	H = new double* [n];
	C = new double* [n];
	P = new double[n];

	for (int i = 0; i < n; i++) {
		H[i] = new double[n];
		C[i] = new double[n];
	}

	for (int i = 0; i < n; i++) {
		P[i] = 0;
		for (int j = 0; j < n; j++) {
			H[i][j] = 0;
			C[i][j] = 0;
		}
	}

	for (int i = 0; i < Data.nE; i++) {
		for (int j = 0; j < 4; j++) {
			ID[j] = elements[i].ID[j];
		}

		for (int x = 0; x < 4; x++) {
			P[ID[x]] += elements[i].P[x];
			for (int y = 0; y < 4; y++) {
				H[ID[x]][ID[y]] += elements[i].H[x][y];
				C[ID[x]][ID[y]] += elements[i].C[x][y];
			}
		}
	}

}