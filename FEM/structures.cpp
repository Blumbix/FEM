#include "structures.h"
#include <iostream>
#include <vector>
using namespace std;

extern GlobalData Data;
extern El_Universal Univ;

void El_Universal::fill() {
	Pc = new double* [Npc * Npc];
	Wc = new double* [Npc * Npc];
	tab_N = new double* [Npc * Npc];
	dN_dKsi = new double* [Npc * Npc];
	dN_dEta = new double* [Npc * Npc];
	for (int i = 0; i < Npc * Npc; i++) {
		Pc[i] = new double[2];
		Wc[i] = new double[2];
		tab_N[i] = new double[Npc * Npc];
		dN_dKsi[i] = new double[Npc * Npc];
		dN_dEta[i] = new double[Npc * Npc];
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
	int x = 1;
	for (int i = 0; i < Data.nE; i++) {
		if (i % (Data.nH - 1) == 0 && i != 0) x++;
		elements[i].ID[0] = i + x;
		elements[i].ID[1] = elements[i].ID[0] + Data.nH;
		elements[i].ID[2] = elements[i].ID[1] + 1;
		elements[i].ID[3] = elements[i].ID[0] + 1;
	}
}