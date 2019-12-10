#pragma once
#include <iostream>
using namespace std;

struct GlobalData {
	double H, W, dX, dY, Cw, Ro, k, a, T0, dTau, time, ambT;
	int nH, nW, nN, nE, Npc;
};

struct El_Universal {
	int Npc;
	double** Pc, ** Wc, ** Pc_BC, ** Wc_BC;
	double** tab_N, ** dN_dKsi, ** dN_dEta, ** tab_N_BC;

	void fill();
	void calculate();
};

struct Node {
	double x, y;
	double T;
	bool BC;
};

struct Element {
	int ID[4];
	double** H;
	double** C;
	double* P;
	void calculate(Node *);
};

struct Grid {
	Node* nodes;
	Element* elements;
	double** H;
	double** C;
	double* P;

	void fillGrid();
	void calculate();
};