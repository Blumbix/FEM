#pragma once
#include <iostream>
using namespace std;

struct GlobalData {
	double H, W, dX, dY, Cw, Ro, k;
	int nH, nW, nN, nE, Npc;
};

struct El_Universal {
	int Npc;
	double** Pc, ** Wc;
	double** tab_N, ** dN_dKsi, ** dN_dEta;

	void fill();
	void calculate();
};

struct Node {
	int x, y;
	double T;
	bool BC;
};

struct Element {
	int ID[4];
	double** H;
	double** C;
	void calculate(Node *);
};

struct Grid {
	Node* nodes;
	Element* elements;

	void fillGrid();
};