#pragma once
#include <iostream>
using namespace std;

struct GlobalData {
	float H, W, dX, dY;
	int nH, nW, nN, nE, Npc;
};

struct El_Universal {
	int Npc;
	double Cw, Ro;
	double** Pc, ** Wc;
	double** tab_N, ** dN_dKsi, ** dN_dEta;

	void fill();
	void calculate();
};

struct Node {
	int x, y;
	double T;
	bool BC;

	double** H;
	double** C;
	void calculate();
};

struct Element {
	int ID[4];
};

struct Grid {
	Node* nodes;
	Element* elements;

	void fillGrid();
};