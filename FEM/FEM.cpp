#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include "structures.h"

using namespace std;

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
	file >> Univ.Cw;
	file >> Univ.Ro;
	Univ.fill();
	for (int i = 0; i < Univ.Npc*Univ.Npc; i++) {
		file >> Univ.Pc[i][0];
		file >> Univ.Wc[i][0];
		file >> Univ.Pc[i][1];
		file >> Univ.Wc[i][1];
	}
	file.close();

	Univ.calculate();
	Data.nN = Data.nH * Data.nW;
	Data.nE = (Data.nH - 1) * (Data.nW - 1);
	Data.dX = Data.W / (Data.nW - 1);
	Data.dY = Data.H / (Data.nH - 1);
}

int main()
{
	readFile();
	Grid g;
	g.fillGrid();

	cout << "24. wezel:\n";
	cout << "x: " << g.nodes[23].x << "\ty: " << g.nodes[23].y << endl;
	cout << "Granica: " << g.nodes[23].BC << endl;
	int x = 0;
	cout << "\n1. element:\n";
	cout << g.elements[x].ID[0] << " " << g.elements[x].ID[1] << " " << g.elements[x].ID[2] << " " << g.elements[x].ID[3] << endl;
	cout << "\n2. punkt calkowania, ksi i eta oraz ich wagi:\n";
	cout << Univ.Pc[1][0] << "  " << Univ.Wc[1][0] << "\n" << Univ.Pc[1][1] << "  " << Univ.Wc[1][1] << endl;
	cout << "\nN1 dla 1. pkt.calk.:\n" << Univ.tab_N[0][0] << endl;
	cout << "\ndN1 po dKsi dla 1. pkt.calk.:\n" << Univ.dN_dKsi[0][0];
	cout << "\ndN1 po dEta dla 1. pkt.calk.:\n" << Univ.dN_dEta[0][0];

	cout << "\n1. Node, H[0][1]:\t" << g.nodes[0].H[0][1];
	cout << "\n1. Node, C[0][1]:\t" << g.nodes[0].C[0][1] << endl;
	return 0;
}
