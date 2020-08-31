#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

const int MATCH = 1;
const int MISMATCH = -1;
const int GAP = -2;
const int EXTEND_GAP = 0;

struct read {
	string name = "";
	string seq = "";
	int len = -1;	
};

void swap(string &a, string &b) {
	string aux = a;
	a = b;
	b = aux;
}

int smithwaterman(string &source, string &target) {
	if ((source.length()) < (target.length())) {
		swap(source, target);
	}

	vector<vector<int>> V(target.length() + 1, vector<int> (source.length() + 1, 0));

	//cout << V.size() << " " << V[0].size() << endl;

	for (int i = 1; i < V.size(); i++) {
		for (int j = 1; j < V[i].size(); j++) {
			int diag = V[i-1][j-1];
			int left = max(V[i][j-1] + GAP, V[i][j-1] + GAP + EXTEND_GAP);
			int up = max(V[i-1][j] + GAP, V[i-1][j] + GAP + EXTEND_GAP);

			if (target[i-1] == source[j-1]) {
				diag += MATCH;
			}
			else {
				diag += MISMATCH;
			}

			int valmax_aux = max(diag, left);
			int valMax = max(valmax_aux, up);

			if (valMax > 0) {
				V[i][j] = valMax;
			}
			else {
				V[i][j] = 0;
			}
		}
	}

	int vmaxvalue = -1;

	for (int i = 1; i < V.size(); i++) {
		for (int j = 1; j < V[i].size(); j++){
			if (V[i][j] > vmaxvalue) {
				vmaxvalue = V[i][j];
			}
		}
	}

	return vmaxvalue;
}

int main (void) {
	ifstream file;
	string line;
	vector<read> fichero;

	file.open("graph-contigs.fasta");

	getline(file, line);

	while (!file.eof()) {
		if (line[0] == '>') {
			read r;
			r.name = line.substr(1);

			bool condicion = true;

			while ((condicion == true) && (!file.eof())) {
				getline(file, line);

				if (line[0] != '>') {
					r.seq += line;
				}
				else {
					condicion = false;
				}
			}

			r.len = r.seq.length();

			fichero.push_back(r);
		}
	}

	file.close();

	/*
	for (int i = 0; i < fichero.size(); i++) {
		cout << fichero[i].name << ": " << fichero[i].len << endl;
	}
	*/

	ofstream outputfile;

	outputfile.open("smithwaterman.csv");
	outputfile << "source,target,sw" << endl;

	for (int i = 0; i < fichero.size(); i++) {
		for (int j = i + 1; j < fichero.size(); j++) {
			string source_name = fichero[i].name;
			string target_name = fichero[j].name;
			string source_seq = fichero[i].seq;
			string target_seq = fichero[j].seq;
			int diff = smithwaterman(source_seq, target_seq);

			outputfile << source_name << "," << target_name << "," << diff << endl;
		}
	}

	outputfile.close();

	return 0;
}