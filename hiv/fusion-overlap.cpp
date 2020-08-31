#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

const int ERRORS = 0;

struct read {
	string name = "";
	string seq = "";
	int len = -1;
	int pos = -1;	
};

bool comparator(const read& lhs, const read& rhs) {
	return lhs.len > rhs.len;
}

int differences(const string &s, const string &t) {
	int d = 0;

	for (int i = 0; i < s.length(); i++) {
		if (s[i] != t[i]) {
			d++;
		}
	}

	return d;
}

int main(void) {
	ifstream inputfile;
	vector<read> secuencias;

	inputfile.open("contigs.sam");

	string line;	
	char delimiter = '\t';

	while (!inputfile.eof()) {
		getline(inputfile, line);

		if (line[0] != '@') {
			vector<string> str;
			string acc = "";

			for (int i = 0; i <= line.size(); i++) {
				if ((line[i] == delimiter) || (i == line.size())) {
			        str.push_back(acc);
			        acc = "";
		    	}
		    	else {
		        	acc += line[i];
	    		}

			}

			int flag = stoi(str[1]);
			int pos = stoi(str[3]);

			if ((pos != 0) && ((flag == 0) || (flag == 16))) {
				read r;
				r.name = str[0];
				r.seq = str[9];
				r.len = r.seq.length();
				r.pos = pos;

				secuencias.push_back(r);
			}
		}
	}

	cout << secuencias.size() << " sequences loaded!" << endl;

	bool iter = true; //true
	int iters = -1;

	while (iter) {
		iters++;

		cout << "iteration " << iters << endl;
		bool fusion = false;

		sort(secuencias.begin(), secuencias.end(), &comparator);

		cout << secuencias[0].len << " " << secuencias[1].len << " " << secuencias[2].len << endl;

		for (int i = 0; i < secuencias.size(); i++) {
			if (secuencias[i].pos != -1) {
				for (int j = i + 1; j < secuencias.size(); j++) {
					if (secuencias[j].pos != -1) {
						string base_seq = secuencias[i].seq;
						string next_seq = secuencias[j].seq;

						int base_start = secuencias[i].pos;
						int base_end = base_start + secuencias[i].len;

						int next_start = secuencias[j].pos;
						int next_end = next_start + secuencias[j].len;

						int index_base_start = -1;
						int index_base_end = -1;
						int index_next_start = -1;
						int index_next_end = -1;

						if ((next_end >= base_start) && (next_end <= base_end) && (next_start <= base_start)) {
							index_next_start = base_start - next_start;
                            index_next_end = next_seq.length();

                            index_base_start = 0;
                            index_base_end = index_next_end - index_next_start;

                            int size_base = index_base_end - index_base_start;
	                        int size_next = index_next_end - index_next_start;                       	

	                        int diff = differences(base_seq.substr(index_base_start, size_base), next_seq.substr(index_next_start, size_next));

                            if (diff == ERRORS) {
                            	fusion = true;

                            	string newseq = next_seq + base_seq.substr(index_base_end, base_seq.length());
	                            int lennewseq = newseq.length();
	                            int indexnewseq = next_start;

	                            secuencias[i].seq = newseq;
	                            secuencias[i].len = lennewseq;
	                            secuencias[i].pos = indexnewseq;

	                        	secuencias[j].seq = "";
	                            secuencias[j].len = -1;
	                            secuencias[j].pos = -1;
                            }
						}
						else {
							if ((base_start <= next_start) && (base_end >= next_end)) {
								index_next_start = 0;
	                            index_next_end = next_seq.length();


	                            index_base_start = next_start - base_start;
	                            index_base_end = index_base_start + next_seq.length();

	                            int size_base = index_base_end - index_base_start;
	                            int size_next = index_next_end - index_next_start;                       	

	                            int diff = differences(base_seq.substr(index_base_start, size_base), next_seq.substr(index_next_start, size_next));

	                            if (diff == ERRORS) {
	                            	fusion = true;

	                            	string newseq = base_seq;
		                            int lennewseq = newseq.length();
		                            int indexnewseq = base_start;

		                            secuencias[i].seq = newseq;
		                            secuencias[i].len = lennewseq;
		                            secuencias[i].pos = indexnewseq;

		                        	secuencias[j].seq = "";
		                            secuencias[j].len = -1;
		                            secuencias[j].pos = -1;
								}
							} else {
								if ((base_start <= next_start) && (base_end >= next_start) && (base_end <= next_end)) {
									index_next_start = 0;
		                            index_next_end = base_end - next_start;

		                            index_base_start = base_seq.length() - index_next_end;
		                            index_base_end = base_seq.length();

		                            int size_base = index_base_end - index_base_start;
			                        int size_next = index_next_end - index_next_start;                       	

			                        int diff = differences(base_seq.substr(index_base_start, size_base), next_seq.substr(index_next_start, size_next));

		                            if (diff == ERRORS) {
		                            	fusion = true;

		                            	string newseq = base_seq + next_seq.substr(index_next_end, next_seq.length());
			                            int lennewseq = newseq.length();
			                            int indexnewseq = base_start;

			                            secuencias[i].seq = newseq;
			                            secuencias[i].len = lennewseq;
			                            secuencias[i].pos = indexnewseq;

			                        	secuencias[j].seq = "";
			                            secuencias[j].len = -1;
			                            secuencias[j].pos = -1;
									}
								}
							}
						}
					}
				}
			}
		}

		if (!fusion) {
			iter = false;
		}
	}

	int count = 0;

	for (int i = 0; i < secuencias.size(); i++) {
		if (secuencias[i].pos != -1) {
			count++;
		}
	}
	
	cout << count << " final sequences" << endl;

	ofstream outputfile;

	outputfile.open("final-contigs.fasta");

	for (int i = 0; i < secuencias.size(); i++) {
		ouputfile << ">" << secuencias[i].name << endl;

		if (i == secuencias.size() - 1){
			outputfile << secuencias[i].seq;
		}
		else {
			outputfile << secuencias[i].seq << endl;
		}
	}

	return 0;
}