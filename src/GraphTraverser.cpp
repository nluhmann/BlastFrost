/*
 * GraphTraverser.cpp
 *
 *  Created on: 26 Feb 2018
 *      Author: nina
 */

#include "GraphTraverser.hpp"

using namespace std;






GraphTraverser::GraphTraverser(ColoredCDBG<UnitigData>& graph) :
	cdbg(graph) {
	cout << "GraphTraverser initialized!" << endl;
}




//unordered_map<size_t,vector<int>> GraphTraverser::search(string query, int k) const{
//
//	vector<Kmer> kmers;
//
//	//split query into sequence of kmers (!!! query can contain a kmer multiple times, do not change the order of the kmers at this point!)
//	for(int i = 0; i< (query.length()-k+1); ++i){
//	    const string kmer = query.substr(i,k);
//	    const char *cstr = kmer.c_str();
//	    Kmer next(cstr);
//	    kmers.push_back(next);
//	}
//
//	const size_t num_kmers = kmers.size();
//
//	unordered_map<size_t,vector<int>> arr;
//
//	//search each kmer in cdbg and return color set
//	int kmer_count = 0;
//	bool first = true;
//	bool wasEmpty = false;
//
//	for (const auto& kmer: kmers){
//
//		const const_UnitigColorMap<UnitigData> map = cdbg.find(kmer);
//
//		if (! map.isEmpty) {
//			const DataAccessor<UnitigData>* da = map.getData();
//			UnitigColors uc = da->getSubUnitigColors(map);
//			UnitigMap<UnitigData> newmap(0, 1, Kmer::k, map.strand);
//			for(UnitigColors::const_iterator it = uc.begin(newmap); it != uc.end(); ++it) {
//				const size_t color = it.getColorID();
//				//note to self: the iterator goes through all colors of the unitig, but we want to only keep the ones that the kmer is really annotated with
//				//if (uc.contains(map, color)){
//				std::unordered_map<size_t,vector<int>>::iterator iter = arr.find(color);
//
//				if (iter == arr.end()){
//					vector<int> vec(num_kmers, 0);
//					arr.insert({color,vec});
//				}
//				arr[color][kmer_count] = 1;
//						//}
//			}
//
//		}
//		kmer_count++;
//	}
//	return arr;
//}


unordered_map<size_t,vector<int>> GraphTraverser::search2(string query, int k){

	vector<Kmer> kmers;

	//split query into sequence of kmers (!!! query can contain a kmer multiple times, do not change the order of the kmers at this point!)
	for(int i = 0; i< (query.length()-k+1); ++i){
	    const string kmer = query.substr(i,k);
	    const char *cstr = kmer.c_str();
	    Kmer next(cstr);
	    kmers.push_back(next);
	}

	//test neighborhood function!
	string test = "ACCA";
	vector<string> neighborhood = compute_neighborhood(test, 1);
	cout << "Neighbohood of: " << test << endl;
	for (auto& elem: neighborhood){
		cout << elem << endl;
	}




	const size_t num_kmers = kmers.size();

	unordered_map<size_t,vector<int>> arr;

	//search each kmer in cdbg and return color set
	int kmer_count = 0;
	bool first = true;
	bool wasEmpty = false;
	UnitigColors* old_uc;

	for (const auto& kmer: kmers){
		UnitigMap<DataAccessor<UnitigData>, DataStorage<UnitigData>, false> map = cdbg.find(kmer);

		if (! map.isEmpty) {
			const DataAccessor<UnitigData>* da = map.getData();
			UnitigColors* uc = da->getUnitigColors(map);

			//ToDo: if this UnitigColors object contains the same colors as the object of the previous kmer (which is likely), then we already know whats happening!
			if (! first) {
				if (uc == old_uc && (! wasEmpty)){
					//we can simply copy the result of the previous kmer!
					for(auto& color : arr){
						arr[color.first][kmer_count] = arr[color.first][kmer_count -1];
					}
				} else {
					//ToDo: check these differences with Guillaume
					//for(UnitigColors::const_iterator it = uc->begin(map); it != uc->end(); it.nextColor()) {
					for (UnitigColors::const_iterator it = uc->begin(map); it != uc->end(); ++it) {

						const size_t color = it.getColorID();
						//note to self: the iterator goes through all colors of the unitig, but we want to only keep the ones that the kmer is really annotated with
						//if (uc -> contains(map, color)){
						if (it.getKmerPosition() == map.dist){
							std::unordered_map<size_t,vector<int>>::iterator iter = arr.find(color);

							if (iter == arr.end()){
								vector<int> vec(num_kmers, 0);
								arr.insert({color,vec});
							}
							arr[color][kmer_count] = 1;
						}
					}
				}
			} else {
				first = false;
				//for(UnitigColors::const_iterator it = uc->begin(map); it != uc->end(); it.nextColor()) {
				for (UnitigColors::const_iterator it = uc->begin(map); it != uc->end(); ++it) {
					const size_t color = it.getColorID();
					//note to self: the iterator goes through all colors of the unitig, but we want to only keep the ones that the kmer is really annotated with
					//if (uc -> contains(map, color)){
					if (it.getKmerPosition() == map.dist){
						std::unordered_map<size_t,vector<int>>::iterator iter = arr.find(color);

						if (iter == arr.end()){
							vector<int> vec(num_kmers, 0);
							arr.insert({color,vec});
						}
						arr[color][kmer_count] = 1;
					}
				}
			}
			old_uc = uc;
			wasEmpty = false;
		} else {
			wasEmpty = true;
		}
		kmer_count++;
	}
	return arr;
}



vector<string> GraphTraverser::compute_neighborhood(string kmer_str, int d){
	vector<char> alphabet = {'A','C'};
	vector<string> neighborhood;
	//string kmer_str = kmer.toString();
	vector<int> firstRow;
	for (int i = 0; i <= kmer_str.size(); ++i){
		firstRow.push_back(i);
	}

	GraphTraverser::searchNextRow("", kmer_str, firstRow, neighborhood, alphabet,d);

	return neighborhood;

}

void GraphTraverser::searchNextRow(string v, string& word, vector<int> lastRow, vector<string>& neighborhood, vector<char>& alphabet, int& d){
	int min = *(std::min_element(lastRow.begin(), lastRow.end()));
	if (min == d){
		int counter = word.size();
		for(auto& elem : lastRow){
			if (elem == d){
				//report v*w^x
				int pos = word.size() - counter;
				string suffix = word.substr(pos);
				neighborhood.push_back(v+suffix);
			}
			counter--;
		}
	} else if(min < d){
		for (char a : alphabet){
			vector<int> nextRow;
			//first entry can only be an insert
			nextRow.push_back(lastRow[0]+1);
			const char* str = word.c_str();
			for(int i = 1; i< lastRow.size(); ++i){
				int ins = lastRow[i]+1;
				int del = nextRow[i-1]+1;
				int sub;
				if (str[i-1] == a){
					sub = lastRow[i-1];
				} else {
					sub = lastRow[i-1]+1;
				}
				nextRow.push_back(std::min(std::min(ins, del), sub));
			}
			string next = v+a;
			GraphTraverser::searchNextRow(next,word,nextRow,neighborhood,alphabet,d);
		}
	}
}


void GraphTraverser::remove_singletonHits(unordered_map<size_t,vector<int>>& hits){

	const int k = cdbg.getK();

	vector<size_t> to_be_removed;

	for (auto& hit : hits){
		const vector<int> seq = hit.second;
		//if first and last appearance of '1' in vector is more than k appart, then everything is ok!
		int cnt_begin = 0;
		for (vector<int>::const_iterator i = seq.begin(); i != seq.end(); ++i){
			if (*i == 1){
				break;
			}
			cnt_begin++;
		}

		int cnt_end = 0;
		for (vector<int>::const_reverse_iterator i = seq.rbegin(); i != seq.rend(); ++i){
			if (*i == 1){
				break;
			}
			cnt_end++;
		}

		if ((seq.size() - cnt_end - cnt_begin) <= k){
			to_be_removed.push_back(hit.first);
		}
	}

	for (auto& elem : to_be_removed){
		hits.erase(elem);
	}

	//cout << "Number of singleton hits removed: " << to_be_removed.size() << endl;
}



//deprecated.
long double GraphTraverser::compute_p(long double& k, long double& x, long double& sigma){
	long double p = 1-pow(1-pow(sigma,-k),x);
	return p;
}


//deprecated.
unordered_map<size_t,double> GraphTraverser::compute_significance(unordered_map<size_t,vector<int>>& hits, long double& p) {
	unordered_map<size_t,double> p_values;

	for (auto& hit: hits){
		int q = hit.second.size();

		//compute number of 1's in vector
		int m = std::count(hit.second.begin(), hit.second.end(), 1);

		double r = gsl_cdf_binomial_Q(m, p, q);

		if (r > 0.05){
			cout << "p-value: " << r << endl;
		}
		p_values[hit.first] = r;
	}

	return p_values;
}


/*
 * _Approximate_ number of matches and mismatches from k-mer hits.
 * Atm, we use the average number of mismatches that can explain a run of 0's of a specific length.
 */
int GraphTraverser::compute_score(const vector<int>& hit){
	const int score_match = 1;
	const int score_mismatch = -2;

	const int l = hit.size();

	int mismatch = 0;
	int cnt = 0;
	for(auto& elem: hit){
		if (elem == 0){
			cnt+=1;
		} else if (cnt != 0){
			//cout << "cnt: " << cnt << endl;
			int local = ceil(float(cnt)/31);
			int local2 = floor(cnt - 31 + 1);
			int avg = (local+local2)/2;
			mismatch += avg;
			cnt = 0;
		}
	}
	if (cnt != 0){
		int local = ceil(float(cnt)/31);
		int local2 = floor(cnt - 31 + 1);
		int avg = (local+local2)/2;
		mismatch += avg;
	}

	int match = l - mismatch;
	if (match <= 0){
		cout << "ERROR! No matches!" << endl;
	}

	const int score = match*score_match+mismatch*score_mismatch;
	return score;
}




long double GraphTraverser::compute_evalue(const int& score, const double& db_size, const int& n) const{
	long double evalue = blast_k*db_size*n*exp(-lambda*score);
	return evalue;
}

long double GraphTraverser::compute_pvalue(const long double& evalue) const{
	long double pvalue = 1-exp(-evalue);
	return pvalue;
}


long double GraphTraverser::compute_log_evalue(const int& score, const double& db_size, const int& n) const{
	long double log_evalue = round(log(blast_k*db_size*n)-lambda*score);
	return log_evalue;

}

long double GraphTraverser::compute_log_pvalue(const long double& log_evalue) const{
	long double evalue = pow(10,log_evalue);
	if (1-exp(-evalue) > 0){
		return round(log(1-exp(-evalue)));
	} else {
		return round(log_evalue);
	}
}



void GraphTraverser::writePresenceMatrix(const unordered_map<size_t,vector<int>>& arr, const string& outfile, const unordered_map<size_t,long double>& pvalues){
	std::ofstream output(outfile,std::ofstream::binary);
	for (auto& elem : arr){
		const string color = cdbg.getColorName(elem.first);
		output << color << "\t" << pvalues.at(elem.first);
		for (auto& p : elem.second){
			output << "\t" << p;
		}
		output << endl;
	}
	output.close();
}





