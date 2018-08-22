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





/*
 * explore the bubble defined by the left and right gene by a depth-first search,
 * stopping each run if a max. threshold is met
 */
void GraphTraverser::exploreBubble(string left, string right, int threshold){

	//for left, find last kmer that is found in the graph --- NOT NICE
	vector<Kmer> kmers;
	int k = cdbg.getK();
	//split query into sequence of kmers
	for(int i = 0; i< (left.length()-k+1); ++i){
		const string kmer = left.substr(i,k);
		const char *cstr = kmer.c_str();
		Kmer next(cstr);
		kmers.push_back(next);
	}

	UnitigColorMap<UnitigData> last;
	for (const auto& kmer: kmers){
		UnitigColorMap<UnitigData> map = cdbg.find(kmer);
		if (! map.isEmpty) {
			last = map;
		}
	}


	//for right, find first kmer that is found in the graph
	kmers.clear();
	//split query into sequence of kmers
	for(int i = 0; i< (right.length()-k+1); ++i){
		const string kmer = right.substr(i,k);
		const char *cstr = kmer.c_str();
		Kmer next(cstr);
		kmers.push_back(next);
	}

	UnitigColorMap<UnitigData> first;
	for (const auto& kmer: kmers){
		UnitigColorMap<UnitigData> map = cdbg.find(kmer);
		if (! map.isEmpty) {
				first = map;
				break;
		}
	}

	//search paths between the last and first unitig, with a search threshold for the total length of paths
	Kmer stop = first.getUnitigHead();

	GraphTraverser::DFS_Iterative(last, stop, threshold);


}


void GraphTraverser::DFS_Iterative(const UnitigColorMap<UnitigData>& start, Kmer& stop, int threshold){


	//record the path we are iterating over, and push it to the result if we reach the stop Kmer

    stack<UnitigColorMap<UnitigData>> stck; // Create stack of unitig to traverse
    UnitigColorMap<UnitigData> ucm_tmp(start); // Create a non-const local copy of unitig given in parameter

    stck.push(ucm_tmp); // Push first unitig to traverse on the stack


    bool stop_reached = false;

    vector<UnitigColorMap<UnitigData>> currentPath;

    while (!stck.empty()){ // While there are unitigs to traverse in the stack
    	cout << currentPath.size() << endl;
    	ucm_tmp = stck.top(); // Get the unitig on top of the stack

    	stck.pop(); // Delete unitig on the top of the stack

    	//we might have peviously found a stop criterion
    	if (stop_reached) {
    		cout << "clear path" << endl;
    		stop_reached = false;


    		vector<UnitigColorMap<UnitigData>> oldPath = currentPath;
    		currentPath.clear();
    		cout << currentPath.size() << endl;

    		for (auto& predecessor : ucm_tmp.getPredecessors()){
    			DataAccessor<UnitigData>* pre_da = predecessor.getData(); // Get DataAccessor from unitig
    			UnitigData* pre_data = pre_da-> getData(predecessor); // Get boolean from DataAccessor
    			if (pre_data->is_visited()){
    				//we found the branching node, re-add old path to current path up to this point!
    				for (auto& unitig : oldPath){
    					if (unitig != predecessor){
    						currentPath.push_back(unitig);
    					}

    				}
    			}
    		}

    		cout << "adjusted currentPath" << endl;
    	}

    	if (currentPath.size() > threshold){
    		cout << "threshold reached!" << endl;
    		stop_reached = true;
    		cout << currentPath[currentPath.size()-1].toString() << endl;
    	}

    	if (ucm_tmp.getUnitigHead() == stop){
    		cout << "Found the end!" << endl;
    		stop_reached = true;
    		//ToDo: mark path found as target path!

    		for (auto& elem : currentPath){
    			DataAccessor<UnitigData>* da = elem.getData(); // Get DataAccessor from unitig
    			UnitigData* data = da-> getData(elem); // Get boolean from DataAccessor
    			data-> set_on_target_path();
    		}


    	} else {

    		DataAccessor<UnitigData>* da = ucm_tmp.getData(); // Get DataAccessor from unitig
    		UnitigData* data = da->getData(ucm_tmp); // Get boolean from DataAccessor

    		//ToDo: if it has already been visited, but is on a target path, we need to go on again! Q: can this result in an endless loop?
    		//HowTo: backtrack when a path is found!
    		if (data->is_not_visited() || data -> is_on_target_path()){ // If boolean indicates the unitig was not visited

    			data->set_visited(); // Set boolean to indicate unitig was visited
    			currentPath.push_back(ucm_tmp);


    			for (auto& successor : ucm_tmp.getSuccessors()) {
    				stck.push(successor);
    			}


    		} else {
    			//we found a dead end, backtrack current path until last merging point -> backtrack to the last time we added more than 1 successor to the stack!
    			stop_reached = true;
    		}
    	}
    }
}












void GraphTraverser::exploreSubgraph(string s){
	//ToDo: check graph structure of left and right border sequences! If we naively only look at last and first k-mer of the borders, we will miss everything with sequencing errors/mutations in there!

	vector<Kmer> kmers;
	int k = cdbg.getK();
	//split query into sequence of kmers
	for(int i = 0; i< (s.length()-k+1); ++i){
		const string kmer = s.substr(i,k);
		const char *cstr = kmer.c_str();
		Kmer next(cstr);
		kmers.push_back(next);
	}

	Kmer old;
	for (const auto& kmer: kmers){
		UnitigColorMap<UnitigData> map = cdbg.find(kmer);

		cout << "Kmer: " << kmer.toString() << endl;
		if (! map.isEmpty) {
			bool new_unitig = false;


			//record if the previous kmer was present on the same unitig
			Kmer head = map.getUnitigHead();
			if (head == old){
			} else {
				cout << "new unitig!" << endl;
				cout << map.toString() << endl;
				cout << map.size << endl;
				old = head;
				new_unitig = true;
			}

			if (new_unitig){
				//count number of predecessors and successors of this unitig
				int succ_count = 0;
				for (auto& successor : map.getSuccessors()){
					succ_count += 1;
				}

				cout << "Current unitig #successors: " << succ_count << endl;

				int pre_count = 0;
				for (auto& predecessor : map.getPredecessors()){
					pre_count += 1;
				}

				cout << "Current unitig #predecessors: " << pre_count << endl;
			}

		}
		else{
			cout << "empty map!" << endl;
		}

	}

}






unordered_map<size_t,vector<int>> GraphTraverser::search(string query, int k, int ndistance){

	vector<Kmer> kmers;

	//split query into sequence of kmers (!!! query can contain a kmer multiple times, do not change the order of the kmers at this point!)
	for(int i = 0; i< (query.length()-k+1); ++i){
	    const string kmer = query.substr(i,k);
	    const char *cstr = kmer.c_str();
	    Kmer next(cstr);
	    kmers.push_back(next);
	}

	//test neighborhood function!
	//string test = "ACA";
	//vector<Kmer> neighborhood = compute_neighborhood(test, 2);
	//cout << "Neighbohood of: " << test << endl;
	//for (auto& elem: neighborhood){
	//	cout << elem.toString() << endl;
	//}
	//cout << neighborhood.size() << endl;

	const size_t num_kmers = kmers.size();

	unordered_map<size_t,vector<int>> arr;

	//search each kmer in cdbg and return color set
	int kmer_count = 0;
	bool first = true;
	bool wasEmpty = false;
	UnitigColors* old_uc;

	for (const auto& kmer: kmers){
		UnitigColorMap<UnitigData> map = cdbg.find(kmer);

		if (! map.isEmpty) {
			const DataAccessor<UnitigData>* da = map.getData();
			UnitigColors* uc = da->getUnitigColors(map);
			bool copy = false;

			//ToDo: if this UnitigColors object contains the same colors as the object of the previous kmer (which is likely), then we already know whats happening!
			if (! first) {
				if (uc == old_uc && (! wasEmpty)){
					//we can simply copy the result of the previous kmer!
					copy = true;
					for(auto& color : arr){
						if (arr[color.first][kmer_count -1] == 1){
							arr[color.first][kmer_count] = 1;
						}

					}
				}
			}

			if (! copy) {
				first = false;

				for(UnitigColors::const_iterator it = uc->begin(map); it != uc->end(); it.nextColor()) {
				//for (UnitigColors::const_iterator it = uc->begin(map); it != uc->end(); ++it) {
					const size_t color = it.getColorID();

					if (uc -> contains(map, color)){ //note to self: the iterator goes through all colors of the unitig, but we want to only keep the ones that the kmer is really annotated with

						std::unordered_map<size_t,vector<int>>::iterator iter = arr.find(color);

						if (iter == arr.end()){
							vector<int> vec(num_kmers, 0);
							arr.insert({color,vec});
						}
						arr[color][kmer_count] = 1;
					}
				}
			}

			//ToDo: refactor!

			//now check the k-mers neighborhood too, but do not overwrite perfect matches!
			if(ndistance > 0){
				vector<Kmer> neighborhood = GraphTraverser::compute_neighborhood(kmer.toString(), ndistance);
				cout << "n-size: " << neighborhood.size() << endl;
				for (auto& nkmer : neighborhood){
					UnitigMap<DataAccessor<UnitigData>, DataStorage<UnitigData>, false> nmap = cdbg.find(nkmer);
					if (! nmap.isEmpty) {
						const DataAccessor<UnitigData>* nda = nmap.getData();
						UnitigColors* nuc = nda->getUnitigColors(nmap);

						for (UnitigColors::const_iterator nit = nuc->begin(nmap); nit != nuc->end(); ++nit) {
							const size_t ncolor = nit.getColorID();
							if (nuc -> contains(nmap, ncolor)){
								std::unordered_map<size_t,vector<int>>::iterator niter = arr.find(ncolor);

								if (niter == arr.end()){
									vector<int> vec(num_kmers, 0);
									arr.insert({ncolor,vec});
									arr[ncolor][kmer_count] = 2;
								} else if (arr[ncolor][kmer_count] == 0){
									arr[ncolor][kmer_count] = 2;
								}
							}
						}
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



vector<Kmer> GraphTraverser::compute_neighborhood(string kmer_str, int d){
	vector<char> alphabet = {'A','C','G','T'};
	vector<Kmer> neighborhood;
	//string kmer_str = kmer.toString();
	vector<int> firstRow;
	for (int i = 0; i <= kmer_str.size(); ++i){
		firstRow.push_back(i);
	}

	GraphTraverser::searchNextRow("", kmer_str, firstRow, neighborhood, alphabet,d);

	return neighborhood;

}

void GraphTraverser::searchNextRow(string v, string& word, vector<int> lastRow, vector<Kmer>& neighborhood, vector<char>& alphabet, int& d){

	int min = *(std::min_element(lastRow.begin(), lastRow.end()));
	if(min <= d && min > 0 && v.length() == word.length()){
		string suffix = v.substr(0,word.length());
		//cout << "suffix: " << suffix << endl;
		Kmer next_kmer(suffix.c_str());
		neighborhood.push_back(next_kmer);
	} else if (min == d){
		//int counter = word.size();
		for(auto& elem : lastRow){
			if (elem == d){
				//report v*w^x
				if (v.length() == word.length()){
					//cout << "v: " << v << endl;
					Kmer next_kmer(v.c_str());
					neighborhood.push_back(next_kmer);
					break;
				}
				else if (v.length() < word.length()) {
					int pos = v.length();
					//cout << pos << endl;
					string suffix = word.substr(pos);
					string concat = v+suffix;
					//cout << "concat: " << concat << endl;
					Kmer next_kmer(concat.c_str());
					neighborhood.push_back(next_kmer);
					break;
				}




				//cout << "concat: " << concat << endl;
				//cout << "----" << endl;

				//we can only search for kmers in neighborhood of same length
				//if(concat.length() == word.length()){
				//	Kmer next_kmer(concat.c_str());
				//	neighborhood.push_back(next_kmer);
				//}
			}
			//counter--;
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


/*
 * Find the unitig that corresponds to this string, and report all colors of this unitig
 */
vector<string> GraphTraverser::getColors(const string& u){
	vector<string> colors;
	const int k = cdbg.getK();
	const string kmer = u.substr(0,k);
	const char *cstr = kmer.c_str();
	Kmer head(cstr);

	UnitigMap<DataAccessor<UnitigData>, DataStorage<UnitigData>, false> map = cdbg.find(head,true);

	if (! map.isEmpty) {
		const DataAccessor<UnitigData>* da = map.getData();
		UnitigColors* uc = da->getUnitigColors(map);
		for (UnitigColors::const_iterator it = uc->begin(map); it != uc->end(); it.nextColor()) {
			const size_t colorID = it.getColorID();
			string color = cdbg.getColorName(colorID);
			colors.push_back(color);
		}

	} else {
		cout << "Error, unitig not found." << endl;
	}

	return colors;
}





