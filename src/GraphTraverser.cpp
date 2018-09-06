/*
 * GraphTraverser.cpp
 *
 *  Created on: 26 Feb 2018
 *      Author: nina
 */

#include "GraphTraverser.hpp"

using namespace std;


GraphTraverser::GraphTraverser(ColoredCDBG<UnitigData>& graph) : cdbg(graph) {
	cout << "GraphTraverser initialized!" << endl;
}

/*
 * explore the bubble defined by the left and right gene by a depth-first search,
 * stopping each run if a max. threshold is met
 */
void GraphTraverser::exploreBubble(const string& left, const string& right, const int threshold){

	const char *left_str = left.c_str();
	const char *right_str = right.c_str();

	UnitigColorMap<UnitigData> first;
	UnitigColorMap<UnitigData> last;


	//for left, find last kmer that is found in the graph --- NOT NICE
	for(KmerIterator it_km(left_str), it_km_end; it_km != it_km_end; ++it_km){

		const UnitigColorMap<UnitigData> map = cdbg.find(it_km->first);

		if (!map.isEmpty) last = map;
	}


	//for right, find first kmer that is found in the graph
	for(KmerIterator it_km(right_str), it_km_end; it_km != it_km_end; ++it_km){

		const UnitigColorMap<UnitigData> map = cdbg.find(it_km->first);

		if (!map.isEmpty){

			first = map;
			break;
		}
	}

	//search ALL paths between the last and first unitig, with a search threshold for the total length of paths
	GraphTraverser::DFS_Iterative(last, first.getUnitigHead(), threshold);
}




//void GraphTraverser::BFS_Iterative(const UnitigColorMap<UnitigData>& start, Kmer& stop, int threshold){
//
//    queue<UnitigColorMap<UnitigData>> q; // Create queue of unitig to traverse
//    UnitigColorMap<UnitigData> ucm_tmp(start); // Create a non-const local copy of unitig given in parameter
//
//    DataAccessor<UnitigData>* da = ucm_tmp.getData(); // Get DataAccessor from unitig
//    UnitigData* data = da->getData(ucm_tmp); // Get boolean from DataAccessor
//
//    data->set_visited(); // Set boolean to indicate unitig was visited
//
//    q.push(ucm_tmp); // Push unitig to traverse on the stack
//
//    while (!q.empty()){ // While they are unitigs to traverse in the stack
//
//        ucm_tmp = q.front(); // Get unitig at the front of the queue
//
//        q.pop(); // Delete unitig at the front of the queue
//
//        bool stop_reached = false;
//
//        if (ucm_tmp.getUnitigHead() == stop){
//        	cout << "found stop!" << endl;
//
//
//        	//ToDo: every time I find a stop here, they are unconnected paths between the two anchoring nodes
//        	//hence for every stop, I then need to backtrace all possibilities to reach this point
//        }
//
//
//        for (auto& successor : ucm_tmp.getSuccessors()){ // Traverse successors
//
//            DataAccessor<UnitigData>* da_succ = successor.getData(); // Get DataAccessor from successor
//            UnitigData* data_succ = da_succ->getData(successor); // Get boolean from DataAccessor
//
//            if (data_succ->is_not_visited()){ // If boolean indicates the successor was not visited
//
//                data_succ->set_visited(); // Set boolean to indicate successor was visited
//
//                q.push(successor); // Traverse neighbors of successor
//            }
//        }
//    }
//}




void GraphTraverser::DFS_Iterative(const UnitigColorMap<UnitigData>& start, const Kmer& stop, const int threshold){
	//record the path we are iterating over, and push it to the result if we reach the stop Kmer

	stack<UnitigColorMap<UnitigData>> stck; // Create stack of unitig to traverse
	UnitigColorMap<UnitigData> ucm_tmp(start); // Create a non-const local copy of unitig given in parameter

	stck.push(ucm_tmp); // Push first unitig to traverse on the stack


	bool stop_reached = false;
	vector<UnitigColorMap<UnitigData>> currentPath;

	while (!stck.empty()){ // While there are unitigs to traverse in the stack
	    	ucm_tmp = stck.top(); // Get the unitig on top of the stack

	    	stck.pop(); // Delete unitig on the top of the stack

	    	currentPath.push_back(ucm_tmp);


	    	//Case1: we found the end
	    	if (ucm_tmp.getUnitigHead() == stop){
	    	    //scout << "Found the end!" << endl;
	    	    stop_reached = true;

	    	    //ToDo: clear marking
	    	    for (auto& unitig : currentPath){
	    	    	DataAccessor<UnitigData>* pathda = unitig.getData(); // Get DataAccessor from unitig
	    	    	UnitigData* pathdata = pathda->getData(unitig); // Get boolean from DataAccessor
	    	    	pathdata->set_on_target_path();
	    	    }
	    	}


	    	//Case2: we reached the limit
	    	if (currentPath.size() > threshold){
	    		//cout << "Threshold reached!" << endl;
	    		stop_reached = true;

//	    		for (auto& unitig : currentPath){
//	    			DataAccessor<UnitigData>* pathda = unitig.getData(); // Get DataAccessor from unitig
//	    			UnitigData* pathdata = pathda->getData(unitig); // Get boolean from DataAccessor
//	    			pathdata->set_seen();
//	    		}
	    	}


	    	if (stop_reached){
	    		//backtrace
	    		stop_reached = false;

	    		//adjust current path such that only the predecessor of the next unitig on the stack is present
	    		UnitigColorMap<UnitigData> ucm_next = stck.top();
	    		DataAccessor<UnitigData>* nextda = ucm_next.getData(); // Get DataAccessor from unitig
	    		UnitigData* nextdata = nextda->getData(ucm_next); // Get boolean from DataAccessor

	    		for (auto& predecessor : ucm_next.getPredecessors()){
	    			for (vector<UnitigColorMap<UnitigData>>::iterator it = currentPath.begin() ; it != currentPath.end(); ++it){
	    				if (predecessor == *it){
	    					currentPath.erase(it+1,currentPath.end());
	    					break;
	    				}
	    			}
	    		}


	    	} else {
	    		//Case 3: just go on
	    		DataAccessor<UnitigData>* da = ucm_tmp.getData(); // Get DataAccessor from unitig
	    		UnitigData* data = da->getData(ucm_tmp); // Get boolean from DataAccessor

	    		data->set_visited(); // Set boolean to indicate unitig was visited


	    		for (auto& successor : ucm_tmp.getSuccessors()) {
	    			DataAccessor<UnitigData>* succda = successor.getData(); // Get DataAccessor from unitig
	    			UnitigData* succdata = succda->getData(successor); // Get boolean from DataAccessor
	    			if(succdata->is_not_visited() || succdata->is_on_target_path()){
	    				stck.push(successor);
	    			}
	    		}
	    	}
	}


}





//void GraphTraverser::DFS_Iterative(const UnitigColorMap<UnitigData>& start, Kmer& stop, int threshold){
//
//
//	//record the path we are iterating over, and push it to the result if we reach the stop Kmer
//
//    stack<UnitigColorMap<UnitigData>> stck; // Create stack of unitig to traverse
//    UnitigColorMap<UnitigData> ucm_tmp(start); // Create a non-const local copy of unitig given in parameter
//
//    stck.push(ucm_tmp); // Push first unitig to traverse on the stack
//
//
//    bool stop_reached = false;
//
//    vector<UnitigColorMap<UnitigData>> currentPath;
//
//    while (!stck.empty()){ // While there are unitigs to traverse in the stack
//
//    	cout << currentPath.size() << endl;
//
//    	ucm_tmp = stck.top(); // Get the unitig on top of the stack
//
//    	stck.pop(); // Delete unitig on the top of the stack
//
//    	//we might have peviously found a stop criterion
//    	if (stop_reached) {
//    		cout << "clear path" << endl;
//    		stop_reached = false;
//
//
//    		vector<UnitigColorMap<UnitigData>> oldPath = currentPath;
//    		currentPath.clear();
//    		//cout << currentPath.size() << endl;
//
//
//    		//Backtrace to last branching point
//    		for (auto& predecessor : ucm_tmp.getPredecessors()){
//    			DataAccessor<UnitigData>* pre_da = predecessor.getData(); // Get DataAccessor from unitig
//    			UnitigData* pre_data = pre_da-> getData(predecessor); // Get boolean from DataAccessor
//    			if (pre_data->is_visited()){
//    				//we found the branching node, re-add old path to current path up to this point!
//    				for (auto& unitig : oldPath){
//    					if (unitig != predecessor){
//    						currentPath.push_back(unitig);
//    					}
//    				}
//    			}
//    		}
//
//    		//cout << "adjusted currentPath" << endl;
//    	}
//
//    	if (currentPath.size() > threshold){
//    		cout << "threshold reached!" << endl;
//    		stop_reached = true;
//    		//cout << currentPath[currentPath.size()-1].toString() << endl;
//    	}
//
//    	if (ucm_tmp.getUnitigHead() == stop){
//    		cout << "Found the end!" << endl;
//    		stop_reached = true;
//
//    		//mark current path found as target path!
//    		for (auto& elem : currentPath){
//    			DataAccessor<UnitigData>* da = elem.getData(); // Get DataAccessor from unitig
//    			UnitigData* data = da-> getData(elem); // Get boolean from DataAccessor
//    			data-> set_on_target_path();
//    		}
//
//    	} else if (! stop_reached) {
//
//    		DataAccessor<UnitigData>* da = ucm_tmp.getData(); // Get DataAccessor from unitig
//    		UnitigData* data = da->getData(ucm_tmp); // Get boolean from DataAccessor
//
//    		if (data->is_not_visited() || data -> is_on_target_path()){ // If boolean indicates the unitig was not visited
//
//    			data->set_visited(); // Set boolean to indicate unitig was visited
//    			currentPath.push_back(ucm_tmp);
//
//
//    			for (auto& successor : ucm_tmp.getSuccessors()) {
//    				stck.push(successor);
//    			}
//
//
//    		} else {
//    			cout << "already visited" << endl;
//    			//we found a dead end, backtrack current path until last merging point -> backtrack to the last time we added more than 1 successor to the stack!
//    			stop_reached = true;
//    		}
//    	}
//    }
//}






void GraphTraverser::exploreSubgraph(const string& s) const {
	//ToDo: check graph structure of left and right border sequences! If we naively only look at last and first k-mer of the borders, we will miss everything with sequencing errors/mutations in there!

	Kmer old;

	for(KmerIterator it_km(s.c_str()), it_km_end; it_km != it_km_end; ++it_km){

		const const_UnitigColorMap<UnitigData> map = cdbg.find(it_km->first);

		cout << "Kmer: " << it_km->first.toString() << endl;

		if (! map.isEmpty) {

			bool new_unitig = false;

			//record if the previous kmer was present on the same unitig
			const Kmer head(map.getUnitigHead());

			if (head != old){

				cout << "new unitig!" << endl;
				cout << map.referenceUnitigToString() << endl;
				cout << map.size << endl;

				old = head;
				new_unitig = true;
			}

			if (new_unitig){
				//count number of predecessors and successors of this unitig
				int succ_count = 0;
				int pre_count = 0;

				for (const auto& successor : map.getSuccessors()) ++succ_count;
				for (const auto& predecessor : map.getPredecessors()) ++pre_count;

				cout << "Current unitig #successors: " << succ_count << endl;
				cout << "Current unitig #predecessors: " << pre_count << endl;
			}

		}
		else cout << "empty map!" << endl;
	}
}

unordered_map<size_t,vector<int>> GraphTraverser::search(const string& query, const int k, const int ndistance) const {

	const size_t num_kmers = query.length() - k + 1;

	//split query into sequence of kmers (!!! query can contain a kmer multiple times, do not change the order of the kmers at this point!)
	const char *query_str = query.c_str();

	//test neighborhood function!
	//string test = "ACA";
	//vector<Kmer> neighborhood = compute_neighborhood(test, 2);
	//cout << "Neighbohood of: " << test << endl;
	//for (auto& elem: neighborhood){
	//	cout << elem.toString() << endl;
	//}
	//cout << neighborhood.size() << endl;

	unordered_map<size_t,vector<int>> arr;


	//search each kmer in cdbg and return color set
	int kmer_count = 0;
	bool first = true;
	bool wasEmpty = false;

	const UnitigColors* old_uc;

	for(KmerIterator it_km(query_str), it_km_end; it_km != it_km_end; ++it_km){

		const const_UnitigColorMap<UnitigData> map = cdbg.find(it_km->first);

		if (!map.isEmpty) {

			const DataAccessor<UnitigData>* da = map.getData();
			const UnitigColors* uc = da->getUnitigColors(map);

			const bool copy = (!first && !wasEmpty && (uc == old_uc));

			if (copy) {

				//ToDo: if this UnitigColors object contains the same colors as the object of the previous kmer (which is likely), then we already know whats happening!
				for(auto& color : arr){

					if (arr[color.first][kmer_count - 1] == 1) arr[color.first][kmer_count] = 1;
				}
			}
			else {

				first = false;

				for (UnitigColors::const_iterator it = uc->begin(map); it != uc->end(); it.nextColor()) {
				//for (UnitigColors::const_iterator it = uc->begin(map); it != uc->end(); ++it) {
						const size_t color = it.getColorID();

						if (uc->contains(map, color)){ //note to self: the iterator goes through all colors of the unitig, but we want to only keep the ones that the kmer is really annotated with

							const std::unordered_map<size_t, vector<int>>::const_iterator iter = arr.find(color);

							if (iter == arr.end()) arr.insert({color, vector<int>(num_kmers, 0)});

								arr[color][kmer_count] = 1;
							}
						}

				}



			//ToDo: refactor!

			//now check the k-mers neighborhood too, but do not overwrite perfect matches!
			if (ndistance > 0){

				const vector<Kmer> neighborhood = GraphTraverser::compute_neighborhood(it_km->first.toString(), ndistance);

				cout << "n-size: " << neighborhood.size() << endl;

				for (const auto& nkmer : neighborhood){

					const const_UnitigColorMap<UnitigData> nmap = cdbg.find(nkmer);

					if (!nmap.isEmpty) {

						const DataAccessor<UnitigData>* nda = nmap.getData();
						const UnitigColors* nuc = nda->getUnitigColors(nmap);

						for (UnitigColors::const_iterator nit = nuc->begin(nmap); nit != nuc->end(); ++nit) {

							const size_t ncolor = nit.getColorID();

							if (nuc -> contains(nmap, ncolor)){

								std::unordered_map<size_t,vector<int>>::iterator niter = arr.find(ncolor);

								if (niter == arr.end()){

									arr.insert({ncolor, vector<int>(num_kmers, 0)});
									arr[ncolor][kmer_count] = 2;

								}
								else if (arr[ncolor][kmer_count] == 0) arr[ncolor][kmer_count] = 2;
							}
						}
					}
				}
			}

			old_uc = uc;
		}
		wasEmpty = map.isEmpty;

		++kmer_count;
	}

	return arr;
}




vector<Kmer> GraphTraverser::compute_neighborhood(const string& kmer_str, const int d) const {

	static const size_t alphabet_sz = 4;
	static const char alphabet[alphabet_sz] = {'A','C','G','T'};

	vector<Kmer> neighborhood;
	//string kmer_str = kmer.toString();
	vector<int> firstRow;

	firstRow.reserve(kmer_str.size());

	for (int i = 0; i <= kmer_str.size(); ++i) firstRow.push_back(i);

	GraphTraverser::searchNextRow(string(""), kmer_str, firstRow, neighborhood, alphabet, alphabet_sz, d);

	return neighborhood;

}

void GraphTraverser::searchNextRow(const string& v, const string& word, const vector<int>& lastRow, vector<Kmer>& neighborhood, const char* alphabet, const size_t alphabet_sz, const int d) const {

	int min = *(std::min_element(lastRow.begin(), lastRow.end()));

	if ((min <= d) && (min > 0) && (v.length() == word.length())) {

		const string suffix(v.substr(0, word.length()));

		neighborhood.push_back(Kmer(suffix.c_str()));

	} else if (min == d){

		//int counter = word.size();
		for(const auto& elem : lastRow){
			if (elem == d){
				//report v*w^x
				if (v.length() == word.length()){
					//cout << "v: " << v << endl;
					neighborhood.push_back(Kmer(v.c_str()));
					break;
				}
				else if (v.length() < word.length()) {

					//cout << pos << endl;
					const string suffix = word.substr(v.length());
					const string concat = v + suffix;
					//cout << "concat: " << concat << endl;
					neighborhood.push_back(Kmer(concat.c_str()));
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
	}
	else if (min < d){

		for (size_t j = 0; j != alphabet_sz; ++j){

			const char* str = word.c_str();

			vector<int> nextRow;
			//first entry can only be an insert
			nextRow.push_back(lastRow[0]+1);

			for (int i = 1; i< lastRow.size(); ++i){

				const int ins = lastRow[i]+1;
				const int del = nextRow[i-1]+1;
				const int sub = lastRow[i-1] + (str[i-1] != alphabet[j]);

				nextRow.push_back(std::min(std::min(ins, del), sub));
			}

			const string next = v + alphabet[j];

			GraphTraverser::searchNextRow(next, word, nextRow, neighborhood, alphabet,alphabet_sz, d);
		}
	}
}


void GraphTraverser::remove_singletonHits(unordered_map<size_t, vector<int>>& hits) const {

	const int k = cdbg.getK();

	vector<size_t> to_be_removed;

	for (const auto& hit : hits){

		const vector<int>& seq = hit.second;
		//if first and last appearance of '1' in vector is more than k appart, then everything is ok!
		int cnt_begin = 0;
		int cnt_end = 0;

		for (vector<int>::const_iterator i = seq.begin(); i != seq.end(); ++i){

			if (*i == 1) break;
			++cnt_begin;
		}

		for (vector<int>::const_reverse_iterator i = seq.rbegin(); i != seq.rend(); ++i){

			if (*i == 1) break;
			++cnt_end;
		}

		if ((seq.size() - cnt_end - cnt_begin) <= k) to_be_removed.push_back(hit.first);
	}

	for (const auto elem : to_be_removed) hits.erase(elem);

	//cout << "Number of singleton hits removed: " << to_be_removed.size() << endl;
}



//deprecated.
unordered_map<size_t,double> GraphTraverser::compute_significance(const unordered_map<size_t,vector<int>>& hits, const long double p) const {
	unordered_map<size_t,double> p_values;

	for (const auto& hit: hits){

		const int q = hit.second.size();

		//compute number of 1's in vector
		const int m = std::count(hit.second.begin(), hit.second.end(), 1);

		const double r = gsl_cdf_binomial_Q(m, p, q);

		if (r > 0.05) cout << "p-value: " << r << endl;

		p_values[hit.first] = r;
	}

	return p_values;
}


/*
 * _Approximate_ number of matches and mismatches from k-mer hits.
 * Atm, we use the average number of mismatches that can explain a run of 0's of a specific length.
 */
int GraphTraverser::compute_score(const vector<int>& hit) const {

	const int score_match = 1;
	const int score_mismatch = -2;
	const int l = hit.size();

	int mismatch = 0;
	int cnt = 0;

	for (const auto elem : hit){

		if (elem == 0) ++cnt;
		else if (cnt != 0){
			//cout << "cnt: " << cnt << endl;
			const int local = ceil(float(cnt)/31);
			const int local2 = floor(cnt - 31 + 1);
			const int avg = (local+local2)/2;

			mismatch += avg;
			cnt = 0;
		}
	}

	if (cnt != 0){
		const int local = ceil(float(cnt)/31);
		const int local2 = floor(cnt - 31 + 1);
		const int avg = (local+local2)/2;
		mismatch += avg;
	}

	const int match = l - mismatch;
	if (match <= 0) cout << "ERROR! No matches!" << endl;


	return match*score_match+mismatch*score_mismatch;

}



//long double GraphTraverser::compute_evalue(const int& score, const double& db_size, const int& n) const{
//	long double evalue = blast_k*db_size*n*exp(-lambda*score);
//	return evalue;
//}
//
//long double GraphTraverser::compute_pvalue(const long double& evalue) const{
//	long double pvalue = 1-exp(-evalue);
//	return pvalue;
//}
//
//
//long double GraphTraverser::compute_log_evalue(const int& score, const double& db_size, const int& n) const{
//	long double log_evalue = round(log(blast_k*db_size*n)-lambda*score);
//	return log_evalue;
//
//}
//
//long double GraphTraverser::compute_log_pvalue(const long double& log_evalue) const{
//	long double evalue = pow(10,log_evalue);
//	if (1-exp(-evalue) > 0){
//		return round(log(1-exp(-evalue)));
//	} else {
//		return round(log_evalue);
//	}
//}



void GraphTraverser::writePresenceMatrix(const unordered_map<size_t,vector<int>>& arr, const string& outfile, const unordered_map<size_t,long double>& pvalues) const {

	std::ofstream output(outfile, std::ofstream::binary);

	for (const auto& elem : arr){

		const string color = cdbg.getColorName(elem.first);

		output << color << "\t" << pvalues.at(elem.first);

		for (const auto& p : elem.second) output << "\t" << p;

		output << endl;
	}

	output.close();
}


/*
 * Find the unitig that corresponds to this string, and report all colors of this unitig
 */

vector<string> GraphTraverser::getColors(const string& u) const {

	vector<string> colors;

	const int k = cdbg.getK();
	const string kmer = u.substr(0,k);
	const char *cstr = kmer.c_str();

	const Kmer head(cstr);

	const const_UnitigColorMap<UnitigData> map = cdbg.find(head, true);

	if (!map.isEmpty) {

		const DataAccessor<UnitigData>* da = map.getData();
		const UnitigColors* uc = da->getUnitigColors(map);

		for (UnitigColors::const_iterator it = uc->begin(map); it != uc->end(); it.nextColor()) {

			colors.push_back(cdbg.getColorName(it.getColorID()));
		}

	}
	else cout << "Error, unitig not found." << endl;

	return colors;
}





