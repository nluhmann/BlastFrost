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


	//ToDo: cannot do that. first, run BlastFrost for left and right string with d>2,



	//for left, find last kmer that is found in the graph --- NOT NICE
	for(KmerIterator it_km(left_str), it_km_end; it_km != it_km_end; ++it_km){

		const UnitigColorMap<UnitigData> map = cdbg.find(it_km->first);

		if (!map.isEmpty){
			last = map;
		}
	}


	//for right, find first kmer that is found in the graph
	for(KmerIterator it_km(right_str), it_km_end; it_km != it_km_end; ++it_km){

		const UnitigColorMap<UnitigData> map = cdbg.find(it_km->first);

		if (!map.isEmpty){

			first = map;
			break;
		}
	}


	last.dist = 0;
	last.len = last.size - last.getGraph() -> getK() + 1;

	//first.dist = 0;
	//first.len = first.size - first.getGraph() -> getK() + 1;

	//search ALL paths between the last and first unitig, with a search threshold for the total length of paths
	vector<vector<UnitigColorMap<UnitigData>>> allPaths = GraphTraverser::DFS_Iterative(last, first.getUnitigHead(), threshold);

	GraphTraverser::printPaths(allPaths);


}


void GraphTraverser::printPaths(const vector<vector<UnitigColorMap<UnitigData>>> allPaths){
	int pathcounter = 0;


	std::ofstream p("paths",std::ofstream::binary);
	for (auto& path : allPaths){
		int unitigCounter = 0;

		for (auto& elem : path) {
			p << ">path" << pathcounter << "_" << unitigCounter << endl;

			p << elem.mappedSequenceToString() << endl;

			const DataAccessor<UnitigData>* da = elem.getData();
			const UnitigColors* uc = da->getUnitigColors(elem);

			p << "#";
			for (UnitigColors::const_iterator it = uc->begin(elem); it != uc->end(); it.nextColor()) {
				const size_t colorID = it.getColorID();
				const string color = cdbg.getColorName(colorID);
				p << color << ",";
			}
			p << endl;
			++unitigCounter;
		}

		++pathcounter;
	}

	p.close();

}




vector<vector<UnitigColorMap<UnitigData>>> GraphTraverser::DFS_Iterative(const UnitigColorMap<UnitigData>& start, const Kmer& stop, const int threshold){
	//record the path we are iterating over, and push it to the result if we reach the stop Kmer

	stack<UnitigColorMap<UnitigData>> stck; // Create stack of unitig to traverse
	UnitigColorMap<UnitigData> ucm_tmp(start); // Create a non-const local copy of unitig given in parameter

	stck.push(ucm_tmp); // Push first unitig to traverse on the stack

	vector<vector<UnitigColorMap<UnitigData>>> allPaths;
	bool stop_reached = false;
	vector<UnitigColorMap<UnitigData>> currentPath;

	while (!stck.empty()){ // While there are unitigs to traverse in the stack

	    	ucm_tmp = stck.top(); // Get the unitig on top of the stack

	    	stck.pop(); // Delete unitig on the top of the stack

	    	if(std::find(currentPath.begin(), currentPath.end(), ucm_tmp) != currentPath.end()) {
	    	    stop_reached = true;
	    	} else {
	    		currentPath.push_back(ucm_tmp);
	    	}



	    	DataAccessor<UnitigData>* da = ucm_tmp.getData(); // Get DataAccessor from unitig
	    	UnitigData* data = da->getData(ucm_tmp); // Get boolean from DataAccessor

    		data->set_visited(); // Set boolean to indicate unitig was visited

	    	//Case1: we found the end
	    	if (ucm_tmp.getUnitigHead() == stop){
	    	    cout << "Found the end!" << endl;
	    	    stop_reached = true;
	    	    allPaths.push_back(currentPath);
	    	    for (auto& unitig : currentPath){
	    	    	DataAccessor<UnitigData>* pathda = unitig.getData(); // Get DataAccessor from unitig
	    	    	UnitigData* pathdata = pathda->getData(unitig); // Get boolean from DataAccessor
	    	    	pathdata->set_on_target_path();
	    	    }
	    	}


	    	//Case2: we reached the limit
	    	if (currentPath.size() > threshold){
	    		stop_reached = true;
	    	}


	    	if (stop_reached){
	    		//backtrace
	    		stop_reached = false;

	    		//adjust current path such that only the predecessor of the next unitig on the stack is present
	    		UnitigColorMap<UnitigData> ucm_next = stck.top();

	    		DataAccessor<UnitigData>* nextda = ucm_next.getData(); // Get DataAccessor from unitig
	    		UnitigData* nextdata = nextda->getData(ucm_next); // Get boolean from DataAccessor

	    		bool broken = false;
	    		for (auto& predecessor : ucm_next.getPredecessors()){
	    			for (vector<UnitigColorMap<UnitigData>>::reverse_iterator it = currentPath.rbegin(); it != currentPath.rend(); ++it ) {
	    				if (predecessor == *(it.base())){
	    					//std::advance(it, 1);
	    					currentPath.erase((--it).base(),currentPath.end());
	    					broken = true;
	    					break;
	    				}
	    			}
	    		}
	    		if (! broken){
	    			currentPath.clear();
	    		}




	    	} else {
	    		//Case 3: just go on

	    		bool succs = false;
	    		for (auto& successor : ucm_tmp.getSuccessors()) {
	    			DataAccessor<UnitigData>* succda = successor.getData(); // Get DataAccessor from unitig
	    			UnitigData* succdata = succda->getData(successor); // Get boolean from DataAccessor
	    			if(succdata->is_not_visited()){
	    				stck.push(successor);
	    				succs = true;
	    			} else if (succdata -> is_on_target_path()) {
	    				//we can stop here, cause the rest of the path has been found already previously!
	    				allPaths.push_back(currentPath);
	    			}

	    		}
	    		if (! succs){
	    			if (! stck.empty()){
	    				UnitigColorMap<UnitigData> ucm_next = stck.top();

	    				DataAccessor<UnitigData>* nextda = ucm_next.getData(); // Get DataAccessor from unitig
	    				UnitigData* nextdata = nextda->getData(ucm_next); // Get boolean from DataAccessor

	    				bool broken = false;
	    				for (auto& predecessor : ucm_next.getPredecessors()){
	    					for (vector<UnitigColorMap<UnitigData>>::reverse_iterator it = currentPath.rbegin(); it != currentPath.rend(); ++it ) {
	    						if (predecessor == *(it.base())){
	    							currentPath.erase((--it).base(),currentPath.end());
	    							broken = true;
	    							break;
	    						}
	    					}
	    				}
	    				if (! broken){
	    					currentPath.clear();
	    				}


	    			}
	    		}


	    	}
	}
	return allPaths;
}



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





void GraphTraverser::extractSubGraph(const string& query, const int k, const int distance) {

	//ToDo: this should be a vector of Kmer!
	unordered_map<size_t,vector<Kmer>> map;

	//1) search for all k-mer hits inside of neighborhood defined by distance, save references to all unitigs covered
	//ToDo: we will only consider exact matches for now!
	const char *query_str = query.c_str();
	int startcounter = 0;

	for(KmerIterator it_km(query_str), it_km_end; it_km != it_km_end; ++it_km){
		startcounter += 1;
		const UnitigColorMap<UnitigData> ucm = cdbg.find(it_km->first);

		if (!ucm.isEmpty) {

			const DataAccessor<UnitigData>* da = ucm.getData();
			const UnitigColors* uc = da->getUnitigColors(ucm);

			//for each color, store all UnitigColorMap references as a set
			//ToDo: we can likely speed this up!
			for (UnitigColors::const_iterator it = uc->begin(ucm); it != uc->end(); it.nextColor()) {
				const size_t color = it.getColorID();

				//check if color already in map
				//ToDo: do we need to do that?
				const std::unordered_map<size_t, vector<Kmer>>::iterator iter = map.find(color);
				Kmer head = ucm.getUnitigHead();

				if (iter == map.end()) {
					cout << color << endl;
					cout << startcounter << endl;

					vector<Kmer> newset;
					map.insert({color, newset});
					map[color].push_back(head);
				} else if (map[color].back() != head){
					map[color].push_back(head);
				}
			}
		}
	}


	//2) for each color, give me a path!
	//ToDo: the path will be the same for some colors, or even subpaths, so they should be discovered simultanously!
//	for(const auto& color : map){
//		//first, give me the numbr of unitigs for each color...
//		cout << "color: " << color.first << endl;
//		cout << "#unitigs: " << color.second.size() << endl;
//
//		//for each UnitigMap in a color, check if any successor is also in the list!
//		for (auto& map : color.second){
//			cout << "seq: " << map.referenceUnitigToString() << endl;
//			for (const auto& successor : map.getSuccessors()){
//				bool colorFound = false;
//
//				const DataAccessor<UnitigData>* da = successor.getData();
//				const UnitigColors* uc = da->getUnitigColors(successor);
//				for (UnitigColors::const_iterator it = uc->begin(successor); it != uc->end(); it.nextColor()) {
//					if (it.getColorID() == color.first){
//						colorFound = true;
//					}
//				}
//
//				if (colorFound){
//					cout << "color found" << endl;
//				}
//
//				if(std::find(color.second.begin(), color.second.end(), successor) != color.second.end()) {
//				    cout << "map in list" << endl;
//				}
//
//			}
//		}
//	}

		//Strategy: for each color, starting from the first unitig in the list, find all successors of the same color
		//if there is only one successor of the same color, add it to the path and pop it from the seed list if possible
		//otherwise, follow each successor path until a unitig from the seed list is found. Ignore the other path(s).
		//stop if the last unitig from the seed list is found

		//ToDo: Afterwards, for each color, we might want to extent the path at the end or the beginning...but how far?
		//We can first extend it following the longest color path? How can I anchor that?


	unordered_map<size_t,vector<Kmer>> all_paths;
	cout << "complete paths" << endl;
	for(const auto& color : map){
		//add this node to the collected path
		vector<Kmer> path;

		vector<Kmer> current = color.second;
		std::reverse(current.begin(),current.end());

		cout << color.first << endl;
		while (! current.empty()){
			//cout << current.size() << endl;
			Kmer first = current.back();
			current.pop_back();

			path.push_back(first);

			//find all successors with the same color
			vector<Kmer> sameCol;

			const UnitigColorMap<UnitigData> ucm = cdbg.find(first);

			for (const auto& successor : ucm.getSuccessors()){
				const DataAccessor<UnitigData>* da = successor.getData();
				const UnitigColors* uc = da->getUnitigColors(successor);

				for (UnitigColors::const_iterator it = uc->begin(successor); it != uc->end(); it.nextColor()) {
					if (it.getColorID() == color.first){
						sameCol.push_back(successor.getUnitigHead());
					}
				}
			}

			if (sameCol.empty()){
				//there is no successor of the same color, this path stops here
				cout << "no successor!" << endl;

				//Note: the while expression will still try to extend other seed unitigs for this color!




			} else if (sameCol.size() == 1){
				if (! current.empty() && sameCol.back() == current.back()){
					//the single successor is already in the seed list
					cout << "seed already in list" << endl;
				} else if (! current.empty()){
					//the successor is not in the list, can fill a gap between seed unitigs
					current.push_back(sameCol.back());
				}
			} else {
				cout << "more than 2!" << endl;
			}
		}
		all_paths[color.first] = path;
	}


	//different strategy: add seed unitigs to a new cdbg
	//add missing unitigs completing paths to new cdbg
	//in the end, the new cdbg should contain one CC
	//write new cdbg to gfa file!

	GraphTraverser::pathLength(all_paths, query.size());

	//take the length of the query as reference
	//previously, for each color, record the length of the prefix of the reference before the first seed hit
	//(record the length of the reference after the last seed hit as well?!)
	//for each color, extend the path of seeds in both ends until the length distance to the reference is minimized



	//at this point, assume we filled the prefix of the paths up already
	//so we can fill up its suffix until we reach the approximate length of the reference

	for(const auto& color : all_paths){





	}





	//to think about: output these unitigs as simple fasta for each color, then rebuild Bifrost graph with smaller k to get gfa file -> look at it in bandage





}



void GraphTraverser::pathLength(const unordered_map<size_t,vector<Kmer>>& all_paths, const int& ref_length) {
	//compute the length of each of the paths
		for (const auto& color : all_paths){
			vector<Kmer> current = color.second;
			int length = 0;
			for(auto& head : current){
				int unitig_length = cdbg.find(head).referenceUnitigToString().size();
				if (length == 0){
					length = unitig_length;
				} else {
					length += unitig_length - 31;
				}
			}

			cout << color.first << endl;
			cout << length << endl;
		}
}




unordered_map<size_t,vector<int>> GraphTraverser::search(const string& query, const int k, const int ndistance, const string& query_name) const {

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


		//const const_UnitigColorMap<UnitigData> map = cdbg.find(it_km->first);
		UnitigColorMap<UnitigData> map = cdbg.find(it_km->first);

		//It.nextColor() iterates over the first k-mer of each color in the mapping
		//++it iterate over all pair (kmer pos, color) in the mapping

		//For each ++it, you check if the k-mer position is the previous k-mer position + 1 or not
		// Then you check the color

		map.dist = 0;
		map.len = map.size - Kmer::k + 1;
		map.strand = true;

		if (!map.isEmpty) {

			const DataAccessor<UnitigData>* da = map.getData();
			const UnitigColors* uc = da->getUnitigColors(map);


			//const UnitigColors sub_uc = da->getSubUnitigColors(map);

			bool copy = (!first && !wasEmpty && (uc == old_uc));

			old_uc = uc;

			if (copy) {

				//ToDo: if this UnitigColors object contains the same colors as the object of the previous kmer (which is likely), then we already know whats happening!
				for(auto& color : arr){

					if (arr[color.first][kmer_count - 1] == 1) arr[color.first][kmer_count] = 1;
				}

			}
			else {

				first = false;

				for (UnitigColors::const_iterator it = uc->begin(map); it != uc->end(); it.nextColor()) {

						const size_t color = it.getColorID();

						//if (uc->contains(map, color)){ //note to self: the iterator goes through all colors of the unitig, but we want to only keep the ones that the kmer is really annotated with

							const std::unordered_map<size_t, vector<int>>::const_iterator iter = arr.find(color);

							if (iter == arr.end()) arr.insert({color, vector<int>(num_kmers, 0)});

								arr[color][kmer_count] = 1;
						//}

				}


			}



			//ToDo: refactor!

			//now check the k-mers neighborhood too, but do not overwrite perfect matches!
			if (ndistance > 0){

				const vector<Kmer> neighborhood = GraphTraverser::compute_neighborhood(it_km->first.toString(), ndistance);

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
//unordered_map<size_t,double> GraphTraverser::compute_significance(const unordered_map<size_t,vector<int>>& hits, const long double p) const {
//	unordered_map<size_t,double> p_values;
//
//	for (const auto& hit: hits){
//
//		const int q = hit.second.size();
//
//		//compute number of 1's in vector
//		const int m = std::count(hit.second.begin(), hit.second.end(), 1);
//
//		const double r = gsl_cdf_binomial_Q(m, p, q);
//
//		if (r > 0.05) cout << "p-value: " << r << endl;
//
//		p_values[hit.first] = r;
//	}
//
//	return p_values;
//}


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



