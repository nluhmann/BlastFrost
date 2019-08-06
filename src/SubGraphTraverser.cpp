/*
 * GraphTraverser.cpp
 *
 *  Created on: 26 Feb 2018
 *      Author: nina
 */

#include "SubGraphTraverser.hpp"

using namespace std;


SubGraphTraverser::SubGraphTraverser(ColoredCDBG<UnitigData>& graph, double& db, QuerySearch& q) : cdbg(graph), db_size(db), que(q){
	cout << "GraphTraverser initialized!" << endl;
	//QuerySearch que(graph);

}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// FUNC 2: Explore subgraph that contains queried sequence, i.e. all variants of a gene queried against the graph //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

unordered_map<size_t,vector<std::string>> SubGraphTraverser::extractSubGraph(const string& query, const int k, const int distance, string& outprefix, string& queryfile) {

	//ToDo: this should be a vector of Kmer!
	QuerySearch::searchResultSubgraph res = que.search_kmers(query, k, distance, db_size);
	unordered_map<size_t,vector<Kmer>> map = res.mapping;

	//if prefix_missing > prefix_offset: we need to extend our search at the prefix
	//if prefix_missing < prefix_offset: we need to cut the first unitig by (prefix_offset - prefix_missing) positions
	//if suffix_missing > suffix_offset: we need to extend our search at the suffix
	//if suffix_missing < suffix_offset: we need to cut the last unitig by (suffix_offset - suffix_missing) positions


//	//1) search for all k-mer hits inside of neighborhood defined by distance, save references to all unitigs covered
//	//ToDo: we will only consider exact matches for now!
//	const char *query_str = query.c_str();
//	unordered_map<size_t,vector<Kmer>> map;
//
//
//	for(KmerIterator it_km(query_str), it_km_end; it_km != it_km_end; ++it_km){
//		UnitigColorMap<UnitigData> ucm = cdbg.find(it_km->first);
//
//
//		ucm.dist = 0;
//		ucm.len = ucm.size - Kmer::k + 1;
//		//ucm.strand = true;
//
//
//		if (!ucm.isEmpty) {
//
//			const DataAccessor<UnitigData>* da = ucm.getData();
//			const UnitigColors* uc = da->getUnitigColors(ucm);
//
//
//			for (UnitigColors::const_iterator it = uc->begin(ucm); it != uc->end(); it.nextColor()) {
//				const size_t color = it.getColorID();
//
//				//check if color already in map
//				//ToDo: do we need to do that?
//				const std::unordered_map<size_t, vector<Kmer>>::iterator iter = map.find(color);
//				//Kmer head = ucm.getUnitigHead();
//				Kmer head = ucm.getMappedHead();
//
//				if (iter == map.end()) {
//					//cout << color << endl;
//					//cout << startcounter << endl;
//
//					vector<Kmer> newset;
//					map.insert({color, newset});
//					map[color].push_back(head);
//
//				} else if (map[color].back() != head){ //check if this unitig has been found already
//					map[color].push_back(head);
//				}
//			}
//		}
//	}
//
//
//	cout << "OLD MAP SIZE: " << map.size() << endl;
//	cout << "NEW MAP SIZE: " << map2.size() << endl;

	//for each color, output last unitig
//	for(const auto& color : map){
//		vector<Kmer> current = color.second;
//		cout << "color: " << cdbg.getColorName(color.first) << endl;
//
//		cout << cdbg.find(current.back()).referenceUnitigToString() << endl;
//		cout << cdbg.find(current.front()).referenceUnitigToString() << endl;
//	}



		//Strategy: for each color, starting from the first unitig in the list, find all successors of the same color
		//if there is only one successor of the same color, add it to the path and pop it from the seed list if possible
		//otherwise, follow each successor path until a unitig from the seed list is found. Ignore the other path(s).
		//stop if the last unitig from the seed list is found

		//ToDo: Afterwards, for each color, we might want to extent the path at the end or the beginning...but how far?
		//We can first extend it following the longest color path? How can I anchor that?


	unordered_map<size_t,vector<std::string>> all_sequences;

	for(const auto& color : map){


		//add this node to the collected path
		vector<Kmer> path;
		vector<Kmer> current = color.second;


//		cout << cdbg.getColorName(color.first) << endl;
//		cout << "CURR: " << current.size() << endl;
//		bool debug = false;
//		if (current.size() <= 2){
//			debug = true;
//		}
//
//		if (debug){
//			for (auto& elem: current){
//				cout << elem.toString() << endl;
//			}
//		}

		std::reverse(current.begin(),current.end());
		int add_count = 0;

		while (! current.empty()){

			if(add_count > 200){
				path.clear();
//				if (cdbg.getColorName(color.first) == "assemblies/SAL_LA1103AA_AS.scaffold.fasta"){
//					cout << "clear path" << endl;
//				}
				cout << "clear path" << endl;
				break;
			}

//			if (cdbg.getColorName(color.first) == "assemblies/SAL_LA1103AA_AS.scaffold.fasta"){
//				cout << add_count << endl;
//			}

			Kmer first = current.back();
			current.pop_back();

			path.push_back(first);



			if (! current.empty()){
				//find all successors with the same color
				vector<Kmer> sameCol;

				UnitigColorMap<UnitigData> ucm = cdbg.find(first);

				for (const auto& successor : ucm.getSuccessors()){

					const DataAccessor<UnitigData>* da = successor.getData();
					const UnitigColors* uc = da->getUnitigColors(successor);
					//Kmer head = successor.getUnitigHead();


					Kmer head = successor.getMappedHead();
					//Kmer alternative = successor.getUnitigHead();

					if (! (head == first)){

						for (UnitigColors::const_iterator it = uc->begin(successor); it != uc->end(); it.nextColor()) {
							if (it.getColorID() == color.first){
								sameCol.push_back(head);
								//sameCol.push_back(alternative);
//								if (cdbg.getColorName(color.first) == "assemblies/SAL_EA7814AA_AS.scaffold.fasta"){
//										cout << head.toString() << endl;
//										cout << alternative.toString() << endl;
//								}
							}
						}
					}
				}

				//if (sameCol.size() == 2 && sameCol[0].toString().compare(current.back().toString()) != 0 && sameCol[1].toString().compare(current.back().toString()) != 0){
				if (sameCol.size() == 1 && sameCol[0].toString().compare(current.back().toString()) != 0){

//					if (cdbg.getColorName(color.first) == "assemblies/SAL_LA1103AA_AS.scaffold.fasta"){
//						cout << sameCol[0].toString() << endl;
//						cout << current.back().toString() << endl;
//					}

					current.push_back(sameCol[0]);
					++add_count;

				} else {

					add_count = 0;
				}

				//			if (sameCol.empty()){
				//				//there is no successor of the same color, this path stops here
				//				cout << "no successor!" << endl;
				//
				//				//Note: the while expression will still try to extend other seed unitigs for this color!
				//
				//
				//			} else if (sameCol.size() == 1){
				//
				////				if ((!current.empty()) && sameCol.back().toString() == current.back().toString()){
				////					//the single successor is already in the seed list
				////				} else if ((!current.empty()) && ( std::find(path.begin(), path.end(), sameCol.back()) == path.end())){
				////					//the successor is not in the list, can fill a gap between seed unitigs
				////					current.push_back(sameCol.back());
				////				}
				//				if ((!current.empty()) && ( std::find(path.begin(), path.end(), sameCol.back()) == path.end())){
				//					cout << "push back" << endl;
				//					//the successor is not in the list, can fill a gap between seed unitigs
				//					current.push_back(sameCol.back());
				////					cout << "Not yet in List" << endl;
				////					cout << sameCol.back().toString() << endl;
				////					cout << path.back().toString() << endl;
				//				} else {
				//					cout << "ALREADY IN LIST" << endl;
				//				}
				//			} else {
				//				cout << "Repeat" << endl;
				//				//We have a repeat, I repeat, we have a repeat!
				//				//cout << "more than 2!" << endl;
				//				//cout << sameCol.size() << endl;
				//			}
			}
		}

//		cout << res.prefix_missing[color.first] << endl;
//		cout << res.prefix_offset[color.first] << endl;
//		cout << res.suffix_missing[color.first] << endl;
//		cout << res.suffix_offset[color.first] << endl;

		int diff_prefix = 0;
		if (res.prefix_missing[color.first] > res.prefix_offset[color.first]){
			//extend search
			//ToDo how to extent to more untigs? orientation?
			//cout << "EXTEND" << endl;
		} else {
			//cut prefix unitig respectively
			diff_prefix = (res.prefix_offset[color.first] - res.prefix_missing[color.first]);
		}

		int diff_suffix = 0;
		if (res.suffix_missing[color.first] > res.suffix_offset[color.first]){
			//extend search
			//ToDo how to extent to more untigs? orientation?
			//cout << "EXTEND" << endl;
		} else {
			//cut suffix unitig respectively
			diff_suffix = (res.suffix_offset[color.first] - res.suffix_missing[color.first]);
		}

		if (path.size() > 0){
			vector<std::string> sequences = pathSequence(path, diff_prefix, diff_suffix);

			all_sequences[color.first] = sequences;

		} else {
			//cout << cdbg.getColorName(color.first) << endl;
		}
	}

//	cout << "all colors done" << endl;


	//SubGraphTraverser::pathLength(all_paths, query.size());
	//SubGraphTraverser::testPath(all_paths);



	return all_sequences;
}





vector<std::string> SubGraphTraverser::pathSequence(const vector<Kmer>& path, const int& diff_prefix, const int& diff_suffix) {
	vector<std::string> path_sequences;
	string result;
	string previous_suffix;
	string previous_seq;

	bool first = true;
	for(auto& head : path){
		UnitigColorMap<UnitigData> ucm = cdbg.find(head);
		ucm.dist = 0;
		ucm.len = ucm.size - Kmer::k + 1;

		string next = ucm.mappedSequenceToString();
		//string next = ucm.referenceUnitigToString();
		string next_prefix = next.substr(0,30);

		if (! previous_suffix.empty()){
			if (previous_suffix.compare(next_prefix) != 0){
				//cout << "ERROR path not continuous, suffix: " << previous_suffix << endl;
				//cout << "ERROR path not continuous, prefix: " << next_prefix << endl;

//				cout << "BREAK" << endl;
//				cout << next_prefix << endl;
//				cout << previous_suffix << endl;
//				cout << result << endl;
//				cout << next << endl;


				//push results with current number to file
				path_sequences.push_back(result);
				result.clear();
				result += next;


			} else {
				//cout << "EXTEND" << endl;
				string suff = next.substr(30);
				result += suff;
			}
		} else {
			if (first && diff_prefix > 0){
				result += next.substr(diff_prefix-1);
				first = false;
			} else {
				result += next;
			}
		}
		previous_seq = next;
		previous_suffix = next.substr(next.size() - 30);

	}

	string result_end;
	if (diff_suffix > 0){
		result_end = result.substr(0,result.length()-diff_suffix-1);
	} else {
		result_end = result;
	}

	path_sequences.push_back(result_end);
	//cout << all_paths_sequences.size() << endl;


	return path_sequences;
}





////////////////////////////////////////////////////////////////////////////////
// HELPER FUNCTIONS
////////////////////////////////////////////////////////////////////////////////


void SubGraphTraverser::testPath(const unordered_map<size_t,vector<Kmer>>& all_paths) {

	//for each path, check if it continuous in the graph, i.e. the unitig overlaps are correct!
	for (const auto& color : all_paths) {
		vector<Kmer> current = color.second;
		string previous_suffix;
		string previous_seq;
		for (auto& head : current){

			UnitigColorMap<UnitigData> ucm = cdbg.find(head);

			ucm.dist = 0;
			ucm.len = ucm.size - Kmer::k + 1;
			//ucm.strand = true;

			string next = ucm.mappedSequenceToString();

			string next_prefix = next.substr(0,30);
			if (! previous_suffix.empty()){
				if (! previous_suffix.compare(next_prefix) == 0){
					cout << "ERROR path not continuous, suffix: " << previous_suffix << endl;
					cout << "ERROR path not continuous, prefix: " << next_prefix << endl;
				}
			}
			previous_seq = next;
			previous_suffix = next.substr(next.size() - 30);
		}
	}
}


void SubGraphTraverser::printPaths(string& outprefix, string& query, unordered_map<size_t,vector<std::string>>& paths, const string& queryID) {
//output paths in fasta format, encode color name and paths length in header

	std::ofstream output(outprefix+"_"+query+"_subgraph.fasta",std::ios_base::app);

	for (const auto& color : paths){
		string color_name = cdbg.getColorName(color.first);
		vector<std::string> p = color.second;
		if (p.size() == 1) {
			output << ">" << queryID << "|" << color_name << "_len_" << color.second[0].length() << endl;
			output << color.second[0] << endl;
		} else {
			int i = 0;
			for(std::vector<int>::size_type i = 0; i != p.size(); i++) {
				output << ">" << queryID << "|" << color_name << "_" << i << "_len_" << p[i].length() << endl;
				output << p[i] << endl;
			}
		}
	}

	output.close();
}




void SubGraphTraverser::pathLength(const unordered_map<size_t,vector<Kmer>>& all_paths, const int& ref_length) {
	//compute the length of each of the paths
		for (const auto& color : all_paths){
			vector<Kmer> current = color.second;
			int length = 0;
			cout << "color:" << color.first << endl;
			for(auto& head : current){
				int unitig_length = cdbg.find(head).referenceUnitigToString().size();

				if (length == 0){
					length = unitig_length;
				} else {
					length += unitig_length - 30;
				}
			}
			cout << length << endl;
		}
}





void SubGraphTraverser::exploreSubgraph(const string& s) const {
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


/*
 * Find the unitig that corresponds to this string, and report all colors of this unitig
 */
vector<string> SubGraphTraverser::getColors(const string& u) const {

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













