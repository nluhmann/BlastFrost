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

unordered_map<size_t,vector<std::string>> SubGraphTraverser::extractSubGraph(const string& query, const int k, const int distance) {

	//ToDo: this should be a vector of Kmer!
	QuerySearch::searchResultSubgraph res = que.search_kmers(query, k, distance, db_size);

	cout << "Found seed kmers" << endl;
	cout << "Extending..." << res.mapping.size() << endl;
	unordered_map<size_t,vector<Kmer>> map = res.mapping;



	unordered_map<size_t,vector<std::string>> all_sequences;

	//int counter = 0;
	for(const auto& color : map){


		//add this node to the collected path
		vector<Kmer> path;
		vector<Kmer> current = color.second;


		std::reverse(current.begin(),current.end());
		int add_count = 0;

		while (! current.empty()){

			if(add_count > 50){

				//remove the last Kmer in current, and remove the last 200 Kmers in the path, then try again
				if (current.size() > 1){
					current.pop_back();
					path.erase(path.end()-50,path.end());
					add_count = 0;
				} else {
					path.clear();
					cout << "clear path" << endl;
					break;
				}

			}


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

							}
						}
					}
				}

				//if (sameCol.size() == 2 && sameCol[0].toString().compare(current.back().toString()) != 0 && sameCol[1].toString().compare(current.back().toString()) != 0){
				if (sameCol.size() == 1 && sameCol[0].toString().compare(current.back().toString()) != 0){

					current.push_back(sameCol[0]);
					++add_count;
//					if (debug) {
//						cout << "extend" << endl;
//					}


				} else {
					add_count = 0;
				}
			}
		}

		int diff_prefix = 0;
		if (res.prefix_missing[color.first] > res.prefix_offset[color.first]){
			//extend search
		} else {
			//cut prefix unitig respectively
			diff_prefix = (res.prefix_offset[color.first] - res.prefix_missing[color.first]);
		}

		int diff_suffix = 0;
		if (res.suffix_missing[color.first] > res.suffix_offset[color.first]){
			//extend search
		} else {
			//cut suffix unitig respectively
			diff_suffix = (res.suffix_offset[color.first] - res.suffix_missing[color.first]);
		}

		if (path.size() > 0){
			vector<std::string> sequences = pathSequence(path, diff_prefix, diff_suffix);

			all_sequences[color.first] = sequences;

		}
	}



	//SubGraphTraverser::pathLength(all_paths, query.size());
	//SubGraphTraverser::testPath(all_paths);



	return all_sequences;
}




SubGraphTraverser::subgraphs SubGraphTraverser::extractSubGraph_intelligent(const string& query, const int k, const int distance) {

	QuerySearch::searchResultSubgraph res = que.search_kmers(query, k, distance, db_size);

	//cout << "Extending..." << res.mapping.size() << endl;
	unordered_map<size_t,vector<Kmer>> map = res.mapping;

	subgraphs result;

	unordered_map<size_t, vector<size_t>> groups = SubGraphTraverser::groupSeedHits(map);

	unordered_map<size_t,vector<std::string>> all_sequences;

	for (auto& group: groups){

		vector<size_t> colors = group.second;
		colors.push_back(group.first);

		while (! colors.empty()){
			size_t current_color = colors.back();
			colors.pop_back();

			//copy vector to keep track of colors
			vector<size_t> color_group = colors;

			//run subgraph extraction for current color, while recording all colors in current group that do not differ!
			vector<Kmer> path;
			vector<Kmer> current = map[current_color];

			std::reverse(current.begin(),current.end());
			int add_count = 0;

			while (! current.empty()){

				if(add_count > 50){
					//remove the last Kmer in current, and remove the last 200 Kmers in the path, then try again
					if (current.size() > 1){
						current.pop_back();
						path.erase(path.end()-50,path.end());
						add_count = 0;
					} else {
						path.clear();
						break;
					}

				}

				Kmer first = current.back();
				current.pop_back();
				path.push_back(first);

				if (! current.empty()){
					vector<Kmer> sameCol;
					UnitigColorMap<UnitigData> ucm = cdbg.find(first);

					for (const auto& successor : ucm.getSuccessors()){
						Kmer head = successor.getMappedHead();

						if (head == current.back()){
							break;
						} else if (! (head == first)){

							const DataAccessor<UnitigData>* da = successor.getData();
							const UnitigColors* uc = da->getUnitigColors(successor);
							for (UnitigColors::const_iterator it = uc->begin(successor); it != uc->end(); it.nextColor()) {
								if (it.getColorID() == current_color && sameCol.empty()){
									sameCol.push_back(head);
									break;
								} else if ((! sameCol.empty()) && it.getColorID() == current_color){
									sameCol.clear();
									break;
								}
							}
						}
					}

					if (sameCol.size() == 1 && sameCol[0].toString().compare(current.back().toString()) != 0){
						current.push_back(sameCol[0]);
						++add_count;

						if (color_group.size() > 1){
							UnitigColorMap<UnitigData> ucm = cdbg.find(sameCol[0]);
							const DataAccessor<UnitigData>* da = ucm.getData();
							const UnitigColors* uc = da->getUnitigColors(ucm);
							for (auto& col : color_group){
								if (! uc->contains(ucm,col)) {
									color_group.erase(std::remove(color_group.begin(), color_group.end(), col), color_group.end());
								}
							}
						}

					} else {
						add_count = 0;
					}
				}
			}

			if (color_group.size() < colors.size()){
				//remove all colors that have been found with this path!
				for (auto& col: color_group){
					colors.erase(std::remove(colors.begin(), colors.end(), col), colors.end());
				}
			} else {
				colors.clear();
			}


			int diff_prefix = 0;
			if (res.prefix_missing[current_color] > res.prefix_offset[current_color]){
				//extend search
			} else {
				//cut prefix unitig respectively
				diff_prefix = (res.prefix_offset[current_color] - res.prefix_missing[current_color]);
			}

			int diff_suffix = 0;
			if (res.suffix_missing[current_color] > res.suffix_offset[current_color]){
				//extend search
			} else {
				//cut suffix unitig respectively
				diff_suffix = (res.suffix_offset[current_color] - res.suffix_missing[current_color]);
			}

			if (path.size() > 0){
				vector<std::string> sequences = pathSequence(path, diff_prefix, diff_suffix);

				result.colors[current_color] = color_group;
				result.sequences[current_color] = sequences;
			}

		}
	}

	return result;

}

unordered_map<size_t, vector<size_t>> SubGraphTraverser::groupSeedHits(unordered_map<size_t,vector<Kmer>>& map){

	unordered_map<size_t,vector<Kmer>> already_seen;

	unordered_map<size_t,vector<size_t>> groups;

	for (auto& new_elem : map){
		bool found = false;
		for(auto& seen_elem : already_seen){
			if(seen_elem.second == new_elem.second){
				groups[seen_elem.first].push_back(new_elem.first);
				found = true;
				break;
			}
		}
		if (! found){
			already_seen[new_elem.first] = new_elem.second;
			vector<size_t> newgroup;
			groups[new_elem.first] = newgroup;
		}
	}


	return groups;
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

				//push results with current number to file
				path_sequences.push_back(result);
				result.clear();
				result += next;


			} else {
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
				if (! (previous_suffix.compare(next_prefix) == 0)){
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
			//int i = 0;
			for(std::vector<int>::size_type i = 0; i != p.size(); i++) {
				output << ">" << queryID << "|" << color_name << "_" << i << "_len_" << p[i].length() << endl;
				output << p[i] << endl;
			}
		}
	}

	output.close();
}

void SubGraphTraverser::printPaths_intelligent(string& outprefix, string& query, SubGraphTraverser::subgraphs result, const string& queryID) {
//output paths in fasta format, encode color name and paths length in header

	std::ofstream output(outprefix+"_"+query+"_subgraph.fasta",std::ios_base::app);

	for (const auto& color : result.colors){
		string color_name = cdbg.getColorName(color.first);
		vector<std::string> p = result.sequences[color.first];
		if (p.size() == 1) {
			output << ">" << queryID << "|" << color_name << "_len_" << p[0].length() << endl;
			output << p[0] << endl;
		} else {
			//int i = 0;
			for(std::vector<int>::size_type i = 0; i != p.size(); i++) {
				output << ">" << queryID << "|" << color_name << "_" << i << "_len_" << p[i].length() << endl;
				output << p[i] << endl;
			}
		}

		vector<size_t> add_colors = color.second;
		for (auto& col : add_colors) {
			string color_name = cdbg.getColorName(col);
			if (p.size() == 1) {
				output << ">" << queryID << "|" << color_name << "_len_" << p[0].length() << endl;
				output << p[0] << endl;
			} else {
				//int i = 0;
				for(std::vector<int>::size_type i = 0; i != p.size(); i++) {
					output << ">" << queryID << "|" << color_name << "_" << i << "_len_" << p[i].length() << endl;
					output << p[i] << endl;
				}
			}
		}

	}

	output.close();
}


void SubGraphTraverser::pathLength(const unordered_map<size_t,vector<Kmer>>& all_paths) {
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













