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




unordered_map<size_t,vector<int>> GraphTraverser::search(string query, int k){

	vector<Kmer> kmers;

	//split query into sequence of kmers (!!! query can contain a kmer multiple times, do not change the order of the kmers at this point!)
	for(int i = 0; i< (query.length()-k+1); ++i){
	    string kmer = query.substr(i,k);
	    //cout << kmer << " " << kmer.length() << endl;
	    const char *cstr = kmer.c_str();
	    Kmer next(cstr);
	    kmers.push_back(next);
	}

	size_t num_kmers = kmers.size();

	unordered_map<size_t,vector<int>> arr;

	//search each kmer in cdbg and return color set
	int kmer_count = 0;
	for (auto& kmer: kmers){
		UnitigMap<DataAccessor<UnitigData>, DataStorage<UnitigData>, false> map = cdbg.find(kmer);
		//set<string> colors;
		if (map.isEmpty) {
			//cout << "kmer not found" << endl;
			//cout << kmer.toString() << endl;
		} else {
			DataAccessor<UnitigData>* da = map.getData();
			UnitigColors* uc = da->getUnitigColors(map);

			for(UnitigColors::const_iterator it = uc->begin(map); it != uc->end(); it.nextColor()) {
			  size_t color = it.getColorID();
				//note to self: the iterator goes through all colors of the unitig, but we want to only keep the ones that the kmer is really annotated with
				if (uc -> contains(map, color)){
					//colors.insert(cdbg.getColorName(color));
					std::unordered_map<size_t,vector<int>>::iterator iter = arr.find(color);

					if (iter == arr.end()){
						vector<int> vec(num_kmers, 0);
						arr.insert({color,vec});
					}

					arr[color][kmer_count] = 1;
				}
			}
		}
		kmer_count++;
	}

	return arr;
}



//void GraphTraverser::writeKmerPresence(vector<pair<Kmer,set<string>>> results, string& resfile) {
//
//	std::ofstream output;
//	output.open(resfile, std::ios_base::app);
//
//	for (auto& res : results){
//		set<string> colors = res.second;
//		for (auto& color : colors){
//			output << res.first.toString() << "\t" << color << endl;
//		}
//	}
//	output.close();
//}



void GraphTraverser::writePresenceMatrix(unordered_map<size_t,vector<int>>& arr, string& outfile){
	std::ofstream output(outfile,std::ofstream::binary);

	for (auto& elem : arr){
		string color = cdbg.getColorName(elem.first);
		output << color;
		for (auto& p : elem.second){
			output << "\t" << p;
		}
		output << endl;
	}
	output.close();
}





