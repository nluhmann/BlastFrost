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




vector<pair<Kmer,set<string>>> GraphTraverser::search(string query, int k){

	vector<Kmer> kmers;

	//split query into sequence of kmers (!!! query can contain a kmer multiple times, do not change the order of the kmers at this point!)
	for(int i = 0; i< (query.length()-k+1); ++i){
	    string kmer = query.substr(i,k);
	    //cout << kmer << " " << kmer.length() << endl;
	    const char *cstr = kmer.c_str();
	    Kmer next(cstr);
	    kmers.push_back(next);
	}

	vector<pair<Kmer,set<string>>> presence;

	//search each kmer in cdbg and return color set
	for (auto& kmer: kmers){
		UnitigMap<DataAccessor<UnitigData>, DataStorage<UnitigData>, false> map = cdbg.find(kmer);
		if (map.isEmpty) {
			cout << "kmer not found" << endl;
			cout << kmer.toString() << endl;
		} else {
			DataAccessor<UnitigData>* da = map.getData();
			UnitigColors* uc = da->getUnitigColors(map);

			set<string> colors;
			for(UnitigColors::const_iterator it = uc->begin(); it != uc->end(); it.nextColor(map.len)) {
				size_t color = it->getColorID(map.len);

				//note to self: the iterator goes through all colors of the unitig, but we want to only keep
				if (uc -> contains(map, color)){
					colors.insert(cdbg.getColorName(color));
				}
			}

			presence.push_back(std::make_pair(kmer,colors));
		}
	}

	return presence;
}



void GraphTraverser::writeKmerPresence(vector<pair<Kmer,set<string>>> results, string& resfile) {

	std::ofstream output;
	output.open(resfile, std::ios_base::app);

	for (auto& res : results){
		set<string> colors = res.second;
		for (auto& color : colors){
			output << res.first.toString() << "\t" << color << endl;
		}

	}
	output.close();
}


