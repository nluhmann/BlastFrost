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




void GraphTraverser::search(string query, int k){

	vector<Kmer> kmers;

	//split query into sequence of kmers (!!! query can contain a kmer multiple times, do not change the order of the kmers at this point!)
	for(int i = 0; i< (query.length()-k+1); ++i){
	    string kmer = query.substr(i,k);
	    //cout << kmer << " " << kmer.length() << endl;
	    const char *cstr = kmer.c_str();
	    Kmer next(cstr);
	    kmers.push_back(next);
	}


	//search each kmer in cdbg and return color set
	for (kmer: kmers){
		UnitigMap<DataAccessor<UnitigData>, DataStorage<UnitigData>, false> map = cdbg.find(kmer);
		if (map.isEmpty) {
			cout << "kmer not found" << endl;
			cout << kmer.toString() << endl;
		} else {
			cout << "kmer found!" << endl;
			DataAccessor<UnitigData>* da = map.getData(); // Get DataAccessor from unitig
			UnitigColors* uc = da->getUnitigColors(map);
		}
	}
}






