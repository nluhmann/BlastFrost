/*
 * GraphTraverser.hpp
 *
 *  Created on: 26 Feb 2018
 *      Author: nina
 */

#ifndef GRAPHTRAVERSER_HPP_
#define GRAPHTRAVERSER_HPP_



#include <bifrost/ColoredCDBG.hpp>
#include "UnitigData.hpp"


#include <string>
#include <unordered_map>
#include <algorithm>
#include <vector>
#include <iterator>
#include <set>
#include <fstream>

class GraphTraverser{
private:

	//reference to the graph that is going to be traversed
	ColoredCDBG<UnitigData>& cdbg;


public:

	GraphTraverser(ColoredCDBG<UnitigData>& cdbg);

	vector<pair<Kmer,set<string>>> search(string query, int k);

	void writeKmerPresence(vector<pair<Kmer,set<string>>> results, string& resfile);

};








#endif /* GRAPHTRAVERSER_HPP_ */
