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





class SubGraphTraverser{

private:

	//reference to the graph that is going to be traversed
	ColoredCDBG<UnitigData>& cdbg;


public:

	SubGraphTraverser(ColoredCDBG<UnitigData>& cdbg);


	void extractSubGraph(const string& query, const int k, const int distance);

	void pathLength(const unordered_map<size_t,vector<Kmer>>& all_paths, const int& ref_length);

	/*
	 * Find the unitig that corresponds to this string, and report all colors of this unitig
	 */
	vector<string> getColors(const string& u) const;

	void exploreSubgraph(const string& s) const;





};








#endif /* GRAPHTRAVERSER_HPP_ */
