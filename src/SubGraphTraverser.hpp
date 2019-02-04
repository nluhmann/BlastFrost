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


	void extractSubGraph(const string& query, const int k, const int distance, string& outprefix, string& queryfile);

	void pathLength(const unordered_map<size_t,vector<Kmer>>& all_paths, const int& ref_length);

	unordered_map<size_t,std::string> pathSequence(const unordered_map<size_t,vector<Kmer>>& all_paths);

	void printPaths(string& outprefix, string& query, unordered_map<size_t,std::string>& paths);

	void testPath(const unordered_map<size_t,vector<Kmer>>& all_paths);



	/*
	 * Find the unitig that corresponds to this string, and report all colors of this unitig
	 */
	vector<string> getColors(const string& u) const;

	void exploreSubgraph(const string& s) const;





};








#endif /* GRAPHTRAVERSER_HPP_ */
