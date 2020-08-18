/*
 * GraphTraverser.hpp
 *
 *  Created on: 26 Feb 2018
 *      Author: nina
 */

#ifndef GRAPHTRAVERSER_HPP_
#define GRAPHTRAVERSER_HPP_



#include <ColoredCDBG.hpp>
#include "UnitigData.hpp"
#include "QuerySearch.hpp"





class SubGraphTraverser{

private:

	//reference to the graph that is going to be traversed
	ColoredCDBG<UnitigData>& cdbg;
	double& db_size;
	QuerySearch& que;


public:

	struct subgraphs {
			unordered_map<size_t,vector<size_t>> colors;
			unordered_map<size_t,vector<std::string>> sequences;
			//searchResultSubgraph() : prefix_offset(0), suffix_offset(0), prefix_missing(0), suffix_missing(0) {};
		};




	SubGraphTraverser(ColoredCDBG<UnitigData>& cdbg, double& db_size, QuerySearch& q);

	unordered_map<size_t,vector<std::string>> extractSubGraph(const string& query, const int k, const int distance);

	unordered_map<size_t, vector<size_t>> groupSeedHits(unordered_map<size_t,vector<Kmer>>& map);

	subgraphs extractSubGraph_intelligent(const string& query, const int k, const int distance);

	void pathLength(const unordered_map<size_t,vector<Kmer>>& all_paths);

	vector<string> pathSequence(const vector<Kmer>& paths, const int& diff_prefix, const int& diff_suffix);

	void printPaths(string& outprefix, string& query, unordered_map<size_t,vector<std::string>>& paths, const string& queryID);

	void printPaths_intelligent(string& outprefix, string& query, subgraphs result, const string& queryID);

	void testPath(const unordered_map<size_t,vector<Kmer>>& all_paths);



	/*
	 * Find the unitig that corresponds to this string, and report all colors of this unitig
	 */
	vector<string> getColors(const string& u) const;

	void exploreSubgraph(const string& s) const;





};








#endif /* GRAPHTRAVERSER_HPP_ */
