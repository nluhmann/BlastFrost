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
#include <mutex>
#include <stack>
#include <queue>

#include <gsl/gsl_cdf.h>
#include <math.h>

class GraphTraverser{

private:

	//reference to the graph that is going to be traversed
	ColoredCDBG<UnitigData>& cdbg;
	const long double lambda = 1.330;
	const float blast_k = 0.621;


public:

	GraphTraverser(ColoredCDBG<UnitigData>& cdbg);

	unordered_map<size_t,vector<int>> search(const string& query, const int k, const int ndistance, const string& query_name) const;

	//void writeKmerPresence(vector<pair<Kmer,set<string>>> results, string& resfile);

	void writePresenceMatrix(const unordered_map<size_t,vector<int>>& arr, const string& outfile, const unordered_map<size_t,long double>& pvalues) const;

	//long double compute_p(long double& k, long double& x, long double& sigma);

	unordered_map<size_t,double> compute_significance(const unordered_map<size_t,vector<int>>& hits, const long double p) const;

	void remove_singletonHits(unordered_map<size_t,vector<int>>& hits) const;

	int compute_score(const vector<int>& hit) const;

	vector<Kmer> compute_neighborhood(const string& kmer_str, const int d) const;

	void searchNextRow(const string& v, const string& word, const vector<int>& lastRow, vector<Kmer>& neighborhood, const char* alphabet, const size_t alphabet_sz, const int d) const;

	void extractSubGraph(const string& query, const int k, const int distance);

	void pathLength(const unordered_map<size_t,vector<Kmer>>& all_paths, const int& ref_length);

	/*
	 * Find the unitig that corresponds to this string, and report all colors of this unitig
	 */
	vector<string> getColors(const string& u) const;

	void exploreSubgraph(const string& s) const;

	void exploreBubble(const string& left, const string& right, const int threshold);

	vector<vector<UnitigColorMap<UnitigData>>> DFS_Iterative(const UnitigColorMap<UnitigData>& start, const Kmer& stop, const int threshold);

	void printPaths(const vector<vector<UnitigColorMap<UnitigData>>> allPaths);

	//void BFS_Iterative(const UnitigColorMap<UnitigData>& start, const Kmer& stop, const int threshold);


	inline long double compute_p(const long double k, const long double x, const long double sigma) const {

		return 1-pow(1-pow(sigma,-k),x);
	}

	inline long double compute_evalue(const int score, const double db_size, const int n) const {

		return blast_k * db_size * n * exp(-lambda * score);
	}

	inline long double compute_pvalue(const long double evalue) const {

		return 1-exp(-evalue);
	}


	inline long double compute_log_evalue(const int score, const double db_size, const int n) const {

		return round(log(blast_k * db_size * n) - lambda * score);
	}

	inline long double compute_log_pvalue(const long double log_evalue) const {

		const long double evalue = pow(10, log_evalue);

		if (1-exp(-evalue) > 0) return round(log(1-exp(-evalue)));

		return round(log_evalue);
	}

};








#endif /* GRAPHTRAVERSER_HPP_ */
