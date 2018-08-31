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

	unordered_map<size_t,vector<int>> search(const string& query, const int k, const int ndistance) const;

	//void writeKmerPresence(vector<pair<Kmer,set<string>>> results, string& resfile);

	void writePresenceMatrix(const unordered_map<size_t,vector<int>>& arr, const string& outfile, const unordered_map<size_t,long double>& pvalues);

	long double compute_p(long double& k, long double& x, long double& sigma);

	unordered_map<size_t,double> compute_significance(unordered_map<size_t,vector<int>>& hits, long double& p);

	void remove_singletonHits(unordered_map<size_t,vector<int>>& hits);

	int compute_score(const vector<int>& hit);

	long double compute_evalue(const int& score, const double& db_size, const int& n) const;

	long double compute_pvalue(const long double& evalue) const;

	long double compute_log_evalue(const int& score, const double& db_size, const int& n) const;

	long double compute_log_pvalue(const long double& log_evalue) const;

	vector<Kmer> compute_neighborhood(const string& kmer_str, const int d) const;

	void searchNextRow(const string v, const string& word, const vector<int>& lastRow, vector<Kmer>& neighborhood, const vector<char>& alphabet, const int d) const;

	/*
	 * Find the unitig that corresponds to this string, and report all colors of this unitig
	 */
	vector<string> getColors(const string& u);

	void exploreSubgraph(const string& s) const;

	void exploreBubble(const string& left, const string& right, const int threshold);

	void DFS_Iterative(const UnitigColorMap<UnitigData>& start, const Kmer& stop, const int threshold);

};








#endif /* GRAPHTRAVERSER_HPP_ */
