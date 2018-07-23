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

	//unordered_map<size_t,vector<int>> search(string query, int k) const;

	unordered_map<size_t,vector<int>> search2(string query, int k);

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

	vector<string> compute_neighborhood(string kmer_str, int d);

	void searchNextRow(string v, string& word, vector<int> lastRow, vector<string>& neighborhood, vector<char>& alphabet, int& d);

};








#endif /* GRAPHTRAVERSER_HPP_ */
