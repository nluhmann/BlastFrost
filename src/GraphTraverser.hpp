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

#include <gsl/gsl_cdf.h>
#include <math.h>

class GraphTraverser{
private:

	//reference to the graph that is going to be traversed
	ColoredCDBG<UnitigData>& cdbg;


public:

	GraphTraverser(ColoredCDBG<UnitigData>& cdbg);

	unordered_map<size_t,vector<int>> search(string query, int k);

	//void writeKmerPresence(vector<pair<Kmer,set<string>>> results, string& resfile);

	void writePresenceMatrix(unordered_map<size_t,vector<int>>& arr, string& outfile);

	long double compute_p(long double& k, long double& x, long double& sigma);

	unordered_map<size_t,double> compute_significance(unordered_map<size_t,vector<int>>& hits, long double& p);

	void remove_singletonHits(unordered_map<size_t,vector<int>>& hits, int k);

	int compute_score(vector<int>& hit);

	long double compute_evalue(int score, double db_size, int n);

	long double compute_pvalue(float evalue);

	long double compute_log_evalue(int score, double db_size, int n);

	long double compute_log_pvalue(long double log_evalue);

};








#endif /* GRAPHTRAVERSER_HPP_ */
