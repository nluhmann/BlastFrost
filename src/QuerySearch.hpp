/*
 * QuerySearch.hpp
 *
 *  Created on: 14 Jan 2019
 *      Author: nina
 */

#ifndef QUERYSEARCH_HPP_
#define QUERYSEARCH_HPP_




#include <ColoredCDBG.hpp>
#include "UnitigData.hpp"


class QuerySearch{


private:

	//reference to the graph that is going to be traversed
	ColoredCDBG<UnitigData>& cdbg;
	const long double lambda = 1.330;
	const float blast_k = 0.621;


public:

	struct searchResultSubgraph {
		unordered_map<size_t,vector<Kmer>> mapping;
		unordered_map<size_t, int> prefix_offset;
		unordered_map<size_t, int> suffix_offset;
		unordered_map<size_t, int> prefix_missing;
		unordered_map<size_t, int> suffix_missing;

		//searchResultSubgraph() : prefix_offset(0), suffix_offset(0), prefix_missing(0), suffix_missing(0) {};
	};


	QuerySearch(ColoredCDBG<UnitigData>& cdbg);


	unordered_map<size_t,vector<int>> search(const string& query, const int k, const int ndistance) const;

	int compute_score(const vector<int>& hit) const;

	vector<Kmer> compute_neighborhood(const string& kmer_str, const int d) const;

	void searchNextRow(const string& v, const string& word, const vector<int>& lastRow, vector<Kmer>& neighborhood, const char* alphabet, const size_t alphabet_sz, const int d) const;

	searchResultSubgraph search_kmers(const string& query, const int k, const int ndistance, const double& db_size) const;


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



#endif /* QUERYSEARCH_HPP_ */
