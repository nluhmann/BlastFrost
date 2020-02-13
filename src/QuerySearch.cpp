/*
 * QuerySearch.cpp
 *
 *  Created on: 14 Jan 2019
 *      Author: nina
 */






#include "QuerySearch.hpp"

using namespace std;


QuerySearch::QuerySearch(ColoredCDBG<UnitigData>& graph) : cdbg(graph) {
	cout << "QuerySearch initialized!" << endl;
}




//////////////////////////////////////////////////////////////////////////////////////////////
// FUNC 1: Search Bifrost graph for queries, consider K-mer neighborhood up to distance d
/////////////////////////////////////////////////////////////////////////////////////////////

unordered_map<size_t,vector<int>> QuerySearch::search(const string& query, const int k, const int ndistance) const {

	const size_t num_kmers = query.length() - k + 1;
	const char *query_str = query.c_str(); //split query into sequence of kmers (!!! query can contain a kmer multiple times, do not change the order of the kmers at this point!)
	unordered_map<size_t,vector<int>> arr;

	int kmer_count = 0;
	bool first = true;
	bool wasEmpty = false;
	const UnitigColors* old_uc;

	for(KmerIterator it_km(query_str), it_km_end; it_km != it_km_end; ++it_km){ //search each kmer in cdbg and return color set

		UnitigColorMap<UnitigData> map = cdbg.find(it_km->first);

		map.dist = 0;
		map.len = map.size - Kmer::k + 1;
		map.strand = true;

		if (!map.isEmpty) {
			const DataAccessor<UnitigData>* da = map.getData();
			const UnitigColors* uc = da->getUnitigColors(map);

			bool copy = (!first && !wasEmpty && (uc == old_uc));
			old_uc = uc;

			if (copy) { //if this UnitigColors object contains the same colors as the object of the previous kmer (which is likely), then we already know whats happening!
				for(auto& color : arr){
					if (arr[color.first][kmer_count - 1] == 1) arr[color.first][kmer_count] = 1;
				}
			}
			else {
				first = false;

				for (UnitigColors::const_iterator it = uc->begin(map); it != uc->end(); it.nextColor()) {
						const size_t color = it.getColorID();
						const std::unordered_map<size_t, vector<int>>::const_iterator iter = arr.find(color);
						if (iter == arr.end()) arr.insert({color, vector<int>(num_kmers, 0)});
						arr[color][kmer_count] = 1;
				}
			}

			//now check the k-mers neighborhood too, but do not overwrite perfect matches!
			if (ndistance > 0){
				const vector<Kmer> neighborhood = QuerySearch::compute_neighborhood(it_km->first.toString(), ndistance);
				for (const auto& nkmer : neighborhood){
					const const_UnitigColorMap<UnitigData> nmap = cdbg.find(nkmer);

					if (!nmap.isEmpty) {
						const DataAccessor<UnitigData>* nda = nmap.getData();
						const UnitigColors* nuc = nda->getUnitigColors(nmap);

						for (UnitigColors::const_iterator nit = nuc->begin(nmap); nit != nuc->end(); ++nit) {
							const size_t ncolor = nit.getColorID();

							if (nuc -> contains(nmap, ncolor)){
								std::unordered_map<size_t,vector<int>>::iterator niter = arr.find(ncolor);

								if (niter == arr.end()){
									arr.insert({ncolor, vector<int>(num_kmers, 0)});
									arr[ncolor][kmer_count] = 2;
								}
								else if (arr[ncolor][kmer_count] == 0) arr[ncolor][kmer_count] = 2;
							}
						}
					}
				}
			}
		}
		wasEmpty = map.isEmpty;

		++kmer_count;
	}
	return arr;
}




vector<Kmer> QuerySearch::compute_neighborhood(const string& kmer_str, const int d) const {

	static const size_t alphabet_sz = 4;
	static const char alphabet[alphabet_sz] = {'A','C','G','T'};

	vector<Kmer> neighborhood;
	vector<int> firstRow;

	firstRow.reserve(kmer_str.size());

	for (unsigned int i = 0; i <= kmer_str.size(); ++i) firstRow.push_back(i);

	QuerySearch::searchNextRow(string(""), kmer_str, firstRow, neighborhood, alphabet, alphabet_sz, d);

	return neighborhood;

}

void QuerySearch::searchNextRow(const string& v, const string& word, const vector<int>& lastRow, vector<Kmer>& neighborhood, const char* alphabet, const size_t alphabet_sz, const int d) const {

	int min = *(std::min_element(lastRow.begin(), lastRow.end()));

	if ((min <= d) && (min > 0) && (v.length() == word.length())) {

		const string suffix(v.substr(0, word.length()));
		neighborhood.push_back(Kmer(suffix.c_str()));

	} else if (min == d){

		for(const auto& elem : lastRow){
			if (elem == d){
				//report v*w^x
				if (v.length() == word.length()){
					neighborhood.push_back(Kmer(v.c_str()));
					break;
				}
				else if (v.length() < word.length()) {
					const string suffix = word.substr(v.length());
					const string concat = v + suffix;
					neighborhood.push_back(Kmer(concat.c_str()));
					break;
				}

				//we can only search for kmers in neighborhood of same length
				//if(concat.length() == word.length()){
				//	Kmer next_kmer(concat.c_str());
				//	neighborhood.push_back(next_kmer);
				//}
			}
		}
	}
	else if (min < d){

		for (size_t j = 0; j != alphabet_sz; ++j){
			const char* str = word.c_str();

			vector<int> nextRow;
			nextRow.push_back(lastRow[0]+1); //first entry can only be an insert

			for (unsigned int i = 1; i< lastRow.size(); ++i){
				const int ins = lastRow[i]+1;
				const int del = nextRow[i-1]+1;
				const int sub = lastRow[i-1] + (str[i-1] != alphabet[j]);

				nextRow.push_back(std::min(std::min(ins, del), sub));
			}
			const string next = v + alphabet[j];
			QuerySearch::searchNextRow(next, word, nextRow, neighborhood, alphabet,alphabet_sz, d);
		}
	}
}




/*
 * _Approximate_ number of matches and mismatches from k-mer hits.
 * Atm, we use the average number of mismatches that can explain a run of 0's of a specific length.
 */
int QuerySearch::compute_score(const vector<int>& hit) const {

	const int score_match = 1;
	const int score_mismatch = -2;
	const int l = hit.size();

	int mismatch = 0;
	int cnt = 0;
	//int cnt_hit = 0;

	for (const auto elem : hit){

		if (elem == 0) ++cnt;
		else if (cnt != 0){
			//cout << "cnt: " << cnt << endl;
			const int local = ceil(float(cnt)/31);
			const int local2 = floor(cnt - 31 + 1);
			const int avg = (local+local2)/2;

			mismatch += avg;
			cnt = 0;
		}
	}

	if (cnt != 0){
		const int local = ceil(float(cnt)/31);
		const int local2 = floor(cnt - 31 + 1);
		const int avg = (local+local2)/2;
		mismatch += avg;
	}

	const int match = l - mismatch;

	if (match <= 0) cout << "ERROR! No matches!" << endl;

	return match*score_match+mismatch*score_mismatch;

}




//////////////////////////////////////////////////////////////////////////////////////////////
// FUNC 2: Search Bifrost graph for queries, evaluate p-value and report hit Kmers subsequently (first step of subgraph extraction)
/////////////////////////////////////////////////////////////////////////////////////////////






QuerySearch::searchResultSubgraph QuerySearch::search_kmers(const string& query, const int k, const int ndistance, const double& db_size) const {

	const size_t num_kmers = query.length() - k + 1;
	const char *query_str = query.c_str(); //split query into sequence of kmers (!!! query can contain a kmer multiple times, do not change the order of the kmers at this point!)
	unordered_map<size_t,vector<int>> arr; //keep this locally to evaluate p-value, return following map if accepted

	searchResultSubgraph result;

	unsigned int kmer_count = 0;
	//bool first_seq = true;
	unordered_map<size_t,bool> first;
	//bool wasEmpty = false;
	//const UnitigColors* old_uc;

	for(KmerIterator it_km(query_str), it_km_end; it_km != it_km_end; ++it_km){ //search each kmer in cdbg and return color set

		UnitigColorMap<UnitigData> map = cdbg.find(it_km->first);

//		if (!map.isEmpty){
//			cout << map.referenceUnitigToString() << endl;
//			cout << map.size << endl;
//			cout << map.dist << endl;
//			cout << "---" << endl;
//		}


		//map.strand = true;


		if (!map.isEmpty) {

			int dist = map.dist;

			map.dist = 0;
			map.len = map.size - Kmer::k + 1;


			const DataAccessor<UnitigData>* da = map.getData();
			const UnitigColors* uc = da->getUnitigColors(map);

			//bool copy = (!first_seq && !wasEmpty && (uc == old_uc));
			bool copy = false;
			//old_uc = uc;

			//ToDo: CHECK HERE!!!
			Kmer head = map.getMappedHead();
			//Kmer head = map.getUnitigHead();

			if (copy) { //if this UnitigColors object contains the same colors as the object of the previous kmer (which is likely), then we already know whats happening!
				for(auto& color : arr){
					if (arr[color.first][kmer_count - 1] == 1) {
						arr[color.first][kmer_count] = 1;

						if (result.mapping[color.first].size() == 0 || result.mapping[color.first].back().toString() != head.toString()){ //check if this unitig has been found already
							result.mapping[color.first].push_back(head);
						}
					}
				}
			}
			else {
				//first_seq = false;

				for (UnitigColors::const_iterator it = uc->begin(map); it != uc->end(); it.nextColor()) {

						const size_t color = it.getColorID();

						if (first.find(color) == first.end()){
							first[color] = false;
							if (map.strand == 0){
								result.prefix_offset[color] = map.size - Kmer::k + 1 - dist;
							} else {
								result.prefix_offset[color] = dist;
							}
						}

						if (map.strand == 1){
							result.suffix_offset[color] = map.size - Kmer::k + 1 - dist;
						} else {
							result.suffix_offset[color] = dist;
						}



						const std::unordered_map<size_t, vector<int>>::const_iterator iter = arr.find(color);
						if (iter == arr.end()){
							arr.insert({color, vector<int>(num_kmers, 0)});

							vector<Kmer> newset;
							result.mapping.insert({color, newset});
						}



						//check if this 1 ends a stretch of 0's, in which case the stretch of 0's has to be at least of length k!
						bool ok = true;
						if (kmer_count > Kmer::k && arr[color][kmer_count-1] == 0){
							for(unsigned int i=2; i<=Kmer::k; i++){
								if (arr[color][kmer_count-i] != 0){
									ok = false;
									break;
								}
							}
						}



						if (ok){
							arr[color][kmer_count] = 1;

							if (result.mapping[color].size() == 0) {
								result.mapping[color].push_back(head);
							} else if (result.mapping[color].back().toString().compare(head.toString()) != 0){
								result.mapping[color].push_back(head);
							}


//							if (result.mapping[color].size() == 0) {
//								result.mapping[color].push_back(head);
//								//arr[color][kmer_count] = 1;
////							} else if (result.mapping[color].back().toString().compare(head.toString()) == 0) {
////								arr[color][kmer_count] = 1;
//							} else if (arr[color][kmer_count-1] == 1){
//
//								if (result.mapping[color].back().toString() != head.toString()){
//
//									Kmer previous = result.mapping[color].back();
//									UnitigColorMap<UnitigData> ucm = cdbg.find(previous);
//
//									for (const auto& successor : ucm.getSuccessors()){
//										const DataAccessor<UnitigData>* da = successor.getData();
//										const UnitigColors* uc = da->getUnitigColors(successor);
//
//										Kmer prev_head = successor.getMappedHead();
//										Kmer alternative = successor.getUnitigHead();
//
//										if (prev_head.toString().compare(head.toString()) == 0 || alternative.toString().compare(head.toString()) == 0){
//											result.mapping[color].push_back(head);
//											if (cdbg.getColorName(color) == "assemblies/SAL_WA5226AA_AS.scaffold.fasta"){
//												cout << "cont" << endl;
//											}
//										} else {
//											if (cdbg.getColorName(color) == "assemblies/SAL_WA5226AA_AS.scaffold.fasta"){
//												cout << "not continuous!" << endl;
//											}
//										}
//									}
//								} else{
//									if (cdbg.getColorName(color) == "assemblies/SAL_WA5226AA_AS.scaffold.fasta"){
//										cout << "same unitig!" << endl;
//									}
//								}
//							}
//							} else {
//								//result.mapping[color].push_back(head);
//								arr[color][kmer_count] = 1;
//							}
						}
				}
			}

			//now check the k-mers neighborhood too, but do not overwrite perfect matches!
			if (ndistance > 0){
				const vector<Kmer> neighborhood = QuerySearch::compute_neighborhood(it_km->first.toString(), ndistance);
				//bool remove_only_once = true;
				for (const auto& nkmer : neighborhood){
					UnitigColorMap<UnitigData> nmap = cdbg.find(nkmer);

					int ndist = nmap.dist;

					nmap.dist = 0;
					nmap.len = nmap.size - Kmer::k + 1;

					if (!nmap.isEmpty) {
						const DataAccessor<UnitigData>* nda = nmap.getData();
						const UnitigColors* nuc = nda->getUnitigColors(nmap);

						Kmer nhead = nmap.getMappedHead();

						for (UnitigColors::const_iterator nit = nuc->begin(nmap); nit != nuc->end(); ++nit) {
							const size_t ncolor = nit.getColorID();

							if (nuc -> contains(nmap, ncolor)){
								std::unordered_map<size_t,vector<int>>::iterator niter = arr.find(ncolor);

//								if (niter == arr.end()){
//									arr.insert({ncolor, vector<int>(num_kmers, 0)});
//
//									vector<Kmer> newset;
//									result.mapping.insert({ncolor, newset});
//								}
//								if (arr[ncolor][kmer_count] == 0){
//									//toDo: we need to check if the current unitig is actually a successor of the previous one, if the previous one had a direct hit...
//									arr[ncolor][kmer_count] = 2;
//
//									if (result.mapping[ncolor].size() == 0 || result.mapping[ncolor].back().toString().compare(nhead.toString()) != 0 ){ //check if this unitig has been found already
//										result.mapping[ncolor].push_back(nhead);
//									}
//
//									if (first.find(ncolor) == first.end()){
//										first[ncolor] = false;
//										if (nmap.strand == 0){
//											result.prefix_offset[ncolor] = nmap.size - Kmer::k + 1 - ndist;
//										} else {
//											result.prefix_offset[ncolor] = ndist;
//										}
//									}
//
//									if (nmap.strand == 1){
//											result.suffix_offset[ncolor] = nmap.size - Kmer::k + 1 - ndist;
//									} else {
//											result.suffix_offset[ncolor] = ndist;
//									}
//								} else if (arr[ncolor][kmer_count] == 2){
//									//todo: check if these kmers hit the same unitig?
//									cout << "same neighborhood kmer hit" << endl;
//								}

								if (niter == arr.end()){
									arr.insert({ncolor, vector<int>(num_kmers, 0)});

									vector<Kmer> newset;
									result.mapping.insert({ncolor, newset});
								}
								if (arr[ncolor][kmer_count] == 0){
									arr[ncolor][kmer_count] = 2;

									bool perfect_hit = true;
									for (auto& elem : arr[ncolor]){
										if (elem == 1){
											perfect_hit = true;
											break;
										}
									}

									if (perfect_hit){

										if (result.mapping[ncolor].size() == 0 || result.mapping[ncolor].back().toString().compare(nhead.toString()) != 0 ){ //check if this unitig has been found already

											if (result.mapping[ncolor].size() != 0 && (arr[ncolor][kmer_count-1] == 1 || arr[ncolor][kmer_count-1] == 2)){
												Kmer previous = result.mapping[ncolor].back();
												UnitigColorMap<UnitigData> ucm = cdbg.find(previous);

												for (const auto& successor : ucm.getSuccessors()){
													//const DataAccessor<UnitigData>* da = successor.getData();
													//const UnitigColors* uc = da->getUnitigColors(successor);

													Kmer prev_head = successor.getMappedHead();
													Kmer alternative = successor.getUnitigHead();

													if (prev_head.toString().compare(nhead.toString()) == 0 || alternative.toString().compare(nhead.toString()) == 0){
														result.mapping[ncolor].push_back(nhead);

													}
												}
											} else if (result.mapping[ncolor].size() == 0) {
												result.mapping[ncolor].push_back(nhead);
											}
										}

										if (first.find(ncolor) == first.end()){
											first[ncolor] = false;
											if (nmap.strand == 0){
												result.prefix_offset[ncolor] = nmap.size - Kmer::k + 1 - ndist;
											} else {
												result.prefix_offset[ncolor] = ndist;
											}
										}

										if (nmap.strand == 1){
											result.suffix_offset[ncolor] = nmap.size - Kmer::k + 1 - ndist;
										} else {
											result.suffix_offset[ncolor] = ndist;
										}
									}

//								} else if (arr[ncolor][kmer_count] == 2){
//									if (remove_only_once && result.mapping[ncolor].size() != 0){
//										result.mapping[ncolor].pop_back();
//										remove_only_once = false;
//									}

								}

							}
						}
					}
				}
			}
		}
		//wasEmpty = map.isEmpty;

		++kmer_count;
	}


	//for each color in arr, compute p-value, remove color from mapping if neccessary
	for (auto& color : arr){

		const int score = compute_score(color.second);
		const int len = color.second.size();


		long double pvalue2 = 0;
		if (score != len){
			const long double log_evalue = compute_log_evalue(score,db_size,len);
			const long double log_pvalue = compute_log_pvalue(log_evalue);
			pvalue2 = pow(10,log_pvalue);
		}


		if (pvalue2 >= 0.05){
			//remove color from mapping
			result.mapping.erase(color.first);

		} else {

			//record number of missed k-mers in beginning and end of query
			//int cnt_prefix = 0;
			for(auto& elem : color.second){
				if (elem == 1){
					break;
				}
				result.prefix_missing[color.first]++;
			}

			//int cnt_suffix = 0;
			std::reverse(color.second.begin(),color.second.end());
			for(auto& elem : color.second){
				if (elem == 1){
					break;
				}
				result.suffix_missing[color.first]++;
			}

		}



	}

	return result;
}





//void GraphTraverser::test_neighborhood_func(){
//	//test neighborhood function!
//	//string test = "ACA";
//	//vector<Kmer> neighborhood = compute_neighborhood(test, 2);
//	//cout << "Neighbohood of: " << test << endl;
//	//for (auto& elem: neighborhood){
//	//	cout << elem.toString() << endl;
//	//}
//	//cout << neighborhood.size() << endl;
//}

//DEPRECATED.
//unordered_map<size_t,double> GraphTraverser::compute_significance(const unordered_map<size_t,vector<int>>& hits, const long double p) const {
//	unordered_map<size_t,double> p_values;
//
//	for (const auto& hit: hits){
//
//		const int q = hit.second.size();
//
//		//compute number of 1's in vector
//		const int m = std::count(hit.second.begin(), hit.second.end(), 1);
//
//		const double r = gsl_cdf_binomial_Q(m, p, q);
//
//		if (r > 0.05) cout << "p-value: " << r << endl;
//
//		p_values[hit.first] = r;
//	}
//
//	return p_values;
//}

