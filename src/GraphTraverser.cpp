/*
 * GraphTraverser.cpp
 *
 *  Created on: 26 Feb 2018
 *      Author: nina
 */

#include "GraphTraverser.hpp"

using namespace std;






GraphTraverser::GraphTraverser(ColoredCDBG<UnitigData>& graph) :
	cdbg(graph) {
	cout << "GraphTraverser initialized!" << endl;
}




unordered_map<size_t,vector<int>> GraphTraverser::search(string query, int k){

	vector<Kmer> kmers;

	//split query into sequence of kmers (!!! query can contain a kmer multiple times, do not change the order of the kmers at this point!)
	for(int i = 0; i< (query.length()-k+1); ++i){
	    string kmer = query.substr(i,k);
	    //cout << kmer << " " << kmer.length() << endl;
	    const char *cstr = kmer.c_str();
	    Kmer next(cstr);
	    kmers.push_back(next);
	}

	size_t num_kmers = kmers.size();

	unordered_map<size_t,vector<int>> arr;

	//search each kmer in cdbg and return color set
	int kmer_count = 0;
	for (auto& kmer: kmers){
		UnitigMap<DataAccessor<UnitigData>, DataStorage<UnitigData>, false> map = cdbg.find(kmer);
		//set<string> colors;
		if (map.isEmpty) {
			//cout << "kmer not found" << endl;
			//cout << kmer.toString() << endl;
		} else {
			DataAccessor<UnitigData>* da = map.getData();
			UnitigColors* uc = da->getUnitigColors(map);

			for(UnitigColors::const_iterator it = uc->begin(map); it != uc->end(); it.nextColor()) {
			  size_t color = it.getColorID();
				//note to self: the iterator goes through all colors of the unitig, but we want to only keep the ones that the kmer is really annotated with
				if (uc -> contains(map, color)){
					//colors.insert(cdbg.getColorName(color));
					std::unordered_map<size_t,vector<int>>::iterator iter = arr.find(color);

					if (iter == arr.end()){
						vector<int> vec(num_kmers, 0);
						arr.insert({color,vec});
					}
					arr[color][kmer_count] = 1;
				}
			}
		}
		kmer_count++;
	}
	return arr;
}


//ToDo: Debug!
void GraphTraverser::remove_singletonHits(unordered_map<size_t,vector<int>>& hits, int k){
	for (auto& hit : hits){
		vector<int> seq = hit.second;
		//if first and last appearance of '1' in vector is more than k appart, then everything is ok!
		int cnt_begin = 0;
		for (vector<int>::iterator i = seq.begin(); i != seq.end(); ++i){
			if (*i == 1){
				break;
			}
			cnt_begin++;
		}

		int cnt_end = 0;
		for (vector<int>::reverse_iterator i = seq.rbegin(); i != seq.rend(); ++i){
			if (*i == 1){
				break;
			}
			cnt_end++;
		}

		//cout << hit.first << endl;
		//cout << cnt_begin << endl;
		//cout << cnt_end << endl;
		//cout << seq.size() << endl;
		if ((seq.size() - cnt_end - cnt_begin) <= k){
			hits.erase(hit.first);
			cout << "erase " << hit.first << endl;
		}
	}
}




long double GraphTraverser::compute_p(long double& k, long double& x, long double& sigma){
	long double p = 1-pow(1-pow(sigma,-k),x);
	return p;
}



unordered_map<size_t,double> GraphTraverser::compute_significance(unordered_map<size_t,vector<int>>& hits, long double& p) {
	unordered_map<size_t,double> p_values;

	for (auto& hit: hits){
		int q = hit.second.size();

		//compute number of 1's in vector
		int m = std::count(hit.second.begin(), hit.second.end(), 1);

		double r = gsl_cdf_binomial_Q(m, p, q);

		if (r > 0.05){
			cout << "p-value: " << r << endl;
		}
		p_values[hit.first] = r;
	}

	return p_values;
}


int GraphTraverser::compute_score(vector<int>& hit){
	int score_match = 1;
	int score_mismatch = -2;


	int l = hit.size();
		//cout << "length: " << l << endl;
	int mismatch = 0;
	int cnt = 0;
	for(auto& elem: hit){
		if (elem == 0){
			cnt+=1;
		} else if (cnt != 0){
			//cout << "cnt: " << cnt << endl;
			int local = ceil(float(cnt)/31);
			int local2 = floor(cnt - 31 + 1);
			int avg = (local+local2)/2;
			mismatch += avg;
			cnt = 0;
		}
	}

	if (cnt != 0){
		int local = ceil(float(cnt)/31);
		int local2 = floor(cnt - 31 + 1);
		int avg = (local+local2)/2;
		mismatch += avg;
	}

	int match = l - mismatch;
	//cout << "match: " << match << endl;
	if (match <= 0){
		cout << "ERROR! No matches!" << endl;
	}

	int score = match*score_match+mismatch*score_mismatch;
	return score;
}




long double GraphTraverser::compute_evalue(int score, double db_size, int n){
	long double lambda = 1.330;
	float k = 0.621;

	long double evalue = k*db_size*n*exp(-lambda*score);

	return evalue;
}

long double GraphTraverser::compute_pvalue(float evalue){
	long double pvalue = 1-exp(-evalue);
	return pvalue;
}


long double GraphTraverser::compute_log_evalue(int score, double db_size, int n){
	long double lambda = 1.330;
	float k = 0.621;

	long double log_evalue = round(log(k*db_size*n)-lambda*score);

	return log_evalue;

}

long double GraphTraverser::compute_log_pvalue(long double log_evalue){
	long double evalue = pow(10,log_evalue);
	if (1-exp(-evalue) > 0){
		return round(log(1-exp(-evalue)));
	} else {
		return round(log_evalue);
	}
}





void GraphTraverser::writePresenceMatrix(unordered_map<size_t,vector<int>>& arr, string& outfile){
	std::ofstream output(outfile,std::ofstream::binary);

	for (auto& elem : arr){
		string color = cdbg.getColorName(elem.first);
		output << color;
		for (auto& p : elem.second){
			output << "\t" << p;
		}
		output << endl;
	}
	output.close();
}





