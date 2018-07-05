
#include <string>
#include <fstream>
#include <cstdio>
#include <thread>
#include <queue>
#include <mutex>
#include <atomic>

#include <bifrost/ColoredCDBG.hpp>
#include "GraphTraverser.hpp"
#include "UnitigData.hpp"


#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <math.h>






using namespace std;



std::mutex mtx_queue;           // mutex for critical section
std::mutex mtx_writing;
std::atomic<int> running(0);



void PrintUsage() {
	cout << "BlastFrost -f <inputFiles> -q <query_sequences> -o <outfile_prefix> " << endl << endl;

	cout << "Mandatory parameters with required argument:" << endl << endl
			<< "  -f,         Input sequence files (FASTA or FASTQ, possibly gziped) and/or graph files (GFA)"
			<< endl
			<< "  -q,         Query sequences (multiple FASTA)"
			<< endl
			<< "  -o,         Prefix for output files (default: 'output')"

			<< endl << endl << "Optional parameters:"<< endl << endl
			<< "  -t,             Number of threads (default is 1)"
			<< endl
			<< "  -k,             Length of k-mers (default is 31, max. is 31)"
			<< endl
			<< "  -i,           Clip tips shorter than k k-mers in length"
			<< endl
			<< "  -v,             Print information messages during construction"
			<< endl << endl;
}


struct searchOptions {
  string graphfile = "";
  vector<string> queryfiles;
  size_t nb_threads = 1;
  bool verbose = false;
  string outprefix = "output";
} ;



void parseArgumentsNew(int argc, char **argv, searchOptions& opt) {

	int oc; //option character

	while ((oc = getopt(argc, argv, "o:t:q:g:v")) != -1) {
		switch (oc) {
		case 'o':
			/* handle -o, set fileprefix */
			opt.outprefix = optarg;
			break;
		case 't':
			/* handle -t, set number of threads */
			opt.nb_threads = atoi(optarg);
			break;
		case 'v':
			//handle -v, verbose mode
			opt.verbose = true;
			break;
		case 'g':
			opt.graphfile = optarg;
			break;
		case 'q':
			/* handle -g, collect query files */
			for (optind--; (optind < argc) && (*argv[optind] != '-'); ++optind) {
				opt.queryfiles.push_back(argv[optind]);
			}
			break;
		case ':':
			cout << "Invalid option" << endl; /* ToDo: proper exception */
			break;
		case '?':
			cout << "Invalid option" << endl;
			break;
		default:
			break; /* ToDo: proper exception */
		}
	}


	//check if precomputed graph prefix is given (-g)
	if (opt.graphfile == ""){
		cout << "No input files given to build graph!" << endl;
		 exit (EXIT_FAILURE);
	}

	//check if output prefix is given (-o), otherwise set default to "output"
	if (opt.outprefix == ""){
		cout << "No outfile prefix given, set default value 'output'" << endl;
		//opt.prefixFilenameOut = "output";
	}

	//check if query file is given (-q)
	if (opt.queryfiles.empty()){
		cout << "No query files given to search graph!" << endl;
		exit (EXIT_FAILURE);
	}
}



int estimate_avg_genomeSize(const vector<string>& files){
	int cnt = 0;
	int size = 0;
	while (cnt < 10){

		ifstream file(files[cnt], ios::binary | ios::ate);
		if (size == 0){
			size = file.tellg();
		} else {
			size = (size + file.tellg())/2;
		}
		cnt++;
	}
	cout << "Avg genome size: " << size << endl;
	return size;
}


string runLength_encode(const vector<int>& v){
	int cnt1 = 0;
	int cnt0 = 0;
	string res;
	string n;
	for (auto& elem : v){
		if (elem == 1){
			if (cnt1 == 0 && cnt0 != 0){
				//we start a new run of 1's
				n = "0:"+std::to_string(cnt0)+",";
				res += n;
				cnt0 = 0;
			}
			cnt1++;
		} else if (elem == 0){
			if (cnt0 == 0 && cnt1 != 0){
				//we start a new run of 0's
				n = "1:"+std::to_string(cnt1)+",";
				res += n;
				cnt1 = 0;
			}
			cnt0++;
		} else {
			cout << "ERROR" << endl;
		}
	}
	if (cnt0 != 0){
		n = "0:"+std::to_string(cnt0)+",";
	} else {
		n = "1:"+std::to_string(cnt1)+",";
	}
	res += n;
	return res;
}


vector<pair<string,string>> parseFasta(const string& queryfile){

	//create and check input stream
	ifstream input(queryfile);
	if(!input.good()){
		cout << "Cannot open query file." << endl;
		exit (EXIT_FAILURE);
	}

	vector<pair<string,string>> fasta;

	string line;
	string name;
	string content;
	while( std::getline(input,line).good() ){
		if(line.empty() || line[0] == '>'){
			if(!name.empty()){
				fasta.push_back(std::make_pair(name,content));
				name.clear();
			}
			if (!line.empty()){
				string token = line.substr(0, line.find(' '));
				name = token.substr(1);
			}
			content.clear();
		} else if(!name.empty()){
			if(line.find(' ') != std::string::npos){
				name.clear();
				content.clear();
			} else {
				content += line;
			}
		}
	}
	if(!name.empty()){
		fasta.push_back(std::make_pair(name,content));
	}

	input.close();
	return fasta;
}



struct searchResult {
  string queryfile;
  long double evalue;
  long double pvalue;
  string query;
  size_t colorID;
  vector<int> hitrun;
} ;




void run_subsample_partialQuery(queue<pair<string,vector<searchResult>>>& q, const double& db_size, const int& k, const vector<pair<string,string>>& fasta, GraphTraverser& tra, int start, int end, const string query){
		running++;
		const int sigma = 4;

		if (end == 0){
			end = fasta.size()-1;
		}

		vector<searchResult> results;
		for(int i = start; i <= end; i++) {

			pair<string,string> seq = fasta[i];
			unordered_map<size_t,vector<int>> res = tra.search(seq.second, k);
			tra.remove_singletonHits(res);


			for (auto& hit : res){
				const int score = tra.compute_score(hit.second);
				const int len = hit.second.size();


				const long double log_evalue = tra.compute_log_evalue(score,db_size,len);
				const long double log_pvalue = tra.compute_log_pvalue(log_evalue);
				const long double evalue2 = pow(10,log_evalue);
				const long double pvalue2 = pow(10,log_pvalue);


				if (pvalue2 <= 0.05){

					//create new searchResult and push to all results
					searchResult result;
					result.evalue = evalue2;
					result.pvalue = pvalue2;
					result.query = seq.first;
					result.colorID = hit.first;
					result.hitrun = hit.second;
					results.push_back(result);

				}

			}
		}

		//we've found all results for this query, push to writing queue!
		mtx_queue.lock();
		q.push(std::make_pair(query,results));
		mtx_queue.unlock();

	running--;
}


void run_subsample_completeQuery(queue<pair<string,vector<searchResult>>>& q, const double& db_size, const int& k, const string queryfile, GraphTraverser& tra, string& outprefix){
	running++;

	const int sigma = 4;

	//ToDo: is it ok if each thread opens a different file?
	//parse fasta file
	vector<pair<string,string>> fasta = parseFasta(queryfile);

	std::string delim = "/";
	auto start = 0U;
	auto end = queryfile.find(delim);
	while (end != std::string::npos){
		start = end + delim.length();
		end = queryfile.find(delim, start);
	}

	string query = queryfile.substr(start, end);

	//delete potentially existing search file for this query
	string f = outprefix+"_"+query+".search";
	if (std::remove(f.c_str()) != 0) {
		perror("Error deleting file");
	} else {
		cout << "File " << f << " removed" << endl;
	}


	vector<searchResult> results;
	for(auto& seq: fasta) {
		unordered_map<size_t,vector<int>> res = tra.search(seq.second, k);
		tra.remove_singletonHits(res);


		for (auto& hit : res){
			const int score = tra.compute_score(hit.second);
			const int len = hit.second.size();


			const long double log_evalue = tra.compute_log_evalue(score,db_size,len);
			const long double log_pvalue = tra.compute_log_pvalue(log_evalue);
			const long double evalue2 = pow(10,log_evalue);
			const long double pvalue2 = pow(10,log_pvalue);


			if (pvalue2 <= 0.05){

				//create new searchResult and push to all results
				searchResult result;
				result.evalue = evalue2;
				result.pvalue = pvalue2;
				result.query = seq.first;
				result.colorID = hit.first;
				result.hitrun = hit.second;
				results.push_back(result);

			}

		}

	}
	//we've found all results for this query, push to writing queue!
	mtx_queue.lock();
	q.push(std::make_pair(query,results));
	mtx_queue.unlock();

	running--;
}



void writeResults(string& outprefix, queue<pair<string,vector<searchResult>>>& q, ColoredCDBG<UnitigData>& cdbg){
	//the main thread will take care of writing the output as soon as it is available


	mtx_writing.lock();
 	while(running > 0 || ! q.empty()){//there are still threads running
		if (! q.empty()){
			mtx_queue.lock();
			pair<string,vector<searchResult>> item = q.front();
			q.pop();
			mtx_queue.unlock();

			//ToDo: If we split the search for a file, we will need to append to an existing file here...
			std::ofstream output(outprefix+"_"+item.first+".search",std::ios_base::app);

			for (auto& result : item.second){
				string color = cdbg.getColorName(result.colorID);
				output << result.query << "\t" << color << "\t" << result.pvalue << "\t";
				string res = runLength_encode(result.hitrun);
				output << res << endl;
			}
			output.close();
		}
	}
	mtx_writing.unlock();
}




int main(int argc, char **argv) {

	if (argc < 2) {
		PrintUsage();
	} else {

		const clock_t load_time = clock();
		searchOptions opt;
		parseArgumentsNew(argc, argv, opt);
		ColoredCDBG<UnitigData> cdbg;

		if (cdbg.read(opt.graphfile, opt.nb_threads, opt.verbose)){
				cout << "Graph loading successful" << endl;
		} else {
			cout << "Graph could not be loaded! Exit." << endl;
			exit (EXIT_FAILURE);
		}


		cout << "Loading took " << (float( clock() - load_time ) /  CLOCKS_PER_SEC) << "sec." << endl;


		////////////////////////////////////////////////////////////////////////////////////////////////////////////

		const clock_t begin_time = clock();


		cout << "---Query graph for input sequences---" << endl;
		GraphTraverser tra(cdbg);

		//const int avg_genomeSize = estimate_avg_genomeSize(opt.filename_seq_in);
		const int avg_genomeSize = cdbg.size();
		cout << avg_genomeSize << endl;
		const size_t numberStrains = cdbg.getNbColors();
		double db_size = numberStrains*avg_genomeSize;

		int k = cdbg.getK();

		//!!! one thread is used by main!!!
		size_t num_threads = opt.nb_threads - 1;

		queue<pair<string,vector<searchResult>>> q;

		//two possibilities: if number of query files given > number of threads, run each query file in a different thread
		//otherwise: distribute threads to query files, split given query file over assigned threads evenly

		if(opt.queryfiles.size() > 1 && num_threads > 0){

			int i = 0;
			while (i < opt.queryfiles.size()-1){
				//run the first t threads, then join them and start the next!
				std::vector<thread> threadList;
				int thread_counter = 0;
				while (thread_counter <= num_threads){
					if (i < opt.queryfiles.size()){
						threadList.push_back(std::thread(run_subsample_completeQuery, std::ref(q), std::ref(db_size), std::ref(k), opt.queryfiles[i], std::ref(tra), std::ref(opt.outprefix)));
						thread_counter++;
						i++;
					} else {
						break;
					}
				}

				//write the results!
				writeResults(opt.outprefix,q,cdbg);

				//Join the threads with the main thread
				for (auto& thread : threadList) {
					thread.join();
				}


			}
		} else {
			//split files and stuff!



			if (num_threads > 0){
				//parse fasta file
				string queryfile = opt.queryfiles[0];

				std::string delim = "/";
				auto start = 0U;
				auto end = queryfile.find(delim);
				while (end != std::string::npos){
					start = end + delim.length();
					end = queryfile.find(delim, start);
				}
				string query = queryfile.substr(start, end);

				//delete potentially existing search file for this query
				string f = opt.outprefix+"_"+query+".search";
				if (std::remove(f.c_str()) != 0) {
					perror( "Error deleting file" );
				} else {
					cout << "File " << f << " removed" << endl;
				}

				vector<pair<string,string>> fasta = parseFasta(queryfile);
				size_t bucketSize = fasta.size() / (num_threads);
				size_t leftOvers = fasta.size() % (num_threads);

				int thread_counter = 0;
				std::vector<thread> threadList;

				for(int i = 0; i <= fasta.size(); i=i+bucketSize) {
					if (thread_counter+1 < num_threads){
						threadList.push_back(std::thread(run_subsample_partialQuery, std::ref(q), std::ref(db_size), std::ref(k), std::ref(fasta),  std::ref(tra), i, (i+bucketSize-1),query));
						thread_counter++;
					} else {
						threadList.push_back(std::thread(run_subsample_partialQuery, std::ref(q), std::ref(db_size), std::ref(k), std::ref(fasta),  std::ref(tra),i,0,query));
						break;
					}
				}

				//write the results!
				cout << "write results!" << endl;
				writeResults(opt.outprefix,q,cdbg);

				//Join the threads with the main thread
				for (auto& thread : threadList) {
					thread.join();
				}
			}
			else {
				//we have only one thread to run, but could hav multiple query files
				for(auto& queryfile : opt.queryfiles){
					run_subsample_completeQuery(q, db_size, k, queryfile, tra, opt.outprefix);
					writeResults(opt.outprefix,q,cdbg);
				}
			}
		}



	cout << "Search took " << (float( clock () - begin_time ) /  CLOCKS_PER_SEC) << "sec." << endl;
	cout << "Goodbye!" << endl;
	}
}










