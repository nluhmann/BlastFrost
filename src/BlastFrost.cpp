
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

void parseArguments(int argc, char **argv, CCDBG_Build_opt& opt, string& queryfile) {

	int oc; //option character

	while ((oc = getopt(argc, argv, "f:g:k:o:t:q:r:iv")) != -1) {
		switch (oc) {
		case 'f':
			/* handle -f, collect input files */
			for (optind--; (optind < argc) && (*argv[optind] != '-');
					++optind) {
			  //cout << argv[optind] << endl;
				opt.filename_seq_in.push_back(argv[optind]);
			}
			break;
		case 'o':
			/* handle -o, set fileprefix */
			opt.prefixFilenameOut = optarg;
			break;
		case 't':
			/* handle -t, set number of threads */
			opt.nb_threads = atoi(optarg);
			break;
		case 'k':
			opt.k = atoi(optarg);
			break;
		case 'g':
			opt.g = atoi(optarg);
			break;
		case 'i':
			// handle -i, tip clipping
			opt.clipTips = true;
			break;
		case 'v':
			//handle -v, verbose mode
			opt.verbose = true;
			break;
		case 'q':
			queryfile = optarg;
			break;
		case 'r':
			queryfile = optarg;
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


	//check if input sequences given (-f)
	if (opt.filename_seq_in.size() == 0){
		cout << "No input files given to build graph!" << endl;
		 exit (EXIT_FAILURE);
	}

	//check if output prefix is given (-o), otherwise set default to "output"
	if (opt.prefixFilenameOut.empty()){
		cout << "No outfile prefix given, set default value 'output'" << endl;
		opt.prefixFilenameOut = "output";
	}

	//check if query file is given (-q)
	if (queryfile.empty()){
		cout << "No query file given to search graph!" << endl;
		exit (EXIT_FAILURE);
	}
}

struct searchOptions {
  string graphfile = "";
  string queryfile = "";
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
		case 'q':
			opt.queryfile = optarg;
			break;
		case 'g':
			opt.graphfile = optarg;
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


	//check if precomputed graph prefx is given (-g)
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
	if (opt.queryfile == ""){
		cout << "No query file given to search graph!" << endl;
		exit (EXIT_FAILURE);
	}
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
  long double evalue;
  long double pvalue;
  string query;
  size_t colorID;
  vector<int> hitrun;
} ;




void run_subsample(queue<vector<searchResult>>& q, const ColoredCDBG<UnitigData>& cdbg, const vector<pair<string,string>>& fasta, GraphTraverser& tra, int start, int end, const int& avg_genomeSize){
	running++;
	//mtx.lock();

	//long double k = 31;
	//const int x = 4500000;
	const int sigma = 4;

	//long double p = tra.compute_p(k, x, sigma);
	const size_t numberStrains = cdbg.getNbColors();

	//ToDo: check db_size parameter
	//double db_size = 714*x;
	const double db_size = numberStrains*avg_genomeSize;

	if (end == 0){
		end = fasta.size()-1;
	}

	for(int i = start; i <= end; i++) {
		vector<searchResult> results;
		pair<string,string> seq = fasta[i];
		unordered_map<size_t,vector<int>> res = tra.search(seq.second, cdbg.getK());
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

			} else {
					//cout << seq.first << endl;
					//cout << "score: " << score << endl;
					//cout << "log e-value: " << log_evalue << endl;
					//cout << "log p-value: " << log_pvalue << endl;
					//cout << "e-value2: " << evalue2 << endl;
					//cout << "p-value2: " << pvalue2 << endl;
			}

		}

		//we've found all results for this query, push to writing queue!
		mtx_queue.lock();
		q.push(results);
		mtx_queue.unlock();

	}

	running--;
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
		//string resultsfile = opt.prefixFilenameOut+".query";

		//build the graph, done for all command variations
		//cout << "---Build de Bruijn graph with BiFrost---" << endl;

		//We have assemblies as input, so need the reference mode (counting singleton kmers) and want to have colors
		//opt.reference_mode = true;
		//opt.outputColors = true;

		//create cdbg instance, k=kmer length (default 31), g=length of minimizer (default 23)
		//ColoredCDBG<UnitigData> cdbg(opt.k);
		//cdbg.build(opt);
		//cdbg.simplify(opt.deleteIsolated, opt.clipTips, opt.verbose);
		//cdbg.mapColors(opt);
		//cdbg.write(opt.prefixFilenameOut, opt.nb_threads, opt.verbose);

		////////////////////////////////////////////////////////////////////////////////////////////////////////////

		const clock_t begin_time = clock();


		cout << "---Query graph for input sequences---" << endl;
		GraphTraverser tra(cdbg);

		/*
		 * for each query in multiple-fasta input file, search all k-mers for presence in the graph
		 * if not present: report 0
		 * if present: report all colors for this kmer
		 */


		//remove resultsfile if existing
		//if (std::remove(resultsfile.c_str()) != 0 ) {
		//    perror( "Error deleting file" );
		//}

		//parse fasta file
		vector<pair<string,string>> fasta = parseFasta(opt.queryfile);
		//const int avg_genomeSize = estimate_avg_genomeSize(opt.filename_seq_in);
		const int avg_genomeSize = cdbg.size();
		cout << avg_genomeSize << endl;


		//!!! one thread is used by main!!!
		size_t& num_threads = opt.nb_threads;

		queue<vector<searchResult>> q;

		cout << "Number threads: " << num_threads << endl;
		if (num_threads > 1){
			cout << "Parallel mode!" << endl;

			size_t bucketSize = fasta.size() / (num_threads-1);
			cout << "Bucket size: " << bucketSize << endl;
			size_t leftOvers = fasta.size() % (num_threads-1);
			cout << "Leftover: " << leftOvers << endl;

			int thread_counter = 0;
			std::vector<thread> threadList;

			for(int i = 0; i <= fasta.size(); i=i+bucketSize) {
				if (thread_counter+1 < num_threads-1){
					threadList.push_back( std::thread(run_subsample, std::ref(q), std::ref(cdbg), std::ref(fasta), std::ref(tra), i, (i+bucketSize-1), std::ref(avg_genomeSize)));
					thread_counter++;
				} else {
					threadList.push_back( std::thread(run_subsample, std::ref(q), std::ref(cdbg), std::ref(fasta), std::ref(tra),i,0,std::ref(avg_genomeSize)));
					break;
				}
			}





			//the main thread will take care of writing the output as soon as it is available
			std::ofstream output(opt.outprefix+".search",std::ofstream::binary);

			mtx_writing.lock();
			while(running > 0 || ! q.empty()){//there are still threads running
				if (! q.empty()){
					mtx_queue.lock();
					vector<searchResult> item = q.front();
					q.pop();
					mtx_queue.unlock();
					for (auto& result : item){
						string color = cdbg.getColorName(result.colorID);
						output << result.query << "\t" << color << "\t" << result.pvalue << "\t";
						string res = runLength_encode(result.hitrun);
						output << res << endl;
					}
				}
			}

			mtx_writing.unlock();
			output.close();


			//Join the threads with the main thread
			for (auto& thread : threadList) {
				thread.join();
			}


		}
		else {
			cout << "Serial mode!" << endl;
			run_subsample(q, cdbg, fasta, tra, 0, fasta.size()-1, avg_genomeSize);
			std::ofstream output(opt.outprefix+".search",std::ofstream::binary);
			while (!q.empty()){
				vector<searchResult> item = q.front();
				q.pop();
				for (auto& result : item){
					string color = cdbg.getColorName(result.colorID);
					output << result.query << "\t" << color << "\t" << result.pvalue << "\t";
					string res = runLength_encode(result.hitrun);
					output << res << endl;
				}
			}
			output.close();
		}


	cout << "Search took " << (float( clock () - begin_time ) /  CLOCKS_PER_SEC) << "sec." << endl;
	cout << "Goodbye!" << endl;
	}
}










