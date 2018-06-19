
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



std::mutex mtx;           // mutex for critical section
std::atomic<int> running(0);



void PrintUsage() {
	cout << "BlastFrost -f <inputFiles> -q <query_sequences> -o <outfile_prefix> " << endl << endl;

	cout << "Mandatory parameters with required argument:" << endl << endl
			<< "  -f,         Input sequence files (FASTA or FASTQ, possibly gziped) and/or graph files (GFA)"
			<< endl
			<< "  -o,         Prefix for output file (GFA output by default)"
			<< endl << endl << "Optional parameters with required argument:"
			<< endl << endl
			<< "  -t,             Number of threads (default is 1)"
			<< endl << endl << "Optional parameters with no argument:" << endl
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


vector<pair<string,string>> parseFasta(string& queryfile){

	//create and check input stream
	ifstream input(queryfile);
	if(!input.good()){
		cout << "Cannot open query file." << endl;
		exit (EXIT_FAILURE);
	}


	vector<pair<string,string>> fasta;

	string line, name, content;
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



void run_subsample(queue<unordered_map<size_t,vector<int>>>& q, vector<pair<string,string>>& fasta, GraphTraverser& tra, int start, int end){
	running++;


	//mtx.lock();
	//cout << "Start thread!" << endl;
	//cout << start << endl;
	//cout << end << endl;


	long double k = 31;
	long double x = 4500000;
	long double sigma = 4;

	long double p = tra.compute_p(k, x, sigma);


	queue<unordered_map<size_t,vector<int>>> local_queue;
	for(int i = start; i <= end; i++) {

		pair<string,string> seq = fasta[i];
		unordered_map<size_t,vector<int>> res = tra.search(seq.second, 31);
		tra.remove_singletonHits(res, 31);

		unordered_map<size_t,vector<int>> filtered;
		for (auto& hit : res){
			int score = tra.compute_score(hit.second);


			int len = hit.second.size();
			double db_size = 714*x;

			long double log_evalue = tra.compute_log_evalue(score,db_size,len);


			long double log_pvalue = tra.compute_log_pvalue(log_evalue);


			long double evalue2 = pow(10,log_evalue);


			long double pvalue2 = pow(10,log_pvalue);

			//if (pvalue2 > 0.05) {
			//	cout << seq.first << endl;
			//	cout << "score: " << score << endl;
			//	cout << "log e-value: " << log_evalue << endl;
			//	cout << "log p-value: " << log_pvalue << endl;
			//	cout << "e-value2: " << evalue2 << endl;
			//	cout << "p-value2: " << pvalue2 << endl;
			//}

			if (pvalue2 <= 0.05){
				filtered.insert(hit);
			}

		}

		//unordered_map<size_t,double> p_values = tra.compute_significance(res,p);


		string out = "search_"+seq.first;
		tra.writePresenceMatrix(filtered,out);




		//local_queue.push(res);


		//if(local_queue.size() > 40){
			//lock general queue
			//mtx.lock();
			//while (!local_queue.empty()){
				//unordered_map<size_t,vector<int>> item = local_queue.front();
				//local_queue.pop();
				//q.push(item);
			//}
			//unlock queue
			//mtx.unlock();
		//}
	}
	//mtx.unlock();
	running--;
}


void test_stats(){
	long double k = 31;
	long double x = 3000000000;
	long double sigma = 4;
	long double y = pow(sigma,-k);
	cout << "y: " << y << endl;
	long double a = pow(1-y,x);
	cout << "a: " << a << endl;
	long double p = 1-pow(1-y,x);
	cout << "p: " << p << endl;


	double b = gsl_ran_binomial_pdf(2, 0.5, 20);
	cout << b << endl;
	b = gsl_ran_binomial_pdf(5, 0.5, 20);
	cout << b << endl;
	b = gsl_ran_binomial_pdf(8, 0.5, 20);
	cout << b << endl;
	b = gsl_ran_binomial_pdf(12, 0.5, 20);
	cout << b << endl;
	b = gsl_ran_binomial_pdf(18, 0.5, 20);
	cout << b << endl;

	cout << "------" << endl;

	double r = gsl_cdf_binomial_Q(2, 0.5, 20);
	cout << r << endl;
	r = gsl_cdf_binomial_Q(5, 0.5, 20);
	cout << r << endl;
	r = gsl_cdf_binomial_Q(8, 0.5, 20);
	cout << r << endl;
	r = gsl_cdf_binomial_Q(12, 0.5, 20);
	cout << r << endl;
	r = gsl_cdf_binomial_Q(13, 0.5, 20);
	cout << r << endl;
	r = gsl_cdf_binomial_Q(14, 0.5, 20);
	cout << r << endl;
	r = gsl_cdf_binomial_Q(150, 0.5, 200);
	cout << r << endl;
	r = gsl_cdf_binomial_Q(16, 0.5, 20);
	cout << r << endl;
	r = gsl_cdf_binomial_Q(17, 0.5, 20);
	cout << r << endl;
	r = gsl_cdf_binomial_Q(18, 0.5, 20);
	cout << r << endl;

	cout << "------" << endl;

	double res = gsl_cdf_binomial_Q(1, p, 1000);
	cout << res << endl;
	res = gsl_cdf_binomial_Q(2, p, 1000);
	cout << res << endl;
	res = gsl_cdf_binomial_Q(5, p, 1000);
	cout << res << endl;
	res = gsl_cdf_binomial_Q(500, p, 1000);
	cout << res << endl;
	res = gsl_cdf_binomial_Q(1000, p, 1000);
	cout << res << endl;
}


int main(int argc, char **argv) {

	if (argc < 2) {
		PrintUsage();
	} else {



		//Parse input arguments, parsing graph parameters and file specifying genome clustering
		CCDBG_Build_opt opt;
		string queryfile;
		parseArguments(argc, argv, opt, queryfile);

		//string resultsfile = opt.prefixFilenameOut+".query";

		//build the graph, done for all command variations
		cout << "---Build de Bruijn graph with BiFrost---" << endl;

		//We have assemblies as input, so need the reference mode (counting singleton kmers) and want to have colors
		opt.reference_mode = true;
		opt.outputColors = true;

		//create cdbg instance, k=kmer length (default 31), g=length of minimizer (default 23)
		ColoredCDBG<UnitigData> cdbg(opt.k);
		cdbg.build(opt);
		cdbg.simplify(opt.deleteIsolated, opt.clipTips, opt.verbose);
		cdbg.mapColors(opt);
		cdbg.write(opt.prefixFilenameOut, opt.nb_threads, opt.verbose);



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
		vector<pair<string,string>> fasta = parseFasta(queryfile);


		queue<unordered_map<size_t,vector<int>>> q;

		//!!! one thread is used by main!!!
		size_t& num_threads = opt.nb_threads;

		size_t bucketSize = fasta.size() / (num_threads-1);
		cout << "Bucket size: " << bucketSize << endl;
		size_t leftOvers = fasta.size() % (num_threads-1);
		cout << "Leftover: " << leftOvers << endl;


		int last = 0;

		int thread_counter = 0;
		std::vector<thread> threadList;



		for(int i = 0; i <= fasta.size(); i=i+bucketSize) {
			if (thread_counter+1 < num_threads-1){
				//cout << i;
				//cout << " ";
				//cout << i+bucketSize-1 << endl;

				threadList.push_back( std::thread(run_subsample, std::ref(q), std::ref(fasta), std::ref(tra), i, (i+bucketSize-1)));
				last = i+bucketSize;
				thread_counter++;
			} else {
				threadList.push_back( std::thread(run_subsample, std::ref(q), std::ref(fasta), std::ref(tra),last,last+bucketSize+leftOvers-1));
				//cout << last;
				//cout << " ";
				//cout << last+bucketSize+leftOvers-1 << endl;
				break;
			}
		}


		//todo: what if there is only 1 thread???




		//the main thread will take care of writing the output as soon as it is available
//		std::ofstream output(opt.prefixFilenameOut+".search",std::ofstream::binary);
//
//		while(running > 0){//there are still threads running
//			while (!q.empty()){
//				mtx.lock();
//				unordered_map<size_t,vector<int>> item = q.front();
//				q.pop();
//				mtx.unlock();
//				for (auto& elem : item){
//					string color = cdbg.getColorName(elem.first);
//					output << color;
//					for (auto& p : elem.second){
//						output << "\t" << p;
//					}
//					output << endl;
//				}
//			}
//		}
//
//
//		output.close();

		//Join the threads with the main thread
		for (auto& thread : threadList) {
			thread.join();
		}


	cout << "Search took " << (float( clock () - begin_time ) /  CLOCKS_PER_SEC) / 100 << "sec." << endl;
	cout << "Goodbye!" << endl;
	}
}










