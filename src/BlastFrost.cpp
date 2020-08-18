#include <string>
#include <fstream>
#include <cstdio>
#include <thread>
#include <queue>
#include <mutex>
#include <atomic>
#include <unordered_map>
#include <sys/stat.h>

#include "ColoredCDBG.hpp"
#include "SubGraphTraverser.hpp"
#include "QuerySearch.hpp"
#include "BubbleExplorer.hpp"
#include "UnitigData.hpp"


using namespace std;



std::mutex mtx_queue;           // mutex for critical section
std::mutex mtx_writing;
std::mutex mtx_deleting;
std::atomic<int> running(0);



void PrintUsage() {
	cout << "BlastFrost -g <BifrostGraph> -f <BifrostColors> -q <query_sequences> -o <outfile_prefix> " << endl << endl;

	cout << "Mandatory parameters with required argument:" << endl << endl
			<< "  -g,         Input Bifrost graph file (GFA format)"
			<< endl
			<< "  -f,         Input Bifrost color file (BFG_COLORS format)"
			<< endl
			<< "  -q,         Query sequences (multiple FASTA)"
			<< endl
			<< "  -o,         Prefix for search result files (default: 'output')"

			<< endl << endl << "Optional parameters:"<< endl << endl
			<< "  -e,             Enable subgraph extraction"
			<< endl
			<< "  -t,             Number of threads (default is 1)"
			<< endl
			<< "  -k,             Length of k-mers (default is 31, max. is 31)"
			<< endl
			<< "  -d,           Compute and search d-neighborhood of queried k-mers"
			<< endl
			<< "  -c,           Enhance gfa file with color information, option can be run without a query search."
			<< endl
			<< "  -s,			Average size of genomes in Bifrost graph in Mb"
			<< endl
			<< "  -a,             Return file containing all colors currently present in Bifrost graph"
			<< endl
			<< "  -v,             Print information messages during construction"
			<< endl << endl;
}


struct searchOptions {
	string graphfile;
	string colorfile;
	vector<string> queryfiles;
	size_t nb_threads;
	bool verbose;
	string outprefix;
	int k;
	int d;
	bool enhanceGFA;
	double avg;
	bool extractSubgraph;
	string s1;
	string s2;
	bool returnColors;

	searchOptions() : nb_threads(1), verbose(false), outprefix("output"), k(31), d(0), enhanceGFA(false), avg(0), extractSubgraph(false), returnColors(false) {};
} ;



struct searchResult {
  string queryfile;
  long double evalue;
  long double pvalue;
  string query;
  size_t colorID;
  vector<int> hitrun;
} ;

inline bool exists_test (const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

void parseArgumentsNew(int argc, char **argv, searchOptions& opt) {

	int oc; //option character

	while ((oc = getopt(argc, argv, "o:t:q:g:f:k:d:l:r:s:evca")) != -1) {
		switch (oc) {
		//mandatory arguments
		case 'g':
			opt.graphfile = optarg;
			break;
		case 'f':
			opt.colorfile = optarg;
			break;
		case 'q':
			for (optind--; (optind < argc) && (*argv[optind] != '-'); ++optind) {
				opt.queryfiles.push_back(argv[optind]);
			}
			break;
		case 'o':
			/* handle -o, set fileprefix */
			opt.outprefix = optarg;
			break;

		//optional arguments
		case 't':
			/* handle -t, set number of threads */
			opt.nb_threads = atoi(optarg);
			break;
		case 'k':
			opt.k = atoi(optarg);
			break;
		case 'd':
			opt.d = atoi(optarg);
			break;
		case 'c':
			opt.enhanceGFA = true;
			break;
		case 'v':
			//handle -v, verbose mode
			opt.verbose = true;
			break;
		case 's':
			opt.avg = atof(optarg);
			break;
		//undocumented testing arguments
		case 'l':
			opt.s1 = optarg;
			break;
		case 'r':
			opt.s2 = optarg;
			break;
		case 'e':
			opt.extractSubgraph = true;
			break;
		case 'a':
			opt.returnColors = true;
			break;

		//invalid option
		case ':':
			cout << "Invalid option" << endl;
			PrintUsage();
			break;
		case '?':
			cout << "Invalid option" << endl;
			PrintUsage();
			break;
		default:
			cout << "Invalid option" << endl;
			PrintUsage();
			break;
		}
	}

	//check for mandatory arguments

	//check if precomputed graph prefix is given (-g)
	if (opt.graphfile.empty()){
		cout << "No input file given to load Bifrost graph!" << endl;
		PrintUsage();
		exit (EXIT_FAILURE);
	}
	else if (opt.colorfile.empty()){
                opt.colorfile = opt.graphfile.substr(0, opt.graphfile.length()-3)+"bfg_colors";
                if (access(opt.colorfile.c_str(), F_OK) == -1) {
                                cout << "No input file given to load Bifrost graph colors!" << endl;
		                PrintUsage();
                		exit (EXIT_FAILURE);

                }
	}
//	} else if (! exists_test(opt.queryfiles[0])){
//		cout << "Cannot read query file." << endl;
//		exit (EXIT_FAILURE);
//	}


}






void augmentGFA(SubGraphTraverser& tra, const string& gfaFile, const size_t& num_colors, bool verbose){
	// open gfa, read line by line, immadiately write to augmented file
	std::ofstream output(gfaFile+"_colored.gfa",std::ofstream::binary);
	if(verbose){
		cout << "### Initialized output file" << endl;
	}

	ifstream input(gfaFile);
	if(!input.good()){
		cout << "Cannot open gfa graph file." << endl;
		exit (EXIT_FAILURE);
	}

	unordered_map<string,int> colorHash;
	unordered_map<int,string> indexHash;
	int largest = 0;
	string line;
	while( std::getline(input,line).good() ){
		if( (! line.empty()) && line[0] == 'S'){
			//assumption: each S line has the same fields: S \t ID \t SEQ \t DZ-Tag
			vector<string> line_arr;
			string delim = "\t";
			//size_t pos = 0;

			auto start = 0U;
			auto end = line.find(delim);
			string sub;
			while (end != std::string::npos) {
				sub = line.substr(start, end - start);
				start = end + delim.length();
				end = line.find(delim, start);
				line_arr.push_back(sub);
			}

			sub = line.substr(start, end - start);
			line_arr.push_back(sub);

			vector<string> colors = tra.getColors(line_arr[2]);
			vector<bool> bset(num_colors, 0);
			for (auto& color: colors){
				if (colorHash.find(color) != colorHash.end()) {
					bset[colorHash[color]] = 1;
				} else {
					colorHash[color] = largest;
					indexHash[largest] = color;
					bset[largest] = 1;
					largest++;
				}
			}
			ostringstream os;
			copy(bset.begin(), bset.end(), ostream_iterator<bool>(os, ""));
			line_arr.push_back("CO:Z:"+os.str());
			output << line_arr[0] << "\t" << line_arr[1] << "\t" << line_arr[2] << "\t" << line_arr[3] << "\t" << line_arr[4] << endl; 			//write augmented line to output

		} else {
			output << line << endl; //write line directly to output
		}
	}

	input.close();
	output.close();

	//write the color index to separate file
	std::ofstream index(gfaFile+".index",std::ofstream::binary);
	std::map<int, string> ordered(indexHash.begin(), indexHash.end());
	for(auto it = ordered.begin(); it != ordered.end(); ++it){
	     index << it->second << " ";
	}
	index.close();
}




int estimate_avg_genomeSize(const vector<string>& files, bool verbose){
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
	if(verbose){
		//cout << "### Avg genome size: " << size << endl;
	}

	return size;
}





const string runLength_encode(const vector<int>& v){
	int cnt = 0;
	string res;
	string n;
	int lastElem = 3;
	for (auto& elem : v){
		if (elem == lastElem){
			cnt++;
		} else { //start a new run
			if (lastElem != 3){
				n = std::to_string(lastElem)+":"+std::to_string(cnt)+",";
				res+= n;
			}
			cnt = 1;
			lastElem = elem;
		}
	}

	n = std::to_string(lastElem)+":"+std::to_string(cnt)+",";
	res += n;
	return res;
}




vector<pair<string,string>> parseFasta(const string& queryfile, const unsigned int& k, bool verbose){
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
	int counter = 0;

	while( std::getline(input,line).good()) {
		if(line.empty() || line[0] == '>'){
			if(!name.empty()){
				if(content.length() >= k){
					fasta.push_back(std::make_pair(name,content));
				} else {
					counter += 1;
				}
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
		content += line;
		fasta.push_back(std::make_pair(name,content));
	}

	input.close();
	if(verbose){
		cout << "### Number of sequences that will not be queried due to length smaller than k: " << counter << endl;
	}
	return fasta;
}






void run_subsample_partialQuery(queue<pair<string,vector<searchResult>>>& q, const double& db_size, const int& k, const int& d, const vector<pair<string,string>>& fasta, QuerySearch& que, int start, int end, const string query, bool verbose){
		running++;
		//const int sigma = 4;

		if (end == 0){
			end = fasta.size();
		}

		vector<searchResult> results;
		for(int i = start; i < end; i++) {
			pair<string,string> seq = fasta[i];

			unordered_map<size_t,vector<int>> res = que.search(seq.second, k, d);
			//tra.remove_singletonHits(res);

			for (auto& hit : res){
				const int score = que.compute_score(hit.second);
				const int len = hit.second.size();

				long double pvalue2 = 0;
				long double evalue2 = 0;
				if (score != len){
					const long double log_evalue = que.compute_log_evalue(score,db_size,len);
					const long double log_pvalue = que.compute_log_pvalue(log_evalue);
					evalue2 = pow(10,log_evalue);
					pvalue2 = pow(10,log_pvalue);
				}



				if (pvalue2 <= 0.05){ //create new searchResult and push to all results
					searchResult result;
					result.evalue = evalue2;
					result.pvalue = pvalue2;
					result.query = seq.first;
					result.colorID = hit.first;
					result.hitrun = hit.second;
					results.push_back(result);

				} else {
					if(verbose){
						cout << "query hit not reported " << seq.first << " " << hit.first << " " << pvalue2 << endl;
					}
				}
			}

			if (results.size() > 30000){
				mtx_queue.lock();
				q.push(std::make_pair(query,results));
				mtx_queue.unlock();
				results.clear();
			}
		}

		//we've found all results for this query, push to writing queue!
		mtx_queue.lock();
		q.push(std::make_pair(query,results));
		mtx_queue.unlock();

	running--;
}


void run_subsample_completeQuery(queue<pair<string,vector<searchResult>>>& q, const double& db_size, const int& k, const int& d, const string queryfile, QuerySearch& que, string& outprefix, bool verbose){
	running++;
	//const int sigma = 4;

	vector<pair<string,string>> fasta = parseFasta(queryfile, k, verbose);

	std::string delim = "/";
	auto start = 0U;
	auto end = queryfile.find(delim);
	while (end != std::string::npos){
		start = end + delim.length();
		end = queryfile.find(delim, start);
	}

	string query = queryfile.substr(start, end);

	mtx_deleting.lock();
	//delete potentially existing search file for this query
	string f = outprefix+"_"+query+".search";
	if (access(f.c_str(), F_OK) != -1) {
		if (std::remove(f.c_str()) != 0) {
			perror("Error deleting file.");
		} else {
			if (verbose){
				cout << "### File " << f << " removed." << endl;
			}
		}
	}
	mtx_deleting.unlock();

	vector<searchResult> results;
	for(auto& seq: fasta) {
		unordered_map<size_t,vector<int>> res = que.search(seq.second, k, d);
		//tra.remove_singletonHits(res);

		for (auto& hit : res){
			const int score = que.compute_score(hit.second);
			const int len = hit.second.size();
			const long double log_evalue = que.compute_log_evalue(score,db_size,len);
			const long double log_pvalue = que.compute_log_pvalue(log_evalue);
			const long double evalue2 = pow(10,log_evalue);
			const long double pvalue2 = pow(10,log_pvalue);

			if (pvalue2 <= 0.05){ //create new searchResult and push to all results
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



void writeResults_singleFile(string& outprefix, string& query, queue<pair<string,vector<searchResult>>>& q, ColoredCDBG<UnitigData>& cdbg){
	//the main thread will take care of writing the output as soon as it is available
	mtx_writing.lock();
	std::ofstream output(outprefix+"_"+query+".search",std::ios_base::app);

 	while(running > 0){//there are still threads running
		if (! q.empty()){
			mtx_queue.lock();
			pair<string,vector<searchResult>> item = q.front();
			q.pop();
			mtx_queue.unlock();

			for (auto& result : item.second){
				string color = cdbg.getColorName(result.colorID);
				output << result.query << "\t" << color << "\t" << result.pvalue << "\t";
				string res = runLength_encode(result.hitrun);
				output << res << endl;
			}
		}
	}

 	while(! q.empty()){ //there is still data to write, but we do not need to lock anymore...
 		pair<string,vector<searchResult>> item = q.front();
 		q.pop();
 		for (auto& result : item.second){
 			string color = cdbg.getColorName(result.colorID);
 			output << result.query << "\t" << color << "\t" << result.pvalue << "\t";
 			string res = runLength_encode(result.hitrun);
 			output << res << endl;
 		}
 	}
 	output.close();
	mtx_writing.unlock();
}


void writeResults_multipleFiles(string& outprefix, queue<pair<string,vector<searchResult>>>& q, ColoredCDBG<UnitigData>& cdbg){
	//the main thread will take care of writing the output as soon as it is available
	mtx_writing.lock();

 	while(running > 0){//there are still threads running
		if (! q.empty()){
			mtx_queue.lock();
			pair<string,vector<searchResult>> item = q.front();
			q.pop();
			mtx_queue.unlock();

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

 	while(! q.empty()){ //there is still data to write, but we do not need to lock anymore...
 		pair<string,vector<searchResult>> item = q.front();
 		q.pop();
 		std::ofstream output(outprefix+"_"+item.first+".search",std::ios_base::app);
 		for (auto& result : item.second){
 			string color = cdbg.getColorName(result.colorID);
 			output << result.query << "\t" << color << "\t" << result.pvalue << "\t";
 			string res = runLength_encode(result.hitrun);
 			output << res << endl;
 		}
 		output.close();
 	}

	mtx_writing.unlock();
}




int main(int argc, char **argv) {

	if (argc < 2) {
		PrintUsage();
	} else {

	  //const clock_t load_time = clock();

		searchOptions opt;
		parseArgumentsNew(argc, argv, opt);


		//init and load pre-computed Bifrost graph
		ColoredCDBG<UnitigData> cdbg;
		if (cdbg.read(opt.graphfile, opt.colorfile, opt.nb_threads, opt.verbose)){
			cout << "Graph loading successful" << endl;
		} else {
			cout << "Graph could not be loaded! Exit." << endl;
			exit (EXIT_FAILURE);
		}


		//cout << "Loading took about" << ((float( clock() - load_time ) / (int)opt.nb_threads) / CLOCKS_PER_SEC) << "sec." << endl;
		////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//const clock_t begin_time = clock();


		//Average genome size
		int avg_genomeSize = 0;
		if(opt.avg == 0){
			avg_genomeSize = cdbg.size(); //overestimate?
			if(opt.verbose){
			  //cout << "Estimate for avg. genome size: " << avg_genomeSize << endl;
			}
		} else {
			avg_genomeSize = opt.avg * 1000000;
			if(opt.verbose){
			  //cout << "Input for avg. genome size: " << avg_genomeSize << endl;
			}
		}

		const size_t numberStrains = cdbg.getNbColors();
		double db_size = numberStrains*avg_genomeSize;


		//check which options are given by the user!
		if (opt.enhanceGFA){

			QuerySearch que(cdbg);
			SubGraphTraverser tra(cdbg, db_size, que);

			if(opt.verbose){
				cout << "--- Augment GFA file with color infos ---" << endl;
			}
			augmentGFA(tra,opt.graphfile,cdbg.getNbColors(), opt.verbose);

		//subgraph extraction
		} else if (opt.extractSubgraph){
			QuerySearch que(cdbg);
			SubGraphTraverser tra(cdbg, db_size, que);


			if(opt.verbose){
				cout << "run subgraph extraction" << endl;
			}

			for(auto& queryfile : opt.queryfiles){
				vector<pair<string,string>> fasta = parseFasta(queryfile, opt.k, opt.verbose);


				std::string delim = "/";
				auto start = 0U;
				auto end = queryfile.find(delim);
				while (end != std::string::npos){
					start = end + delim.length();
					end = queryfile.find(delim, start);
				}
				string query = queryfile.substr(start, end);

				//delete potentially existing search file for this query
				string f = opt.outprefix+"_"+query+"_subgraph.fasta";
                                if (access(f.c_str(), F_OK) != -1) {
					if (std::remove(f.c_str()) != 0) {
						perror( "Error deleting file" );
					} else {
						cout << "File " << f << " removed" << endl;
					}
				}

				int c = 0;
				for(auto& seq: fasta) {
					if (opt.verbose){
						cout << "Processing..." << c << endl;
						++c;
					}

					SubGraphTraverser::subgraphs result = tra.extractSubGraph_intelligent(seq.second, opt.k, opt.d);
					tra.printPaths_intelligent(opt.outprefix, query, result, seq.first);
				}
			}

			//cout << "Subgraph extraction took " << (float( clock() - begin_time ) /  CLOCKS_PER_SEC) << "sec." << endl;


		//TESTING/ IN DEVELOPMENT
		} else if (opt.s1 != "") {
			if(opt.verbose){
				cout << "---Explore subgraph between anchor sequences---" << endl;
			}

			BubbleExplorer bub(cdbg);

			vector<pair<string,string>> seq1 = parseFasta(opt.s1, opt.k, opt.verbose);
			//cout << seq1.size() << endl;
			//tra.exploreSubgraph(seq1[0].second);
			//cout << "---------------------------" << endl;

			vector<pair<string,string>> seq2 = parseFasta(opt.s2, opt.k, opt.verbose);
			//cout << seq2.size() << endl;
			//cout << seq2[0].second << endl;
			//tra.exploreSubgraph(seq2[0].second);

			//assume last k-mer of first and first k-mer of last as start and end unitig
			bub.exploreBubble(seq1[0].second, seq2[0].second, 1000);


		//return all colors currently present in the graph
		} else if (opt.returnColors){

			if(opt.verbose){
				cout << "return colors currently present in given graph" << endl;
			}

			vector<string> colors = cdbg.getColorNames();

			std::ofstream output(opt.graphfile+"_current_colors",std::ofstream::out);
			for(auto& col: colors){
				output << col << endl;
			}

			output.close();


			//Standard graph query
		} else {

			QuerySearch que(cdbg);

			if(opt.verbose){
				cout << "---Query graph for input sequences---" << endl;
			}


			int k = cdbg.getK(); //overwrite parameter given by user
			size_t num_threads = opt.nb_threads - 1; //!!! one thread is used by main!!!
			queue<pair<string,vector<searchResult>>> q;


			//two possibilities: if number of query files given > number of threads, run each query file in a different thread
			//otherwise: distribute threads to query files, split given query file over assigned threads evenly
			if(opt.queryfiles.size() > 1 && num_threads > 0){
				unsigned int i = 0;
				while (i < opt.queryfiles.size()-1){

					std::vector<thread> threadList;
					unsigned int thread_counter = 0;
					while (thread_counter <= num_threads){
						if (i < opt.queryfiles.size()){
							threadList.push_back(std::thread(run_subsample_completeQuery, std::ref(q), std::ref(db_size), std::ref(k), std::ref(opt.d), opt.queryfiles[i], std::ref(que), std::ref(opt.outprefix), opt.verbose));
							thread_counter++;
							i++;
						} else {
							break;
						}
					}

					writeResults_multipleFiles(opt.outprefix,q,cdbg); //write the results!

					for (auto& thread : threadList) { //Join the threads with the main thread
						thread.join();
					}
				}
			} else { //split files and stuff!
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
				        if (access(f.c_str(), F_OK) != -1) {
						if (std::remove(f.c_str()) != 0) {
                       					perror("Error deleting file.");
                				} else {
                                			cout << "### File " << f << " removed." << endl;
                				}
       				 	}

					vector<pair<string,string>> fasta = parseFasta(queryfile, opt.k, opt.verbose);
					size_t bucketSize = 0;
					if (fasta.size() < num_threads ){
						bucketSize = fasta.size();
					} else {
						bucketSize = fasta.size() / (num_threads);
						//size_t leftOvers = fasta.size() % (num_threads);
					}

					unsigned int thread_counter = 0;
					std::vector<thread> threadList;

					for(unsigned int i = 0; i < fasta.size(); i=i+bucketSize) {

						if (thread_counter+1 < num_threads){
							threadList.push_back(std::thread(run_subsample_partialQuery, std::ref(q), std::ref(db_size), std::ref(k), std::ref(opt.d), std::ref(fasta),  std::ref(que), i, (i+bucketSize),query, opt.verbose));
							thread_counter++;
						} else {
							threadList.push_back(std::thread(run_subsample_partialQuery, std::ref(q), std::ref(db_size), std::ref(k), std::ref(opt.d), std::ref(fasta),  std::ref(que),i,0,query, opt.verbose));
							break;
						}
					}

					//write the results!
					if(opt.verbose){
						cout << "Start writing..." << endl;
					}
					writeResults_singleFile(opt.outprefix,query,q,cdbg);

					//Join the threads with the main thread
					for (auto& thread : threadList) {
						thread.join();
					}
				}
				else { //we have only one thread to run, but could have multiple query files
					for(auto& queryfile : opt.queryfiles){
						run_subsample_completeQuery(q, db_size, k, opt.d, queryfile, que, opt.outprefix, opt.verbose);
						writeResults_multipleFiles(opt.outprefix,q,cdbg);

					}
				}
			}

			//cout << "Search took " << (float( clock () - begin_time ) /  CLOCKS_PER_SEC) << "sec." << endl;
		}
		cout << "Goodbye!" << endl;
	}
}










