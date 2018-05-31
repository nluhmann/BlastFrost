
#include <string>
#include <fstream>
#include <cstdio>
#include <array>

#include <bifrost/ColoredCDBG.hpp>
#include "GraphTraverser.hpp"
#include "UnitigData.hpp"



using namespace std;

void PrintUsage() {
	cout << "BlastFrost -f <inputFiles> -o <outfile_graph> -q <query_sequences> -r <query_result>" << endl << endl;

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

void parseArguments(int argc, char **argv, CCDBG_Build_opt& opt, string& queryfile, string& resultfile) {

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
			resultfile = optarg;
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


}



void writeOutput(unordered_map<string,vector<int>>& arr, string& outfile){
	std::ofstream output(outfile,std::ofstream::binary);

	for (auto& elem : arr){
		output << elem.first;
		for (auto& p : elem.second){
			output << "\t" << p;
		}
		output << endl;
	}


	output.close();
}


void transformOutput(vector<pair<Kmer,set<string>>> res){
	//traverse through vector once to get set of all colors identified at least once
	set<string> colors;
	int numberKmer = 0;
	for (auto& p : res){
		numberKmer++;
		for (auto& color : p.second){
			colors.insert(color);
		}
	}

	//next, for each strain in colors, create an int[] with length = numberKmers, initialized to 0
	unordered_map<string,vector<int>> arr;
	for (auto& color: colors){
		vector<int> vec(numberKmer, 0);
		arr.insert({color,vec});
	}

	//traverse through vector res once more, set 1 for all colors for which the Kmer has been found
	int counter = 0;
	for (auto& p : res){
		for (auto& color : p.second){
			arr[color][counter] = 1;
		}
		counter++;
	}

	string outfile = "out.matrix";
	writeOutput(arr,outfile);
}





//ToDo: only one query supported atm


int main(int argc, char **argv) {

	if (argc < 2) {
		PrintUsage();
	} else {

		//Parse input arguments, parsing graph parameters and file specifying genome clustering
		CCDBG_Build_opt opt;
		string queryfile;
		string resultsfile;
		parseArguments(argc, argv, opt, queryfile, resultsfile);


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
		if (std::remove(resultsfile.c_str()) != 0 ) {
		    perror( "Error deleting file" );
		}

		//create and check input stream
		ifstream input(queryfile);
		if(!input) {
			std::cout << "Cannot open input file." << std::endl;
			return 1;
		}

		//ToDo: this is not robust to line breaks in fasta sequence
		//read file line by line, parse header and sequence, run search as soon as both are present
		bool header_bool = true;
		string header;
		string line;
		while (getline(input, line)) {
			if(header_bool){
				header = line;
				header_bool = false;
			} else {
				//we are reading a seq
				vector<pair<Kmer,set<string>>> res = tra.search(line, opt.k);
				tra.writeKmerPresence(res, resultsfile);
				transformOutput(res);
				header_bool = true;
			}


		//output file 1: for each kmer in query, report strains present
		//output file 2: for each strain found, report all kmer present

		}

		input.close();


	cout << "Search took " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << "sec." << endl;
	cout << "Goodbye!" << endl;
	}
}


