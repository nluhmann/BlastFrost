
#include <string>
#include <fstream>
#include <cstdio>

#include <bifrost/ColoredCDBG.hpp>
#include "GraphTraverser.hpp"
#include "UnitigData.hpp"



using namespace std;

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






//ToDo: only one query supported atm


int main(int argc, char **argv) {

	if (argc < 2) {
		PrintUsage();
	} else {

		//Parse input arguments, parsing graph parameters and file specifying genome clustering
		CCDBG_Build_opt opt;
		string queryfile;
		parseArguments(argc, argv, opt, queryfile);

		string resultsfile = opt.prefixFilenameOut+".query";

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
				unordered_map<size_t,vector<int>> res = tra.search(line, opt.k);
				tra.writePresenceMatrix(res, resultsfile);

				header_bool = true;
			}


		}

		input.close();


	cout << "Search took " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << "sec." << endl;
	cout << "Goodbye!" << endl;
	}
}


