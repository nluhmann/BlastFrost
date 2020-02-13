/*
 * BubbleExplorer.cpp
 *
 *  Created on: 14 Jan 2019
 *      Author: nina
 */



#include "BubbleExplorer.hpp"

using namespace std;


BubbleExplorer::BubbleExplorer(ColoredCDBG<UnitigData>& graph) : cdbg(graph) {
	cout << "BubbleExplorer initialized!" << endl;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// FUNC 3: Explore bubble defined by two genes (left and right border)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/*
 * IN DEVELOPMENT
 * explore the bubble defined by the left and right gene by a depth-first search,
 * stopping each run if a max. threshold is met
 */
void BubbleExplorer::exploreBubble(const string& left, const string& right, const int threshold){

	const char *left_str = left.c_str();
	const char *right_str = right.c_str();

	UnitigColorMap<UnitigData> first;
	UnitigColorMap<UnitigData> last;


	//ToDo: cannot do that. first, run BlastFrost for left and right string with d>2,



	//for left, find last kmer that is found in the graph --- NOT NICE
	for(KmerIterator it_km(left_str), it_km_end; it_km != it_km_end; ++it_km){

		const UnitigColorMap<UnitigData> map = cdbg.find(it_km->first);

		if (!map.isEmpty){
			last = map;
		}
	}


	//for right, find first kmer that is found in the graph
	for(KmerIterator it_km(right_str), it_km_end; it_km != it_km_end; ++it_km){

		const UnitigColorMap<UnitigData> map = cdbg.find(it_km->first);

		if (!map.isEmpty){

			first = map;
			break;
		}
	}


	last.dist = 0;
	last.len = last.size - last.getGraph() -> getK() + 1;

	//first.dist = 0;
	//first.len = first.size - first.getGraph() -> getK() + 1;

	//search ALL paths between the last and first unitig, with a search threshold for the total length of paths
	vector<vector<UnitigColorMap<UnitigData>>> allPaths = BubbleExplorer::DFS_Iterative(last, first.getUnitigHead(), threshold);

	BubbleExplorer::printPaths(allPaths);
}




/*
 * Iterative Depth-first search.
 */
vector<vector<UnitigColorMap<UnitigData>>> BubbleExplorer::DFS_Iterative(const UnitigColorMap<UnitigData>& start, const Kmer& stop, const unsigned int threshold){
	//record the path we are iterating over, and push it to the result if we reach the stop Kmer

	stack<UnitigColorMap<UnitigData>> stck; // Create stack of unitig to traverse
	UnitigColorMap<UnitigData> ucm_tmp(start); // Create a non-const local copy of unitig given in parameter

	stck.push(ucm_tmp); // Push first unitig to traverse on the stack

	vector<vector<UnitigColorMap<UnitigData>>> allPaths;
	bool stop_reached = false;
	vector<UnitigColorMap<UnitigData>> currentPath;

	while (!stck.empty()){ // While there are unitigs to traverse in the stack

	    	ucm_tmp = stck.top(); // Get the unitig on top of the stack
	    	stck.pop(); // Delete unitig on the top of the stack

	    	if(std::find(currentPath.begin(), currentPath.end(), ucm_tmp) != currentPath.end()) {
	    	    stop_reached = true;
	    	} else {
	    		currentPath.push_back(ucm_tmp);
	    	}

	    	DataAccessor<UnitigData>* da = ucm_tmp.getData(); // Get DataAccessor from unitig
	    	UnitigData* data = da->getData(ucm_tmp); // Get boolean from DataAccessor

    		data->set_visited(); // Set boolean to indicate unitig was visited

	    	//CASE 1: we found the end
	    	if (ucm_tmp.getUnitigHead() == stop){
	    	    cout << "Found the end!" << endl;
	    	    stop_reached = true;
	    	    allPaths.push_back(currentPath);
	    	    for (auto& unitig : currentPath){
	    	    	DataAccessor<UnitigData>* pathda = unitig.getData(); // Get DataAccessor from unitig
	    	    	UnitigData* pathdata = pathda->getData(unitig); // Get boolean from DataAccessor
	    	    	pathdata->set_on_target_path();
	    	    }
	    	}

	    	//CASE 2: we reached the limit
	    	if (currentPath.size() > threshold){
	    		stop_reached = true;
	    	}

	    	if (stop_reached){ //backtrace
	    		stop_reached = false;

	    		//adjust current path such that only the predecessor of the next unitig on the stack is present
	    		UnitigColorMap<UnitigData> ucm_next = stck.top();

	    		//DataAccessor<UnitigData>* nextda = ucm_next.getData(); // Get DataAccessor from unitig
	    		//UnitigData* nextdata = nextda->getData(ucm_next); // Get boolean from DataAccessor

	    		bool broken = false;
	    		for (auto& predecessor : ucm_next.getPredecessors()){
	    			for (vector<UnitigColorMap<UnitigData>>::reverse_iterator it = currentPath.rbegin(); it != currentPath.rend(); ++it ) {
	    				if (predecessor == *(it.base())){
	    					//std::advance(it, 1);
	    					currentPath.erase((--it).base(),currentPath.end());
	    					broken = true;
	    					break;
	    				}
	    			}
	    		}
	    		if (! broken){
	    			currentPath.clear();
	    		}
	    	} else { //CASE 3: just go on

	    		bool succs = false;
	    		for (auto& successor : ucm_tmp.getSuccessors()) {
	    			DataAccessor<UnitigData>* succda = successor.getData(); // Get DataAccessor from unitig
	    			UnitigData* succdata = succda->getData(successor); // Get boolean from DataAccessor
	    			if(succdata->is_not_visited()){
	    				stck.push(successor);
	    				succs = true;
	    			} else if (succdata -> is_on_target_path()) { //we can stop here, cause the rest of the path has been found already previously!
	    				allPaths.push_back(currentPath);
	    			}

	    		}
	    		if (! succs){
	    			if (! stck.empty()){
	    				UnitigColorMap<UnitigData> ucm_next = stck.top();

	    				//DataAccessor<UnitigData>* nextda = ucm_next.getData(); // Get DataAccessor from unitig
	    				//UnitigData* nextdata = nextda->getData(ucm_next); // Get boolean from DataAccessor

	    				bool broken = false;
	    				for (auto& predecessor : ucm_next.getPredecessors()){
	    					for (vector<UnitigColorMap<UnitigData>>::reverse_iterator it = currentPath.rbegin(); it != currentPath.rend(); ++it ) {
	    						if (predecessor == *(it.base())){
	    							currentPath.erase((--it).base(),currentPath.end());
	    							broken = true;
	    							break;
	    						}
	    					}
	    				}
	    				if (! broken){
	    					currentPath.clear();
	    				}
	    			}
	    		}
	    	}
	}
	return allPaths;
}






void BubbleExplorer::printPaths(const vector<vector<UnitigColorMap<UnitigData>>> allPaths){
	int pathcounter = 0;

	std::ofstream p("paths",std::ofstream::binary);
	for (auto& path : allPaths){
		int unitigCounter = 0;

		for (auto& elem : path) {
			p << ">path" << pathcounter << "_" << unitigCounter << endl;
			p << elem.mappedSequenceToString() << endl;

			const DataAccessor<UnitigData>* da = elem.getData();
			const UnitigColors* uc = da->getUnitigColors(elem);

			p << "#";
			for (UnitigColors::const_iterator it = uc->begin(elem); it != uc->end(); it.nextColor()) {
				const size_t colorID = it.getColorID();
				const string color = cdbg.getColorName(colorID);
				p << color << ",";
			}
			p << endl;
			++unitigCounter;
		}
		++pathcounter;
	}
	p.close();
}

