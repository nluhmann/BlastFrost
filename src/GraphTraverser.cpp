/*
 * GraphTraverser.cpp
 *
 *  Created on: 26 Feb 2018
 *      Author: nina
 */

#include "GraphTraverser.hpp"

using namespace std;

/*
 * ToDo: Load gfa from 120k Salmonella, traverse graph with assemblies, store each assembly as sequence of unitigIDs
 * Load gfa from 120k Salmonella, traverse graph with all cgMLST loci/alleles, store each loci as a sequence of unitigIDs
 * Load gfa from 120k Salmonella, traverse graph with all wgMLST loci/alleles, store each of them as a sequence of unitigIDs
 * Stats: Load gfa from 120k Salmonella, traverse graph with all cgMLST/wgMLST and count the number of connected components
 * Question: how does recombination influence the graph?
 */




GraphTraverser::GraphTraverser(ColoredCDBG<UnitigData>& graph) :
	cdbg(graph) {
	cout << "GraphTraverser initialized!" << endl;
}



/*
 *
 */
int GraphTraverser::countConnectedComponents() {
	unordered_map<string, int> unitig_ids;
	int counter = 0;
	for (ColoredCDBG<UnitigData>::iterator it = cdbg.begin(); it != cdbg.end(); ++it) {
		Kmer km = GraphTraverser::getKmer(*it);
		unitig_ids[km.toString()] = counter;
		counter++;
	}

	UnionFind uf(counter);
	for (ColoredCDBG<UnitigData>::iterator it = cdbg.begin(); it != cdbg.end(); ++it) {
		Kmer left = GraphTraverser::getKmer(*it);

		//visit all successors of this unitig
		ForwardCDBG<DataAccessor<UnitigData> , DataStorage<UnitigData> , false>
				forward = it->getSuccessors();
		for (auto forw = forward.begin(); forw != forward.end(); ++forw) {
			Kmer right = GraphTraverser::getKmer(*forw);
			uf.merge(unitig_ids[left.toString()], unitig_ids[right.toString()]);
		}

		BackwardCDBG<DataAccessor<UnitigData> , DataStorage<UnitigData> , false>
				back = it->getPredecessors();
		for (auto backw = back.begin(); backw != back.end(); ++backw) {
			Kmer right = GraphTraverser::getKmer(*backw);
			uf.merge(unitig_ids[left.toString()], unitig_ids[right.toString()]);
		}
	}
	return uf.count();
}





/*
 * Count all the bubbles! Does not work yet. Debug submethod findSuperBubble().
 */
int GraphTraverser::countSuperBubbles() {
	int counter = 0;
	for (auto& unitig : cdbg) {
		vector<UnitigColorMap<UnitigData>> bubblePair = GraphTraverser::findSuperBubble(unitig);
		if (!bubblePair.empty()) {
			counter++;
		}
	}
	return counter;
}


//ToDo: debug!!!
vector<UnitigColorMap<UnitigData>> GraphTraverser::findSuperBubble(
		const UnitigColorMap<UnitigData>& unitig) {
	vector<UnitigColorMap<UnitigData>> S;
	vector<UnitigColorMap<UnitigData>> pair;
	int seen_but_unvisited = 0;

	cout << "Start" << endl;
	Kmer km = GraphTraverser::getKmer(unitig);
	S.push_back(unitig);
	while (!S.empty()) {
		cout << "Size of S: " << S.size() << endl;
		//get last element from S and remove it from vector
		UnitigColorMap<UnitigData> active = S.back();
		S.pop_back();

		DataAccessor<UnitigData>* da = active.getData(); // Get DataAccessor from unitig
		UnitigData* data = da->getData(active); // Get boolean from DataAccessor
		data->set_visited();
		seen_but_unvisited--;

		ForwardCDBG<DataAccessor<UnitigData> , DataStorage<UnitigData> , false>
				forward = active.getSuccessors();
		auto it = forward.begin();
		if (it != forward.end()) {
			for (it; it != forward.end(); ++it) {
				cout << "Child here" << endl;
				Kmer child = GraphTraverser::getKmer(*it);
				if (child.toString() == km.toString()) {
					//cycle including s, abort!
					cout << "Cycle,abort" << endl;
					return pair;
				} else {
					cout << "Child seen" << endl;
					DataAccessor < UnitigData > *da = it->getData(); // Get DataAccessor from unitig
					UnitigData* data = da->getData(*it); // Get boolean from DataAccessor
					data->set_seen();
					seen_but_unvisited++;

					bool all_visited = true;
					BackwardCDBG<DataAccessor<UnitigData> , DataStorage<
							UnitigData> , false> backward =
							it->getPredecessors();
					for (auto back = backward.begin(); back != backward.end(); ++back) {
						cout << "Check predecessor" << endl;
						DataAccessor<UnitigData>* pre = back->getData(); // Get DataAccessor from unitig
						UnitigData* predata = pre->getData(*back); // Get boolean from DataAccessor

						if (!predata->is_visited()) {
							all_visited = false;
							break;
						}
					}
					if (all_visited) {
						cout << "Push child to S" << endl;
						S.push_back(*it);
					}
				}
			}
		} else {
			//tip, abort!
			cout << "Tip,abort" << endl;
			return pair;

		}

		if (S.size() == 1 && seen_but_unvisited == 1) {
			//todo: need to check if the remaining node is not a direct successor of start!
			pair.push_back(active);
			pair.push_back(S[0]);
			cout << "Pair found!" << endl;
		}

	}
	return pair;
}



/*
 * testfunction: for all simple bubbles, if hamming distance of alternatives is 1 and joint color set equals the colors of the anchors, remove both unitigs and replace by sequence represented by more strains
 * for each bubble, test if this reduces the number of unitigs by exactly 1
 * TO IMPLEMENT!
 */
void GraphTraverser::removeSimpleBubbles(){
	//iterate over unitigs without specific order
		for (ColoredCDBG<UnitigData>::iterator it = cdbg.begin(); it != cdbg.end(); ++it) {

			//ToDo: if not visited yet!

			DataAccessor<UnitigData>* da = it->getData(); // Get DataAccessor from unitig
			UnitigData* data = da->getData(*it); // Get boolean from DataAccessor
			data->set_visited();
			ForwardCDBG<DataAccessor<UnitigData> , DataStorage<UnitigData> , false>
					forward = it->getSuccessors();
			size_t size = distance(forward.begin(), forward.end());

			//Note to me: If a unitig has only one predecessor and only one successor, then they should be the same! It is a self loop...
			vector<Kmer> kmerStore;
			if (size == 2) {
				Kmer kmi = GraphTraverser::getKmer(*(it));
				for (auto forwardit = forward.begin(); forwardit != forward.end(); ++forwardit) {
					DataAccessor<UnitigData>* forward_da = forwardit->getData();
					UnitigData* forward_data = da->getData(*forwardit);
					forward_data->set_visited();
					ForwardCDBG<DataAccessor<UnitigData> ,
							DataStorage<UnitigData> , false> nextforward =
							forwardit->getSuccessors();
					size_t size = distance(nextforward.begin(), nextforward.end());

					if (size != 1) {
						continue;
					} else {
						//we need the first kmer of the end unitig
						Kmer k = GraphTraverser::getKmer(*(nextforward.begin()));
						kmerStore.push_back(k);
					}
				}
				//check if we finished with the same unitig
				if (kmerStore.size() != 2) {
					//cout << kmerStore.size() << endl;
				} else {
					if (kmerStore[0] == kmerStore[1]) {
						if (kmerStore[0].toString() != kmi.toString()) {
							//remove both successors of the initial unitig, replace them
							//ToDo
						}
					}
				}
			}
		}
}





/*
 * Given the cDBG graph, count the number of simple bubbles -> first unitig is branching into two nodes, both nodes unite again in the same fourth untig
 *
 * Input: void, (the cdBG is a global object in this class)
 * Output: count of simple bubbles
 */
int GraphTraverser::countSimpleBubbles() {
	int bubbleCounter = 0;

	//iterate over unitigs without specific order
	for (ColoredCDBG<UnitigData>::iterator it = cdbg.begin(); it != cdbg.end(); ++it) {

		//ToDo: if not visited yet!

		DataAccessor<UnitigData>* da = it->getData(); // Get DataAccessor from unitig
		UnitigData* data = da->getData(*it); // Get boolean from DataAccessor
		data->set_visited();
		ForwardCDBG<DataAccessor<UnitigData> , DataStorage<UnitigData> , false>
				forward = it->getSuccessors();
		size_t size = distance(forward.begin(), forward.end());

		//Note to me: If a unitig has only one predecessor and only one successor, then they should be the same! It is a self loop...
		vector<Kmer> kmerStore;
		if (size == 2) {
			Kmer kmi = GraphTraverser::getKmer(*(it));
			for (auto forwardit = forward.begin(); forwardit != forward.end(); ++forwardit) {
				DataAccessor<UnitigData>* forward_da = forwardit->getData();
				UnitigData* forward_data = da->getData(*forwardit);
				forward_data->set_visited();
				ForwardCDBG<DataAccessor<UnitigData> ,
						DataStorage<UnitigData> , false> nextforward =
						forwardit->getSuccessors();
				size_t size = distance(nextforward.begin(), nextforward.end());

				if (size != 1) {
					continue;
				} else {
					//we need the first kmer of the end unitig
					Kmer k = GraphTraverser::getKmer(*(nextforward.begin()));
					kmerStore.push_back(k);
				}
			}
			//check if we finished with the same unitig
			if (kmerStore.size() != 2) {
				//cout << kmerStore.size() << endl;
			} else {
				if (kmerStore[0] == kmerStore[1]) {
					if (kmerStore[0].toString() != kmi.toString()) {
						bubbleCounter++;
					}
				}
			}
		}
	}
	return bubbleCounter;
}




ColoredCDBG<UnitigData>& GraphTraverser::getGraph() {
	return cdbg;
}



/*
 * Helperfunction
 * Evaluate given UnitigMap, comparing first and last kmer in unitig
 * return smaller kmer
 */
Kmer GraphTraverser::getKmer(auto& map) {
	Kmer head = map.getUnitigHead().rep();
	Kmer tail = map.getUnitigTail().rep();
	Kmer km = (head < tail ? head : tail);

	return km;
}




//check the number of unitigs in the graph
size_t GraphTraverser::sizeOfGraph() {
	size_t x = cdbg.size();
	cout << "Size of resulting graph determined by cdbg.size(): " << x << endl;
	return x;
}



//Dummy testfunction
void GraphTraverser::traverseGraph() {
for (auto& unitig : cdbg) { // Iterate over unitigs of a colored de Bruijn graph
	DataAccessor<UnitigData>* da = unitig.getData(); // Get DataAccessor from unitig
	//UnitigData* data = da->getData(unitig); // Get boolean from DataAccessor
	//UnitigColors* uc = da->getUnitigColors(unitig);

//	set<string> colors;
//	for(UnitigColors::const_iterator it = uc->begin(); it != uc->end(); ++it) {
//		size_t color = it->getColorID(unitig.len);
//		colors.insert(cdbg.getColorName(color));
//	}
//
//
//
//	if (colors.size() >= 2) {
//		cout << unitig.toString() << endl;
//		cout << colors.size() << endl;
//		cout << "------" << endl;
//	}
}
}


void GraphTraverser::BFS_Recursive(const UnitigMap<DataAccessor<UnitigData>, DataStorage<UnitigData>, false>& ucm, int depth, UnionFind& uf, unordered_map<string,int>& unitig_ids){
	if (depth > 0){

    	//reduce depth of traversal
    	depth--;


		Kmer left = getKmer(ucm);

		for (auto& successor : ucm.getSuccessors()){ // Iterate over successors of a unitig

			DataAccessor<UnitigData>* da = successor.getData(); // Get DataAccessor from unitig successor
			UnitigData* data = da->getData(successor); // Get boolean from DataAccessor

				Kmer right = getKmer(successor);
				uf.merge(unitig_ids[left.toString()],unitig_ids[right.toString()]);
		}

		for (auto& predecessor : ucm.getPredecessors()){ // Iterate over predecessors of a unitig

			DataAccessor<UnitigData>* da = predecessor.getData(); // Get DataAccessor from unitig predecessor
			UnitigData* data = da->getData(predecessor); // Get boolean from DataAccessor

				Kmer right = getKmer(predecessor);
				uf.merge(unitig_ids[left.toString()],unitig_ids[right.toString()]);
		}



    	// Traverse successors
    	for (auto& successor : ucm.getSuccessors()) BFS_Recursive(successor, depth, uf, unitig_ids);
    	// Traverse predecessors
    	for (auto& predecessor : ucm.getPredecessors()) BFS_Recursive(predecessor, depth, uf, unitig_ids);
    }
}


void GraphTraverser::BFS_Neighborhood(const UnitigMap<DataAccessor<UnitigData>, DataStorage<UnitigData>, false>& ucm, int depth, vector<UnitigMap<DataAccessor<UnitigData>, DataStorage<UnitigData>, false>>& neighbors){

	if (depth > 0){

    	//reduce depth of traversal
    	depth--;

		for (auto& successor : ucm.getSuccessors()){ // Iterate over successors of a unitig

			DataAccessor<UnitigData>* da = successor.getData(); // Get DataAccessor from unitig successor
			UnitigData* data = da->getData(successor); // Get boolean from DataAccessor

			if (data->is_not_visited()){ // If boolean indicates the unitig was not visited

				data->set_visited(); // Set boolean to indicate unitig was visited
				//Kmer right = getKmer(successor);
				//neighbors.insert(right);
				neighbors.push_back(successor);
			}
		}

		for (auto& predecessor : ucm.getPredecessors()){ // Iterate over predecessors of a unitig

			DataAccessor<UnitigData>* da = predecessor.getData(); // Get DataAccessor from unitig predecessor
			UnitigData* data = da->getData(predecessor); // Get boolean from DataAccessor

			if (data->is_not_visited()){ // If boolean indicates the unitig was not visited

				data->set_visited(); // Set boolean to indicate unitig was visited
				//Kmer right = getKmer(predecessor);
				//neighbors.insert(right);
				neighbors.push_back(predecessor);
			}
		}



    	// Traverse successors
    	for (auto& successor : ucm.getSuccessors()) BFS_Neighborhood(successor, depth, neighbors);
    	// Traverse predecessors
    	for (auto& predecessor : ucm.getPredecessors()) BFS_Neighborhood(predecessor, depth, neighbors);
    }
}

void GraphTraverser::cleanMarking(){

    for (auto& unitig : cdbg){

        DataAccessor<UnitigData>* da_ucm = unitig.getData();
        UnitigData* data_ucm = da_ucm->getData(unitig);

        data_ucm->set_not_seen_visited();
    }
}

