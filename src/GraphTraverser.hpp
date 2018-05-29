/*
 * GraphTraverser.hpp
 *
 *  Created on: 26 Feb 2018
 *      Author: nina
 */

#ifndef GRAPHTRAVERSER_HPP_
#define GRAPHTRAVERSER_HPP_



#include <bifrost/ColoredCDBG.hpp>
#include "UnitigData.hpp"
#include "UnionFind.hpp"


#include <string>
#include <unordered_map>
#include <algorithm>
#include <vector>
#include <iterator>
#include <set>

class GraphTraverser{
private:

	vector<UnitigColorMap<UnitigData>> findSuperBubble(const UnitigColorMap<UnitigData>& unitig);

	//reference to the graph that is going to be traversed
	ColoredCDBG<UnitigData>& cdbg;


public:

	GraphTraverser(ColoredCDBG<UnitigData>& cdbg);

	Kmer getKmer(auto& map);

	void traverseGraph();

	size_t sizeOfGraph();

	int countSimpleBubbles();

	int countSuperBubbles();

	int countConnectedComponents();

	void removeSimpleBubbles();

	ColoredCDBG<UnitigData>& getGraph();

	void BFS_Recursive(const UnitigMap<DataAccessor<UnitigData>, DataStorage<UnitigData>, false>& ucm, int depth, UnionFind& uf, unordered_map<string,int>& unitig_ids);

	void BFS_Neighborhood(const UnitigMap<DataAccessor<UnitigData>, DataStorage<UnitigData>, false>& ucm, int depth, vector<UnitigMap<DataAccessor<UnitigData>, DataStorage<UnitigData>, false>>& neighbors);

	void cleanMarking();


};








#endif /* GRAPHTRAVERSER_HPP_ */
