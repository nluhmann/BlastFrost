/*
 * BubbleExplorer.hpp
 *
 *  Created on: 14 Jan 2019
 *      Author: nina
 */

#ifndef BUBBLEEXPLORER_HPP_
#define BUBBLEEXPLORER_HPP_




#include <ColoredCDBG.hpp>
#include "UnitigData.hpp"

#include <stack>

class BubbleExplorer{


private:

	//reference to the graph that is going to be traversed
	ColoredCDBG<UnitigData>& cdbg;


public:

	BubbleExplorer(ColoredCDBG<UnitigData>& cdbg);


	void exploreBubble(const string& left, const string& right, const int threshold);

	vector<vector<UnitigColorMap<UnitigData>>> DFS_Iterative(const UnitigColorMap<UnitigData>& start, const Kmer& stop, const unsigned int threshold);

	void printPaths(const vector<vector<UnitigColorMap<UnitigData>>> allPaths);























};



#endif /* BUBBLEEXPLORER_HPP_ */
