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


#include <string>
#include <unordered_map>
#include <algorithm>
#include <vector>
#include <iterator>
#include <set>

class GraphTraverser{
private:

	//reference to the graph that is going to be traversed
	ColoredCDBG<UnitigData>& cdbg;


public:

	GraphTraverser(ColoredCDBG<UnitigData>& cdbg);

	void search(string query, int k);

};








#endif /* GRAPHTRAVERSER_HPP_ */
