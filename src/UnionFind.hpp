/*
 * UnionFind.hpp
 *
 *  Created on: 27 Feb 2018
 *      Author: nina
 */

#ifndef UNIONFIND_HPP_
#define UNIONFIND_HPP_

#include <unordered_map>
#include <iostream>
#include <vector>

class UnionFind {

public:

	int *id, cnt, *sz;

	// Create an empty union find data structure with N isolated sets.
	UnionFind(int N);
	~UnionFind();

	// Return the id of component corresponding to object p.
	int find(int p);

	// Replace sets containing x and y with their union.
	void merge(int x, int y);

	// Are objects x and y in the same set?
	bool connected(int x, int y);

	// Return the number of disjoint sets.
	int count();

	// Return the number of connected components of a specific size (size=number of elements)
	std::unordered_map<int,int> countSize(int N);

	std::unordered_map<int,std::vector<int>> compileComponents(int N);

private:
	std::unordered_map<int,int> histogram(std::unordered_map<int,int> app);

};







#endif /* UNIONFIND_HPP_ */
