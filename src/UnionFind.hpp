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

	// Create an empty union find data structure with N isolated sets.
	UnionFind(const int N);
	~UnionFind();

	// Return the id of component corresponding to object p.
	int find(int p);

	// Replace sets containing x and y with their union.
	void merge(const int x, const int y);

	// Are objects x and y in the same set?
	inline bool connected(const int x, const int y) {

		return UnionFind::find(x) == UnionFind::find(y);
	}

	// Return the number of connected components
	inline int count() const {

		return cnt;
	}

	// Return the number of connected components of a specific size (size=number of elements)
	std::unordered_map<int,int> countSize(const int N) const;

	std::unordered_map<int,std::vector<int>> compileComponents(const int N) const;

private:

	int *id, cnt, *sz;

	std::unordered_map<int,int> histogram(const std::unordered_map<int,int>& app) const;

};







#endif /* UNIONFIND_HPP_ */
