/*
 * UnionFind.cpp
 *
 *  Created on: 27 Feb 2018
 *      Author: nina
 */

#include "UnionFind.hpp"

using namespace std;

// ToDo: need to change constructor, initialize each UnitigMap pointer as a set.
UnionFind::UnionFind(int N) {
    cnt = N;
    id = new int[N];
    sz = new int[N];
    for (int i = 0; i<N; i++){
    	id[i] = i;
    	sz[i] = 1;
    }
}

UnionFind::~UnionFind() {
	delete[] id;
	delete[] sz;
}

// Return the id of component corresponding to object p.
int UnionFind::find(int p) {
    int root = p;
    while (root != id[root]){
    	root = id[root];
    }
    while (p != root) {
    	int newp = id[p];
    	id[p] = root;
    	p = newp;
    }
    return root;
}

// Replace sets containing x and y with their union.
void UnionFind::merge(int x, int y) {
    int i = UnionFind::find(x);
    int j = UnionFind::find(y);
    if (i == j) return;
    // make smaller root point to larger one
    if (sz[i] < sz[j]) {
    	id[i] = j;
    	sz[j] += sz[i];
    } else {
    	id[j] = i;
    	sz[i] += sz[j];
    }
    cnt--;
}

// Are objects x and y in the same set?
bool UnionFind::connected(int x, int y) {
	return UnionFind::find(x) == UnionFind::find(y);
}

// Return the number of connected components
int UnionFind::count() {
	return cnt;
}

// Return the number of connected components of a specific size
unordered_map<int,int> UnionFind::countSize(int N) {
	unordered_map<int,int> size;

	for(int i = 0; i<N; i++){
		//find the root for element i
		int root = i;
		while (root != id[root]){
			root = id[root];
		}
		//we found another element of the component represented by root
		int& x = size[root];
		x++;
	}

	unordered_map<int,int> y = UnionFind::histogram(size);
	return y;
}


unordered_map<int,int> UnionFind::histogram(unordered_map<int,int> app) {
	unordered_map<int,int> hist;
	for(auto& it : app){
		hist[it.second]++;
	}
	return hist;
}


// Return connected components as lists of IDs belonging to the same CC
unordered_map<int,vector<int>> UnionFind::compileComponents(int N){
	unordered_map<int,vector<int>> ccs;
	for(int i = 0; i<N; i++){
		//find the root for element i
		int root = i;
		while (root != id[root]){
			root = id[root];
		}
		//we found another element of the component represented by root
		ccs[root].push_back(i);
	}

	return ccs;
}





