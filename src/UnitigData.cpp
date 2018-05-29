/*
 * MyBool.cpp
 *
 *  Created on: 19 Mar 2018
 */

#include "UnitigData.hpp"


UnitigData::UnitigData() : b(NOT_VISITED_SEEN),a(NOT_ANCHOR),c(NOT_CORE),clusterID(0) {} // Initiate the boolean to "not visited"

// Join method for ColoredCDBG
void UnitigData::join(const UnitigColorMap<UnitigData>& um_dest, const UnitigColorMap<UnitigData>& um_src){

	// When joining the unitig matching um_src to the unitig matching um_dest,
    // we set um_dest to "not visited" because it will be a new unitig in the graph.

    DataAccessor<UnitigData>* da = um_dest.getData(); // Get DataAccessor from unitig matching um_dest
    UnitigData* data = da->getData(um_dest); // Get boolean from DataAccessor

    data->set_not_seen_visited(); // Set the unitig to "not visited"
    data->set_not_anchor();
    data->set_not_core();
    data->setClusterID(0);
}


// Sub method for ColoredCDBG
void UnitigData::sub(UnitigData* new_data, const UnitigColors& uc_dest, const UnitigMapBase& um_dest,
		const UnitigColorMap<UnitigData>& um_src, const bool last_extraction){

	// This function creates a new unitig which is a sub-unitig from um_src
	// The new unitig created is set to "not visited" as a measure of precaution
	// (it is already initialed by default to "not visited" in the constructor)

    new_data->set_not_seen_visited();
    new_data->set_not_anchor();
    new_data->set_not_core();
    new_data->setClusterID(0);
}


void UnitigData::toString() const {
            cout << "Unitig visited = " << (is_visited() ? "true" : "false") << endl;
            cout << "Unitig seen = " << (is_seen() ? "true" : "false") << endl;
            cout << "Unitig is anchor = " << (is_anchor() ? "true" : "false") << endl;
            cout << "Unitig is core = " << (is_core() ? "true" : "false") << endl;
            cout << "ClusterID of unitig if anchor: " << clusterID << endl;
}


void UnitigData::setClusterID(int x) {
	clusterID = x;
}

int UnitigData::getClusterID() {
	return clusterID;
}




