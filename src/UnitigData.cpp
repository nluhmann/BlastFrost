/*
 * MyBool.cpp
 *
 *  Created on: 19 Mar 2018
 */

#include "UnitigData.hpp"


UnitigData::UnitigData() : b(NOT_VISITED_SEEN),a(NOT_ANCHOR),c(NOT_CORE),t(NOT_TARGET) {} // Initiate the boolean to "not visited"


void UnitigData::clear(const UnitigColorMap<UnitigData>&  um_dest){

	set_not_seen_visited(); // Set the unitig to "not visited"
	set_not_anchor();
	set_not_core();
	set_not_on_target_path();
}




void UnitigData::concat(const UnitigColorMap<UnitigData>& um_dest, const UnitigColorMap<UnitigData>& um_src){


    //DataAccessor<UnitigData>* da = um_dest.getData(); // Get DataAccessor from unitig matching um_dest
    //UnitigData* data = da->getData(um_dest); // Get boolean from DataAccessor

    //data->set_not_seen_visited(); // Set the unitig to "not visited"
    //data->set_not_anchor();
    //data->set_not_core();
    //data->setClusterID(0);

    set_not_seen_visited(); // Set the unitig to "not visited"
    set_not_anchor();
    set_not_core();
    set_not_on_target_path();
}



void UnitigData::extract(const UnitigColors* uc_dest, const UnitigColorMap<UnitigData>& um_src, const bool last_extraction ){

    set_not_seen_visited();
    set_not_anchor();
    set_not_core();
    set_not_on_target_path();
}


void UnitigData::merge(const UnitigColors& uc_dest, const UnitigColorMap<UnitigData>& um_dest, const const_UnitigColorMap<UnitigData> & um_src){
	//DataAccessor<UnitigData>* da = um_dest.getData(); // Get DataAccessor from unitig matching um_dest
	//UnitigData* data = da->getData(um_dest); // Get boolean from DataAccessor

	//data->set_not_seen_visited(); // Set the unitig to "not visited"
	//data->set_not_anchor();
	//data->set_not_core();

	set_not_seen_visited();
	set_not_anchor();
	set_not_core();
	set_not_on_target_path();
}



void UnitigData::toString() const {
            cout << "Unitig visited = " << (is_visited() ? "true" : "false") << endl;
            cout << "Unitig seen = " << (is_seen() ? "true" : "false") << endl;
            cout << "Unitig is anchor = " << (is_anchor() ? "true" : "false") << endl;
            cout << "Unitig is core = " << (is_core() ? "true" : "false") << endl;
}




