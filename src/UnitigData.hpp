/*
 * UnitigData.hpp
 *
 *  Created on: 19 Mar 2018
 */

#ifndef UNITIGDATA_HPP_
#define UNITIGDATA_HPP_


#include <ColoredCDBG.hpp>



class UnitigData : public CCDBG_Data_t<UnitigData>, CDBG_Data_t<UnitigData> {

    public:

		UnitigData();

		void clear(const UnitigColorMap<UnitigData>&  um_dest);

        // Join method for ColoredCDBG
		void concat(const UnitigColorMap<UnitigData>& um_dest, const UnitigColorMap<UnitigData>& um_src);

        // Sub method for ColoredCDBG
		void extract(const UnitigColors* uc_dest, const UnitigColorMap<UnitigData>& um_src, const bool last_extraction );

		void merge(const UnitigColors& uc_dest, const UnitigColorMap<UnitigData>& um_dest, const const_UnitigColorMap<UnitigData> & um_src);

        void toString() const;

        inline void set_visited() { b = VISITED; } // Set the boolean to "visited"
        inline void set_seen() { b = SEEN; } // Set the boolean to "seen"
        inline void set_not_seen_visited() { b = NOT_VISITED_SEEN; } // Set the boolean to "not seen and not visited"


        inline void set_on_target_path() { t = TARGET; } //set boolean to "visited and on target path"
        inline void set_not_on_target_path() { t = NOT_TARGET; }

        inline void set_anchor() { a = ANCHOR; }
        inline void set_not_anchor() {a = NOT_ANCHOR; }

        inline void set_core() { c = CORE; }
        inline void set_not_core() {c = NOT_CORE; }

        inline bool is_visited() const { return (b == VISITED); } // return if the boolean is "visited"
        inline bool is_not_visited() const { return (b != VISITED); } // return if the boolean is "not visited"

        inline bool is_seen() const { return (b == SEEN); } // return if the boolean is "seen"
        inline bool is_not_seen() const { return (b != SEEN); } // return if the boolean is "not seen"

        inline bool is_on_target_path() { return (t == TARGET); }
        inline bool is_not_on_target_path() { return (t != TARGET); }

        inline bool is_anchor() const {return (a == ANCHOR); }

        inline bool is_core() const {return (c == CORE); }

        void setClusterID(int x);
        int getClusterID();

    private:

        const static uint8_t NOT_VISITED_SEEN = 0x0;
        const static uint8_t VISITED = 0x1;
        const static uint8_t SEEN = 0x2;


        const static uint8_t TARGET = 0x0;
        const static uint8_t NOT_TARGET = 0x1;

        const static uint8_t ANCHOR = 0x0;
        const static uint8_t NOT_ANCHOR = 0x1;

        const static uint8_t CORE = 0x0;
        const static uint8_t NOT_CORE = 0x1;

        uint8_t b;

        uint8_t a;

        uint8_t c;

        uint8_t t;

        int clusterID;
};









#endif /* UNITIGDATA_HPP_ */
