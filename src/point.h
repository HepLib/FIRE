#pragma once

#include "common.h"
typedef uint32_t virt_t;
#include <iomanip>
class point_fast;

typedef vector<pair<vector<COEFF>, point_fast > > ibp_type;


// sector variants permutation product sum coefficient indices powers
// vector - because there can be multiple
typedef map<sector_count_t, vector<pair<vector<t_index>, vector<pair<vector<pair<COEFF, point_fast>>, t_index>>>>> ibases_t;


// sector permutation product sum coefficient indices powers
typedef map<sector_count_t, pair<vector<t_index>, vector<pair<vector<pair<COEFF, point_fast>>, t_index>>>> dbases_t;

#ifndef SMALL_POINT
static const size_t POINT_ALIGN = 8;  //Default size is 24, so align by 8
static const size_t POINT_SIZE = (sizeof(sector_count_t) + MAX_IND + 7) & -POINT_ALIGN;  //Bit magic for rounding to the nearest multiple of POINT_ALIGN
#else
static const size_t POINT_SIZE = 16;
static const size_t POINT_ALIGN = 16;//Small-point size is 16, let's alighn by 16 too
#endif



static const size_t POINT_SECTOR_OFFSET_FROM_BACK = sizeof(sector_count_t);

//As per the comments in point, the byte at MAX_IND-1 is 0 for virtual,
//and the previous 4 bytes (virt_t is a uint32) are where the virtual number is stored
static const size_t VIRT_ZERO_OFFSET = MAX_IND-1;
static const size_t VIRT_OFFSET = VIRT_ZERO_OFFSET-sizeof(virt_t);

/**
 * @brief Ordered point class.
 *
 * Contains point's degrees already multiplied by ordering matrix and calculated sector number.
 * There also are virtual points - there are not real integrals, but some masking expressions of tails in lower sectors.
 * Point is designed to fit into the 24 = 3*8 byte size
 */
class point {
public:
    static vector<set<vector<t_index> > > preferred;
    static vector<set<point_fast> > preferred_fast;
    static bool print_g;
    point();
    point(const point &p);
    point(const point_fast &pf, SECTOR ssector = -1);
    point(const point &p, const point_fast &v, SECTOR ssector, bool not_preferred);
    point &operator=(const point &);
    point(const vector<t_index> &v, virt_t virt = 0, SECTOR ssector = -1);
    point(char *buf);
    char ww[POINT_SIZE];

    sector_count_t h1() const {
      return *reinterpret_cast<const sector_count_t*>(ww + std::size(ww) - POINT_SECTOR_OFFSET_FROM_BACK);
    };

    sector_count_t* h1p() {
        return reinterpret_cast<sector_count_t*>(ww + std::size(ww)  - POINT_SECTOR_OFFSET_FROM_BACK);
    };

    string number() const;

    sector_count_t s_number() const { return (h1() >> 1); }

    bool virt() const { return (ww[VIRT_ZERO_OFFSET] == 0); }

    int level() const;

    void safe_string(char *buf) const;

    vector<t_index> get_vector() const;

    bool is_zero() const;

    friend bool operator==(const point &p1, const point &p2) {
      //All we care about is exact equality between the buffers, so just check for that
      return memcmp(p1.ww,p2.ww,std::size(p1.ww))==0;
    }

    friend bool operator!=(const point &p1, const point &p2) {
        #ifndef SMALL_POINT
            return (!(p1 == p2));
        #else
            return reinterpret_cast<const uint128_t *>(p1.ww)[0] ^ reinterpret_cast<const uint128_t *>(p2.ww)[0];
        #endif
    }

    friend bool operator<=(const point &p1, const point &p2) {
        #ifndef SMALL_POINT
            return (p1 == p2 || p1 < p2);
        #else
            return reinterpret_cast<const uint128_t *>(p1.ww)[0] <= reinterpret_cast<const uint128_t *>(p2.ww)[0];
        #endif
    }

    friend bool operator<(const point &p1, const point &p2) {     
#ifndef SMALL_POINT
      //This should maybe become a static const, but I believe almost all other use cases
      //have been switched over to using memcmp, memset, or similar byte-level manips
      size_t ww_as_int64_size = std::size(p1.ww)/sizeof(uint64_t);
      //Other sections of the process rely on this specific order of comparisons, so
      //can't just do memcmp
      for(size_t i = 1;i<ww_as_int64_size;++i){
	if (reinterpret_cast<const uint64_t *>(p1.ww)[ww_as_int64_size-i] ^ reinterpret_cast<const uint64_t *>(p2.ww)[ww_as_int64_size-i]) {
            return reinterpret_cast<const uint64_t *>(p1.ww)[ww_as_int64_size-i] < reinterpret_cast<const uint64_t *>(p2.ww)[ww_as_int64_size-i];
        }
	
      }
      return reinterpret_cast<const uint64_t *>(p1.ww)[0] < reinterpret_cast<const uint64_t *>(p2.ww)[0];
#else
            return reinterpret_cast<const uint128_t *>(p1.ww)[0] < reinterpret_cast<const uint128_t *>(p2.ww)[0];
        #endif
    }


    friend ostream &operator<<(ostream &out, const point &p);

    static ibases_t ibases;

    static dbases_t dbases;

    static vector<ibp_type> ibps;
}
__attribute__((aligned (POINT_ALIGN))); 

class point_fast {
public:
    ~point_fast();

    point_fast();

    point_fast(const point_fast &v);

    point_fast(const point &p);

    point_fast(const vector<t_index> &v);

    point_fast &operator=(const point_fast &);

    friend bool operator==(const point_fast &p1, const point_fast &p2) {
        return !memcmp(p1.buf, p2.buf, MAX_IND);
    }

    friend bool operator!=(const point_fast &p1, const point_fast &p2) {
        return memcmp(p1.buf, p2.buf, MAX_IND);
    }

    friend bool operator<=(const point_fast &p1, const point_fast &p2) {
        return (memcmp(p1.buf, p2.buf, MAX_IND) <= 0);
    }

    friend bool operator<(const point_fast &p1, const point_fast &p2) {
        return (memcmp(p1.buf, p2.buf, MAX_IND) < 0);
    }

    friend point_fast operator+(const point_fast &p1, const point_fast &p2) {
        point_fast result;
        const t_index *pos1 = p1.buf;
        const t_index *pos2 = p2.buf;
        for (unsigned short i = 0; i != MAX_IND; ++i, ++pos1, ++pos2) result.buf[i] = *pos1 + *pos2;
        return result;
    }

    point_fast degree() const;

    SECTOR sector_fast() const;

    /** Buffer for point's indices. */
    t_index buf[MAX_IND];
};


bool over_fast(const point_fast &l, const point_fast &r);

set<point_fast> level_points_fast(point_fast s, const unsigned int pos, const unsigned int neg);

typedef pair<point,COEFF> pc_pair;
typedef list<pc_pair> pc_pair_list;
typedef vector<pc_pair> pc_pair_vec;
typedef std::shared_ptr<const pc_pair> pc_pair_ptr;
typedef vector<pc_pair_ptr> pc_pair_ptr_vec;
typedef list<pc_pair_ptr> pc_pair_ptr_lst;

template<typename P, typename C>
inline pc_pair_ptr make_pc_ptr(P && p, C && c) {
    return std::make_shared<const pc_pair>(std::forward<P>(p), std::forward<C>(c));
    //return pc_pair_ptr(new pc_pair(std::forward<P>(p), std::forward<C>(c)));
}

template<typename PC>
inline pc_pair_ptr make_pc_ptr(PC && pc) {
    return std::make_shared<const pc_pair>(std::forward<PC>(pc));
    //return pc_pair_ptr(new pc_pair(std::forward<PC>(pc)));
}

//#include "parallel_hashmap/phmap.h"
//namespace std {
//    template<> struct hash<point>
//    {
//        std::size_t operator()(point const &p) const
//        {
//            phmap::HashState hc;
//            std::size_t seed = 0;
//            for(int i=0; i<POINT_SIZE; i++) {
//                seed = hc.combine(seed, p.ww[i]);
//            }
//            return seed;
//        }
//    };
//    template<> struct hash<pair<unsigned int, unsigned int>>
//    {
//        std::size_t operator()(pair<unsigned int, unsigned int> const &p) const
//        {
//            phmap::HashState hc;
//            std::size_t seed = 0;
//            hc.combine(seed, p.first);
//            hc.combine(seed, p.second);
//            return seed;
//        }
//    };
//}
//typedef phmap::parallel_node_hash_map<pair<unsigned int, unsigned int>, bool> umap_type;
//typedef phmap::parallel_node_hash_map<point, pair<unsigned char, pc_pair_ptr_vec>> pmap_type;

//#include "sparsepp/spp.h"
//namespace std
//{
//    template<> struct hash<point>
//    {
//        std::size_t operator()(point const &p) const
//        {
//            std::size_t seed = 0;
//            for(int i=0; i<POINT_SIZE; i++) {
//                spp::hash_combine(seed, p.ww[i]);
//            }
//            return seed;
//        }
//    };
//    template<> struct hash<pair<unsigned int, unsigned int>>
//    {
//        std::size_t operator()(pair<unsigned int, unsigned int> const &p) const
//        {
//            std::size_t seed = 0;
//            spp::hash_combine(seed, p.first);
//            spp::hash_combine(seed, p.second);
//            return seed;
//        }
//    };
//}
//typedef spp::sparse_hash_map<pair<unsigned int, unsigned int>, bool> umap_type;
//typedef spp::sparse_hash_map<point, pair<unsigned char, pc_pair_ptr_vec>> pmap_type;

typedef map<pair<unsigned int, unsigned int>, bool> umap_type;
typedef map<point, pair<unsigned char, pc_pair_ptr_vec>> pmap_type;

struct DB {
    bool need_write = false;
    bool keep = false;
    int open_mode = 0;
    int opened = 0; // open times
    int status = 0; // 1: forward done, -1: substitue done
    int64_t lower_size;
    vector<point> lower;
    umap_type umap;
    pmap_type pmap;
    void clear();
};
extern map<int,DB> DBM; // DB Manager

bool database_exists(int number);
void open_database(int number, int mode=0, const set<point> *pset=NULL);
void close_database(int number);
int read_status(int sn);
