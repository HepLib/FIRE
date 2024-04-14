#pragma once

#include "common.h"
typedef uint32_t virt_t;
#include <iomanip>
class point_fast;

//#include "parallel_hashmap/phmap.h"

typedef vector<pair<vector<COEFF>, point_fast > > ibp_type;

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
    #ifndef SMALL_POINT
    char ww[24];
    #else
    char ww[16];
    #endif

    unsigned short h1() const {
        return *reinterpret_cast<const unsigned short*>(ww + sizeof(point) - 2);
    };

    unsigned short* h1p() {
        return reinterpret_cast<unsigned short*>(ww + sizeof(point) - 2);
    };

    string number() const;

    unsigned short s_number() const { return (h1() >> 1); }

    bool virt() const { return (ww[sizeof(point) - 3] == 0); }

    int level() const;

    void safe_string(char *buf) const;

    vector<t_index> get_vector() const;

    bool is_zero() const;

    friend bool operator==(const point &p1, const point &p2) {
        #ifndef SMALL_POINT
            if (reinterpret_cast<const uint64_t *>(p1.ww)[2] ^ reinterpret_cast<const uint64_t *>(p2.ww)[2]) return false;
            if (reinterpret_cast<const uint64_t *>(p1.ww)[1] ^ reinterpret_cast<const uint64_t *>(p2.ww)[1]) return false;
            if (reinterpret_cast<const uint64_t *>(p1.ww)[0] ^ reinterpret_cast<const uint64_t *>(p2.ww)[0]) return false;
            return true;
        #else
            return !(reinterpret_cast<const uint128_t *>(p1.ww)[0] ^ reinterpret_cast<const uint128_t *>(p2.ww)[0]);
        #endif
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
        if (reinterpret_cast<const uint64_t *>(p1.ww)[2] ^ reinterpret_cast<const uint64_t *>(p2.ww)[2]) {
            return reinterpret_cast<const uint64_t *>(p1.ww)[2] < reinterpret_cast<const uint64_t *>(p2.ww)[2];
        }
        if (reinterpret_cast<const uint64_t *>(p1.ww)[1] ^ reinterpret_cast<const uint64_t *>(p2.ww)[1]) {
            return reinterpret_cast<const uint64_t *>(p1.ww)[1] < reinterpret_cast<const uint64_t *>(p2.ww)[1];
        }
        return reinterpret_cast<const uint64_t *>(p1.ww)[0] < reinterpret_cast<const uint64_t *>(p2.ww)[0];
        #else
            return reinterpret_cast<const uint128_t *>(p1.ww)[0] < reinterpret_cast<const uint128_t *>(p2.ww)[0];
        #endif
    }


    friend ostream &operator<<(ostream &out, const point &p);

    static map<unsigned short,
            vector< // because there can be multiple
                    pair<vector<t_index>,
                            vector<
                                    pair<
                                            vector<
                                                    pair<COEFF, point_fast>
                                            >,
                                            t_index>
                            >
                    >
            >
    > ibases;

    static map<unsigned short,
            pair<vector<t_index>,
                    vector<
                            pair<
                                    vector<
                                            pair<COEFF, point_fast>
                                    >,
                                    t_index>
                    >
            >
    > dbases;

    static vector<ibp_type> ibps;
    
//    friend size_t hash_value(const point &p) {
//        #ifndef SMALL_POINT
//        return phmap::HashState().combine(0, reinterpret_cast<const uint64_t *>(p.ww)[2], reinterpret_cast<const uint64_t *>(p.ww)[1], reinterpret_cast<const uint64_t *>(p.ww)[0]);
//        #else
//        return std::has(reinterpret_cast<const uint128_t *>(p.ww)[0]);
//        #endif
//    }
}
#ifndef SMALL_POINT
__attribute__((aligned (8))); // size is 24, so align by 8
#else
__attribute__((aligned (16))); // size is 16, let's alighn by 16 too
#endif

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

//typedef phmap::parallel_node_hash_map<pair<unsigned int, unsigned int>, bool> umap_type;
//typedef phmap::parallel_node_hash_map<point, pair<unsigned char, pc_pair_ptr_vec>> pmap_type;
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
