/** @file point.cpp
 *  @author Alexander Smirnov
 *
 *  This file is a part of the FIRE package
 *  It contains the basic point class corresponding to a Feynman integral
 */

#include "point.h"

vector<set<vector<t_index> > > point::preferred;
vector<set<point_fast> > point::preferred_fast;
// numbered with sector numbers

bool point::print_g = false;

// we had to move them here to avoid circular dependencies
map<unsigned short, vector<pair<vector<t_index>, vector<pair<vector<pair<COEFF, point_fast> >, t_index> > > > > point::ibases;
//        sector       variants    permutation    product     sum         coefficient   indices       powers

map<unsigned short, pair<vector<t_index>, vector<pair<vector<pair<COEFF, point_fast> >, t_index> > > > point::dbases;
//        sector           permutation    product     sum         coefficient   indices       powers

vector<vector<pair<vector<COEFF>, point_fast > > > point::ibps;

point_fast::point_fast() {
    memset(buf, 0, MAX_IND);
}

point_fast::~point_fast() = default;

point_fast &point_fast::operator=(const point_fast &p) {
    memcpy(buf, p.buf, MAX_IND);
    return *this;
}

point_fast::point_fast(const point &p) {
    memset(buf, 0, MAX_IND);
    vector<t_index> v = p.get_vector();
    t_index *pos = buf;
    for (auto itr = v.begin(); itr != v.end(); ++itr, ++pos) *pos = *itr;
}

point_fast::point_fast(const vector<t_index> &v) {
    memset(buf, 0, MAX_IND);
    t_index *pos = buf;
    for (auto itr = v.begin(); itr != v.end(); ++itr, ++pos) *pos = *itr;
}

point_fast::point_fast(const point_fast &p) {
    memcpy(buf, p.buf, MAX_IND);
}

SECTOR point_fast::sector_fast() const {
    SECTOR result = 0;
    const t_index *pos_begin = this->buf;
    const t_index *pos_end = this->buf + common::dimension;
    while (pos_begin != pos_end) {
        if ((*pos_begin) > 0) {
            result = ((result << 1) ^ 1);
        } else {
            result = result << 1;
        }
        ++pos_begin;
    }
    return result;
}

point_fast point_fast::degree() const {
    point_fast result;
    const t_index *pos_old = this->buf;
    t_index *pos_new = result.buf;
    for (size_t i = 0; i != common::dimension; ++i, ++pos_new, ++pos_old)
        if ((*pos_old) > 0) {
            *pos_new = (*pos_old) - 1;
        } else {
            *pos_new = -(*pos_old);
        }
    return result;
}

bool over_fast(const point_fast &l, const point_fast &r) {
    for (unsigned int i = 0; i != common::dimension; ++i) {
        if (l.buf[i] < r.buf[i]) return false;
    }
    return true;
}

point::point(const point &p) {
    #ifndef SMALL_POINT
        reinterpret_cast<uint64_t *>(this->ww)[2] = reinterpret_cast<const uint64_t *>(p.ww)[2];
        reinterpret_cast<uint64_t *>(this->ww)[1] = reinterpret_cast<const uint64_t *>(p.ww)[1];
        reinterpret_cast<uint64_t *>(this->ww)[0] = reinterpret_cast<const uint64_t *>(p.ww)[0];
    #else
        reinterpret_cast<uint64_t *>(this->ww)[1] = reinterpret_cast<const uint64_t *>(p.ww)[1];
        reinterpret_cast<uint64_t *>(this->ww)[0] = reinterpret_cast<const uint64_t *>(p.ww)[0];
        //points in kyotocabinet are not aligned, so when taking data out we cannot rely on alignment
        //reinterpret_cast<uint128_t *>(this->ww)[0] = reinterpret_cast<const uint128_t *>(p.ww)[0];
    #endif
}

point &point::operator=(const point &p) {
    #ifndef SMALL_POINT
        reinterpret_cast<uint64_t *>(this->ww)[2] = reinterpret_cast<const uint64_t *>(p.ww)[2];
        reinterpret_cast<uint64_t *>(this->ww)[1] = reinterpret_cast<const uint64_t *>(p.ww)[1];
        reinterpret_cast<uint64_t *>(this->ww)[0] = reinterpret_cast<const uint64_t *>(p.ww)[0];
    #else
        reinterpret_cast<uint64_t *>(this->ww)[1] = reinterpret_cast<const uint64_t *>(p.ww)[1];
        reinterpret_cast<uint64_t *>(this->ww)[0] = reinterpret_cast<const uint64_t *>(p.ww)[0];
        //points in kyotocabinet are not aligned, so when taking data out we cannot rely on alignment
        //reinterpret_cast<uint128_t *>(this->ww)[0] = reinterpret_cast<const uint128_t *>(p.ww)[0];
    #endif
    return *this;
}

point::point() {
    #ifndef SMALL_POINT
        reinterpret_cast<uint64_t *>(this->ww)[2] = 0;
        reinterpret_cast<uint64_t *>(this->ww)[1] = 0;
        reinterpret_cast<uint64_t *>(this->ww)[0] = 0;
    #else
        reinterpret_cast<uint128_t *>(this->ww)[0] = 0;
    #endif
}


int point::level() const {
    return positive_index(common::ssectors[s_number()]);
}

#ifdef SMALL_POINT
    void shift_bits(unsigned short &bit_start, unsigned short &bit_end, unsigned short &current_byte, bool& split_byte) {
        bit_start += BITS_PER_INDEX;
        bit_end += BITS_PER_INDEX;
        if (bit_start >= 8) {
            bit_start -=8;
            ++current_byte;
        }
        if (bit_end >=8) {
            bit_end -= 8;
        }
        split_byte = (bit_start > bit_end);
    }

    inline void set_ww_and_shift(char* ww,char source, unsigned short &bit_start, unsigned short &bit_end, unsigned short &current_byte, bool& split_byte) {
        ww[sizeof(point) - 3 - current_byte] |= ((static_cast<unsigned char>(source) << (8 - BITS_PER_INDEX)) >> bit_start);
        if (split_byte) {
            ww[sizeof(point) - 3 - current_byte - 1] += (source << (7 - bit_end));
        }
        shift_bits(bit_start, bit_end, current_byte, split_byte);
    }

    inline char get_ww_and_shift(const char* ww, unsigned short &bit_start, unsigned short &bit_end, unsigned short &current_byte, bool& split_byte) {
        char result = static_cast<unsigned char>(static_cast<unsigned char>(ww[sizeof(point) - 3 - current_byte]) << bit_start) >> (8 - BITS_PER_INDEX);
        if (split_byte) {
            result |= (static_cast<unsigned char>(ww[sizeof(point) - 3 - current_byte - 1])) >> (7 - bit_end);
        }
        shift_bits(bit_start, bit_end, current_byte, split_byte);
        return result;
    }
#endif

point::point(const vector<t_index> &v, virt_t virt, SECTOR ssector) {
    #ifndef SMALL_POINT
        reinterpret_cast<uint64_t *>(ww)[2] = 0;
        reinterpret_cast<uint64_t *>(ww)[1] = 0;
        reinterpret_cast<uint64_t *>(ww)[0] = 0;
    #else
        reinterpret_cast<uint128_t *>(ww)[0] = 0;
    #endif


    SECTOR s = ssector;
    if (s == static_cast<SECTOR>(-1)) s = sector_fast(v);
    unsigned short sn;
    if (s == static_cast<SECTOR>(-2)) {
        sn = 1;
    } else {
        sn = common::sector_numbers_fast[s];
    }

    *h1p() = sn << 1;
    if (sn == 0) {
        cout << s << endl;
        cout << endl << "Sector 0" << endl << "Possible error in symmetries" << endl;
        abort();
    }

    if (virt != 0) {
        *reinterpret_cast<virt_t *>(ww + sizeof(point) - 7) = virt;
        ww[sizeof(point) - 3] = 0;
        return;
    }

    if ((sn == 1) || (sn == common::virtual_sector)) {
        #ifdef SMALL_POINT
            unsigned short bit_start = {0};
            unsigned short bit_end = {BITS_PER_INDEX - 1};
            unsigned short current_byte = {0};
            bool split_byte = false;
        #endif
        for (unsigned int i = 0; i != common::dimension; ++i) {
            #ifndef SMALL_POINT
                ww[sizeof(point) - i - 3] = 128 + v[i];
            #else
                set_ww_and_shift(ww, (1u<<(BITS_PER_INDEX-1)) + v[i], bit_start, bit_end, current_byte, split_byte);
            #endif
        }
    } else {
        auto ordering_now = common::orderings_fast[s].get();
        t_index degrees[MAX_IND];
        t_index *pos = static_cast<t_index *>(degrees);
        for (auto i = v.begin(); i != v.end(); ++pos, ++i) {
            if ((*i) > 0) {
                *pos = (*i) - 1;
            } else {
                *pos = -(*i);
            }
        }
        #ifdef SMALL_POINT
            unsigned short bit_start = {0};
            unsigned short bit_end = {BITS_PER_INDEX - 1};
            unsigned short current_byte = {0};
            bool split_byte = false;
        #endif
        for (unsigned int i = 0; i != common::dimension; ++i) {
            char pr = 1;  // it should be at least 1 to be higher than virtual
            for (unsigned int j = 0; j != common::dimension; ++j) {
                pr += ordering_now[i * common::dimension + j] * (degrees[j]);
            }
            #ifndef SMALL_POINT
                ww[sizeof(point) - 3 - i] = pr;
            #else
                set_ww_and_shift(ww, pr, bit_start, bit_end, current_byte, split_byte);
            #endif
        }
        if (preferred[sn].find(v) == preferred[sn].end()) {
            *h1p() |= 1;
        }
    }
}

point::point(const point_fast &pf, SECTOR ssector) {
    #ifndef SMALL_POINT
        reinterpret_cast<uint64_t *>(ww)[2] = 0;
        reinterpret_cast<uint64_t *>(ww)[1] = 0;
        reinterpret_cast<uint64_t *>(ww)[0] = 0;
    #else
        reinterpret_cast<uint128_t *>(ww)[0] = 0;
    #endif

    SECTOR s = ssector;
    if (s == static_cast<SECTOR>(-1)) { s = pf.sector_fast(); }
    unsigned short sn = common::sector_numbers_fast[s];
    *h1p() = sn << 1;
    if (sn == 0) {
        cout << s << endl;
        cout << endl << "Sector 0" << endl << "Possible error in symmetries" << endl;
        abort();
    }

    if (sn == common::virtual_sector) { // we never have sector 1 here
        #ifdef SMALL_POINT
            unsigned short bit_start = {0};
            unsigned short bit_end = {BITS_PER_INDEX - 1};
            unsigned short current_byte = {0};
            bool split_byte = false;
        #endif
        for (unsigned int i = 0; i != common::dimension; ++i) {
            #ifndef SMALL_POINT
                ww[MAX_IND - i - 1] = 128 + pf.buf[i];
            #else
                set_ww_and_shift(ww, (1u<<(BITS_PER_INDEX-1)) + pf.buf[i], bit_start, bit_end, current_byte, split_byte);
            #endif
        }
    } else {
        auto ordering_now = common::orderings_fast[s].get();
        point_fast degrees = pf.degree();

        // size should be equal to common::dimension in this case
        t_index *pos = ordering_now;

        #ifdef SMALL_POINT
            unsigned short bit_start = {0};
            unsigned short bit_end = {BITS_PER_INDEX - 1};
            unsigned short current_byte = {0};
            bool split_byte = false;
        #endif
        for (unsigned int i = 0; i != common::dimension; ++i) {
            #ifndef SMALL_POINT
                ww[MAX_IND - i - 1] = 1; // to be higher than virtual anyway even for preferred
                for (unsigned int j = 0; j != common::dimension; ++j) {
                    if (*pos++) ww[MAX_IND - i - 1] += (degrees.buf[j]);
                }
            #else
                char res = 1;
                for (unsigned int j = 0; j != common::dimension; ++j) {
                    if (*pos++) res += (degrees.buf[j]);
                }
                set_ww_and_shift(ww, res, bit_start, bit_end, current_byte, split_byte);
            #endif
        }
        if (preferred_fast[sn].find(pf) == preferred_fast[sn].end()) {
            *h1p() |= 1;
        }
    }
}

point::point(const point &p, const point_fast& v, SECTOR ssector, bool not_preferred) {
    #ifndef SMALL_POINT
        reinterpret_cast<uint64_t *>(this)[2] = reinterpret_cast<const uint64_t *>(p.ww)[2];
        reinterpret_cast<uint64_t *>(this)[1] = reinterpret_cast<const uint64_t *>(p.ww)[1];
        reinterpret_cast<uint64_t *>(this)[0] = reinterpret_cast<const uint64_t *>(p.ww)[0];
    #else
        reinterpret_cast<uint128_t *>(this)[0] = 0;
        *h1p() = p.h1();
    #endif

    auto ordering_now = common::orderings_fast[ssector].get();
    SECTOR bit = 1<<(common::dimension-1);
    t_index *pos = ordering_now;
    #ifndef SMALL_POINT
        for (unsigned int j = 0; j != common::dimension; ++j) {
            if (v.buf[j]) {
                t_index shift = (bit & ssector) ? v.buf[j] : -v.buf[j];
                for (unsigned int i = 0; i != common::dimension; ++i) {
                    if (*pos) ww[MAX_IND - i - 1] += shift;
                    pos += common::dimension;
                }
                pos -= common::dimension*common::dimension;
            }
            bit >>= 1;
            ++pos;
        }
    #else
        char temp_ww[MAX_IND];
        {
            unsigned short bit_start = {0};
            unsigned short bit_end = {BITS_PER_INDEX - 1};
            unsigned short current_byte = {0};
            bool split_byte = false;
            for (int i = 0; i!=common::dimension; ++i) {
                temp_ww[i] = get_ww_and_shift(p.ww, bit_start, bit_end, current_byte, split_byte);
            }
        }
        for (unsigned int j = 0; j != common::dimension; ++j) {
            if (v.buf[j]) {
                t_index shift = (bit & ssector) ? v.buf[j] : -v.buf[j];
                for (unsigned int i = 0; i != common::dimension; ++i) {
                    if (*pos) temp_ww[i] +=shift;
                    pos += common::dimension;
                }
                pos -= common::dimension*common::dimension;
            }
            bit >>= 1;
            ++pos;
        }
        {
            unsigned short bit_start = {0};
            unsigned short bit_end = {BITS_PER_INDEX - 1};
            unsigned short current_byte = {0};
            bool split_byte = false;
            for (int i = 0; i!=common::dimension; ++i) {
                set_ww_and_shift(ww, temp_ww[i], bit_start, bit_end, current_byte, split_byte);
            }
        }
    #endif

    *h1p() &= ~1u;
    if (not_preferred) {
        *h1p() |= 1;
    }

}


point::point(char *buf) {
    for (unsigned int i = 0; i != sizeof(point); ++i) {
        ww[i] = (buf[i] << 4) ^ (buf[i + sizeof(point)] & 15);
    }
}

void point::safe_string(char *buf) const {
    for (unsigned int i = 0; i != sizeof(point); ++i) {
        buf[i] = (ww[i] >> 4) | 64;
        buf[i + sizeof(point)] = (ww[i] & 15) | 64;
    }
}


string point::number() const {
    stringstream ss(stringstream::out);
    ss << h1();
    #ifdef SMALL_POINT
        unsigned short bit_start = {0};
        unsigned short bit_end = {BITS_PER_INDEX - 1};
        unsigned short current_byte = {0};
        bool split_byte = false;
    #endif
    for (int i = 0; i != common::dimension; ++i) {
        #ifndef SMALL_POINT
            ss << std::setfill('0') << std::setw(3) << (int(static_cast<unsigned char>(ww[MAX_IND - i - 1])) - 1);
        #else
            unsigned short sn = s_number();
            vector<t_index> result;
            if ((sn == 1) || (sn == common::virtual_sector)) {
                ss << std::setfill('0') << std::setw(3) << (int(static_cast<unsigned char>(
                    get_ww_and_shift(ww, bit_start, bit_end, current_byte, split_byte) - 1 - (1u<<(BITS_PER_INDEX-1)) +128)
                ));
            } else
                ss << std::setfill('0') << std::setw(3) << (int(static_cast<unsigned char>(get_ww_and_shift(ww, bit_start, bit_end, current_byte, split_byte) - 1)));
        #endif
    }
    return ss.str();
}

bool point::is_zero() const {
    #ifndef SMALL_POINT
        return  (reinterpret_cast<const uint64_t *>(this)[2] == 0) &&
            (reinterpret_cast<const uint64_t *>(this)[1] == 0) &&
            (reinterpret_cast<const uint64_t *>(this)[0] == 0);
    #else
        return reinterpret_cast<const uint128_t *>(this)[0] == 0;
    #endif
}

/**
 * Update output string stream, basically printing the point
 * @param out output ostream
 * @param p point to be printed
 * @return reference to updated output stream
 */
ostream &operator<<(ostream &out, const point &p) {
    out << (point::print_g ? "G[" : "{") << common::global_pn << ",{";
    if (p.virt()) {
        out << p.s_number() << ",";
        out << reinterpret_cast<const virt_t *>(p.ww + sizeof(point) - 7)[0];
    } else {
        unsigned int it = 0;
        vector<t_index> v = p.get_vector();
        out << int(v[it]);
        for (it++; (it != v.size()); it++) out << "," << int(v[it]);
    }
    out << "}" << (point::print_g ? "]" : "}");
    return out;
}

vector<t_index> point::get_vector() const {
    unsigned short sn = s_number();
    vector<t_index> result;
    if ((sn == 1) || (sn == common::virtual_sector)) {
        #ifdef SMALL_POINT
            unsigned short bit_start = {0};
            unsigned short bit_end = {BITS_PER_INDEX - 1};
            unsigned short current_byte = {0};
            bool split_byte = false;
        #endif
        for (int i = 0; i != common::dimension; i++) {
            #ifndef SMALL_POINT
                result.push_back(ww[MAX_IND - i - 1] - 128);
            #else
                char res = get_ww_and_shift(ww, bit_start, bit_end, current_byte, split_byte);
                result.push_back(res - (1u<<(BITS_PER_INDEX-1)));
            #endif
        }
    } else {
        vector<vector<t_index> > &iordering = common::iorderings[sn];
        vector<t_index> ssector = common::ssectors[sn];
        for (int i = 0; i != common::dimension; ++i) {
            #ifdef SMALL_POINT
                unsigned short bit_start = {0};
                unsigned short bit_end = {BITS_PER_INDEX - 1};
                unsigned short current_byte = {0};
                bool split_byte = false;
            #endif
            t_index z = 0;
            if (!virt()) {
                for (int j = 0; j != common::dimension; ++j) {
                    #ifndef SMALL_POINT
                        if (iordering[i][j] == 1) z += (ww[MAX_IND - j - 1] - 1);
                        else if (iordering[i][j] == -1) z -= (ww[MAX_IND - j - 1] - 1);
                    #else
                        char val = get_ww_and_shift(ww, bit_start, bit_end, current_byte, split_byte);
                        if (iordering[i][j] == 1) z += (val - 1);
                        else if (iordering[i][j] == -1) z -= (val - 1);
                    #endif
                }
            }
            if (ssector[i] == 1) z++; else z = -z;
            result.push_back(z);
        }
    }
    return result;
}


set<point_fast> level_points_fast(point_fast s, const unsigned int pos, const unsigned int neg) {
    if ((pos == 0) && (neg == 0)) {
        set<point_fast> r;
        r.insert(s);
        return (r);
    }
    if (neg > 0) {
        set<point_fast> old = level_points_fast(s, pos, neg - 1);
        set<point_fast> r;
        for (const auto & p : old) {
            for (unsigned int i = 0; i < common::dimension; ++i) {
                if (p.buf[i] <= 0) {
                    point_fast p2 = p;
                    p2.buf[i]--;
                    r.insert(p2);
                }
            }
        }
        return r;
    }
    set<point_fast> old = level_points_fast(s, pos - 1, neg);
    set<point_fast> r;
    for (const auto & p : old) {
        for (unsigned int i = 0; i < common::dimension; ++i) {
            if (p.buf[i] > 0) {
                point_fast p2 = p;
                p2.buf[i]++;
                r.insert(p2);
            }
        }
    }
    return r;
}

//=============================================================
// Memory DataBase
//=============================================================
inline bool file_exists(string fn) {
    return (access(fn.c_str(),F_OK)!=-1);
}
inline void point_output(const point & p, mpz_t & z, ostringstream &oss, void (*free_mp_str)(void *, size_t)) {
    #ifndef SMALL_POINT
    int wwn = 24;
    #else
    int wwn = 16;
    #endif
    auto & ww = p.ww;
    mpz_import(z, wwn, 1, sizeof(ww[0]), 0, 0, ww);
    char * str = mpz_get_str(NULL, 62, z);
    oss << str << endl;
    free_mp_str(str, strlen(str)+1);
}
inline void db_export(unsigned int sn) {
    auto & dbn = DBM[sn];
    ostringstream oss;
    // status
    oss << dbn.status << endl;
    // lower_size
    oss << dbn.lower_size << endl;
    mpz_t z;
    mpz_init(z);
    void (*free_mp_str)(void *, size_t);
    mp_get_memory_functions(NULL,NULL,&free_mp_str);
    // vector<point> lower;
    for(int i=0; i<dbn.lower_size; i++) point_output(dbn.lower[i],z,oss,free_mp_str);
    // map<pair<unsigned int, unsigned int>, bool> umap;
    oss << dbn.umap.size() << endl;
    for(auto const & kv : dbn.umap) {
        oss << kv.first.first << " " << kv.first.second << " " << kv.second << endl;
    }
    // map<point, pair<unsigned char, pc_pair_vec>> pmap;
    oss << dbn.pmap.size() << endl;
    for(auto & kv : dbn.pmap) {
        point_output(kv.first,z,oss,free_mp_str);
        int c = kv.second.first; // char to int
        oss << c << endl;
        oss << kv.second.second.size() << endl;
        for(auto & pcr : kv.second.second) {
            point_output(pcr->first,z,oss,free_mp_str);
            pcr->second.export_to(oss, 62);
        }
    }
    oss << "#DONE#" << endl;
    string oss_str = oss.str();
    fstream db(common::path+"/"+to_string(sn)+".map", fstream::out|std::ios::binary);
    db.write(oss_str.c_str(), oss_str.size());
    db.close();
    mpz_clear(z);
}

inline void str_z_point(const string &s, mpz_t &z, point &p) {
    #ifndef SMALL_POINT
    int wwn = 24;
    #else
    int wwn = 16;
    #endif
    mpz_set_str(z, s.c_str(), 62);
    int numb = 8*sizeof(p.ww[0]);
    int count = (mpz_sizeinbase(z,2)+numb-1)/numb;
    mpz_export(p.ww+(wwn-count), NULL, 1, sizeof(p.ww[0]), 0, 0, z);
}
inline void db_import(unsigned int sn, int mode, const set<point> *pset) {
    mpz_t z;
    mpz_init(z);
    auto & dbn = DBM[sn];
    fstream db(common::path+"/"+to_string(sn)+".map", fstream::in);
    size_t n;
    // status
    db >> dbn.status;
    // lower_size
    db >> dbn.lower_size;
    // vector<point> lower;
    dbn.lower.clear();
    dbn.lower.reserve(dbn.lower_size);
    for(int i=0; i<dbn.lower_size; i++) {
        string s;
        db >> s;
        point p;
        str_z_point(s,z,p);
        dbn.lower.push_back(p);
    }
    if(mode==0) {
        auto uend(end(dbn.umap));
        db >> n;
        // map<pair<unsigned int, unsigned int>, bool> umap;
        for(int i=0; i<n; i++) {
            unsigned int k1, k2;
            bool b;
            db >> k1 >> k2 >> b;
            uend = dbn.umap.insert(uend,make_pair(make_pair(k1,k2),b));
        }
        auto pend(end(dbn.pmap));
        db >> n;
        // map<point, pair<unsigned char, pc_pair_ptr_vec>> pmap;
        for(int i=0; i<n; i++) {
            string s;
            db >> s;
            point p;
            str_z_point(s,z,p);
            
            bool needed = true;
            if(pset!=NULL && pset->find(p)==pset->end()) needed = false;
            int ci; char c;
            db >> ci; c = ci;
            size_t nn;
            db >> nn;
            pc_pair_ptr_vec vec;
            if(needed) vec.reserve(nn);
            for(int j=0; j<nn; j++) {
                string s;
                db >> s; // point
                if(needed) {
                    point pp;
                    str_z_point(s,z,pp);
                    COEFF c;
                    c.import_from(db,62);
                    vec.push_back(make_pc_ptr(pp,c));
                } else {
                    for(int l=0; l<COEFF::io_lines; l++) db >> s; // need to read and skip
                }
            }
            if(needed) pend = dbn.pmap.insert(pend,make_pair(p, make_pair(c, std::move(vec))));
        }
        string ok;
        db >> ok;
        if(ok!="#DONE#") {
            cout << "expect #DONE#, but get '" << ok << "'" << endl;
            cout << "db_import failed @ " << sn << ".map" << endl;
            abort();
        }
    }
    db.close();
    mpz_clear(z);
}

map<int, DB> DBM;
static mutex db_mutex;
static set<int> dbs;
static map<int,mutex> dbn_mutex;

int read_status(int sn) {
    if(true) {
        lock_guard<mutex> guard(db_mutex);
        auto itr = dbs.find(sn);
        if(itr == dbs.end()) {
            dbs.insert(sn);
            dbn_mutex[sn];
            DBM[sn];
        }
    }
    string fn = common::path+"/"+to_string(sn)+".map";
    if(!file_exists(fn)) return 0;
    int status;
    fstream db(fn, fstream::in);
    db >> status;
    db.close();
    return status;
}

bool database_exists(int number) {
    lock_guard<mutex> guard(db_mutex);
    auto itr = dbs.find(number);
    if(itr != dbs.end()) return true;
    else return false;
}

void open_database(int number, int mode, const set<point> *pset) {
    if(!common::run_sector) {
        lock_guard<mutex> guard(db_mutex);
        dbn_mutex[number];
        auto & dbn = DBM[number];
        auto itr = dbs.find(number);
        if(itr == dbs.end()) {
            dbs.insert(number);
            if(common::run_mode!=2) {
                dbn.opened++; // note ++
                return;
            }
        }
        if(!file_exists(common::path+"/"+to_string(number)+".map")) {
            dbn.opened++; // note ++
            return;
        }
        // still need to import when run_mode>1
    }
    if(common::run_mode==0) return;
    lock_guard<mutex> guard(dbn_mutex[number]);
    auto & dbn = DBM[number];
    if(dbn.keep) return;
    if(dbn.opened<1) {
        if(!file_exists(common::path+"/"+to_string(number)+".map")) {
            cout << "DB: " << (common::path+"/"+to_string(number)+".map") << " Not Found!" << endl;
            if(common::run_sector) cout << "Please run with '-sector 0' first to generate the DB files." << endl;
            abort();
        }
        db_import(number, mode, pset);
        dbn.open_mode = mode;
        dbn.need_write = false; // initial to false
    } else if(mode != dbn.open_mode) { // other thread already call open_database, but different mode
        cout << "something should be wrong in open_database, mode: " << mode << " vs " << dbn.open_mode << endl;
        abort();
    }
    dbn.opened++;
}

// assuming only ONE thread call close for a specific sector
void close_database(int number) {
    if(common::run_mode==0) return;
    lock_guard<mutex> guard(dbn_mutex[number]);
    auto & dbn = DBM[number];
    if(dbn.keep) return;
    if(dbn.need_write && dbn.opened!=1) {
        cout << "dbn.need_write = true &  dbn.opened = " << dbn.opened << endl;
        abort();
    }
    if(dbn.open_mode==1 && dbn.need_write) {
        cout << "something should be wrong in close_database." << endl;
        abort();
    }
    dbn.opened--;
    if(dbn.opened>0) return;
    if(dbn.need_write) db_export(number);
    dbn.clear();
}

void DB::clear() {
    lower_size = 0;
    vector<point> vp;
    lower.swap(vp);
    umap_type um;
    std::swap(umap, um);
    pmap_type pm;
    std::swap(pmap, pm);
}
