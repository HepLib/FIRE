/** @file equation.cpp
 *  @author Alexander Smirnov
 *
 *  This file is a part of the FIRE package.
 */

#include "equation.h"
#include <mutex>
#include <condition_variable>
#include <chrono>

//initialization of static class members
map<vector<t_index>, point> equation::initial;

// get monoms from a point (Feynman integral). database access used
vector<point> p_get_monoms(const point &p, sector_count_t fixed_database_sector) {
    sector_count_t dsector = (fixed_database_sector == 0) ? p.s_number() : fixed_database_sector;
    auto const & pm = DBM[dsector].pmap;
    auto itr = pm.find(p);
    vector<point> result;
    if (itr == pm.end()) return result;
    for(auto const & item : itr->second.second) result.emplace_back(item->first);
    return result;
}

bool p_is_empty(const point &p, sector_count_t fixed_database_sector) {
    sector_count_t dsector = (fixed_database_sector == 0) ? p.s_number() : fixed_database_sector;
    auto const & pm = DBM[dsector].pmap;
    auto itr = pm.find(p);
    if (itr == pm.end()) return true;
    return itr->second.second.size() == 0;
}

void p_get(const point &p, pc_pair_ptr_vec &terms, sector_count_t fixed_database_sector) {
    sector_count_t dsector = (fixed_database_sector == 0) ? p.s_number() : fixed_database_sector;
    auto const & pm = DBM[dsector].pmap;
    auto itr = pm.find(p);
    pc_pair_ptr_vec empty_vec;
    terms.swap(empty_vec);
    if (itr == pm.end()) return;
    
    auto len = itr->second.second.size();
    if(len==0) return;
    
    terms.reserve(len);
    for(auto const & item : itr->second.second) terms.emplace_back(item);
}

void p_get(const point &p, pc_pair_ptr_lst &terms, sector_count_t fixed_database_sector) {
    terms.clear();
    sector_count_t dsector = (fixed_database_sector == 0) ? p.s_number() : fixed_database_sector;
    auto const & pm = DBM[dsector].pmap;
    auto itr = pm.find(p);
    if (itr == pm.end()) return;
    
    auto len = itr->second.second.size();
    if(len==0) return;
    
    for(auto const & item : itr->second.second) terms.emplace_back(item);
}

void p_set(const point &p, pc_pair_ptr_vec && terms, unsigned char level, sector_count_t fixed_database_sector) {
    sector_count_t dsector = (fixed_database_sector == 0) ? p.s_number() : fixed_database_sector;
    auto & dbn = DBM[dsector];
    dbn.need_write = true;
    auto & pm = dbn.pmap;
    auto itr = pm.find(p);
    if (itr == pm.end()) {
        pm.emplace(p, make_pair(level, std::move(terms)));
    } else {
        auto old_level = itr->second.first;
        auto old_n = itr->second.second.size();
        auto n = terms.size();
        
        if ((old_n != 0) && (n == 0)) return;
        if ((old_level == 127) && (n) && (level != 127)) old_level = 126;
        auto new_level = old_level;
        if (level > new_level) new_level = level;
        itr->second.first = new_level;
        itr->second.second.swap(terms);
    }
}

bool is_lower_in_orbit(const vector<t_index> &lhs, const vector<t_index> &rhs) {
    if (lhs == rhs) return false;
    vector<t_index> s1 = _sector_(lhs);
    vector<t_index> s2 = _sector_(rhs);
    if (s1 != s2) {
        sector_count_t sn1 = common::sector_numbers_fast[sector_fast(s1)];
        sector_count_t sn2 = common::sector_numbers_fast[sector_fast(s2)];
        if (sn1 < sn2) return (true);
        if (sn1 > sn2) return (false);
    }

    t_index *ordering_now = common::orderings_fast[sector_fast(lhs)].get();
    vector<t_index> d1 = degree(lhs);
    vector<t_index> d2 = degree(rhs);
    unsigned int n = lhs.size();
    for (unsigned int i = 0; i != n; ++i) {
        int pr = 0;
        for (unsigned int j = 0; j != n; ++j) {
            if (ordering_now[i * common::dimension + j] == 1) {
                pr += (d1[j] - d2[j]);
            }
        }
        if (pr < 0) return true;
        if (pr > 0) return false;
    }
    return false;
}


// sort pairs of points and coefficients using only points
bool pair_point_COEFF_smaller(const pc_pair &lhs, const pc_pair &rhs) {
    return lhs.first < rhs.first;
}


bool point_fast_smaller_in_sector(const point_fast & pf1, const point_fast & pf2, SECTOR s) {
    if (pf1 == pf2) return false;
    t_index *ordering_now = common::orderings_fast[s].get();
    point_fast d1 = pf1.degree();
    point_fast d2 = pf2.degree();

    for (unsigned int i = 0; i != common::dimension; ++i) {
        int pr = 0;
        for (unsigned int j = 0; j != common::dimension; ++j) {
            if (ordering_now[i * common::dimension + j] == 1) pr += (d1.buf[j] - d2.buf[j]);
        }
        if (pr < 0) return true;
        if (pr > 0) return false;
    }
    cout << "impossible point compare" << endl;
    for (unsigned int i = 0; i != common::dimension; ++i) cout << int(pf1.buf[i]) << ";";
    cout << endl;
    for (unsigned int i = 0; i != common::dimension; ++i) cout << int(pf2.buf[i]) << ";";
    cout << endl;
    cout << "impossible point compare" << endl;
    abort();// this should not happen
}


// point reference version without std
point point_reference_fast(const point_fast &v) {
    SECTOR ssector = v.sector_fast();
    sector_count_t sn = common::sector_numbers_fast[ssector];
    if (sn == 0) {
        return point();
    }
    if (sn == 1) {
        return point(v, ssector);
    }

    if (common::symmetries.size() > 1) {  //there are symmetries
        vector<vector<vector<t_index> > > &sym = common::symmetries;
        point_fast best = v;
        SECTOR best_sector = ssector;
        sector_count_t best_sn = sn;
        for (const auto &values : sym) {
            const vector<t_index> &permutation = values[0];
            point_fast p_new;

            for (unsigned int i = 0; i != common::dimension; ++i) { p_new.buf[i] = v.buf[permutation[i] - 1]; }

            // we only use the first part of symmetries, but I do not even know whether the other parts used to work properly
            // part 3 can be added only at parser time and means something related to part 2 = conditional symmetries
            // part 1 i odd symmetries, but they did not work properly even in earlier versions

            // now we need to compare the points and choose the lowest
            // best_sector is either common::virtual_sector or some good sector;

            SECTOR new_sector = p_new.sector_fast();
            sector_count_t new_sn = common::sector_numbers_fast[new_sector];
            if (new_sn == common::virtual_sector) {
                continue; // it is a higher virtual point
            }

            if (best_sn == common::virtual_sector) { // the old one was virtual, but now a real point comes
                best_sn = new_sn;
                best_sector = new_sector;
                best = p_new;
                continue;
            }

            // here we are left with the case when both new and old are virtual; this means they are in the same sector
            if (point_fast_smaller_in_sector(p_new, best, best_sector)) {
                best = p_new;
            }
        }
        return point(best, best_sector);
    } else {
        return point(v, ssector);
    }
}


// get the right symmetry point by the vector of coordinates
point point_reference(const vector<t_index> &v) {
    SECTOR ssector = sector_fast(v);
    sector_count_t sn = common::sector_numbers_fast[ssector];
    if (sn == 0) {
        return point();
    }
    if (sn == 1) {
        return point(v, 0, ssector);
    }

    if (common::symmetries.size() > 1) {  //there are symmetries
        vector<vector<vector<t_index> > > &sym = common::symmetries;
        vector<vector<t_index> > orbit;
        Orbit(v, orbit, sym);
        vector<t_index> *lowest = &(*orbit.begin());
        auto itr = orbit.begin();
        itr++;
        while (itr != orbit.end()) {
            if (is_lower_in_orbit(*itr, *lowest)) {
                lowest = &(*itr);
            }
            itr++;
        }
        return point(*lowest);
    } else {
        return point(v, 0, ssector);
    }
}


