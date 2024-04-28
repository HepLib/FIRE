#pragma once

#include "point.h"
#include "gateToFermat.h"
#include <semaphore.h>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <chrono>

//typedef pc_pair_ptr MONOM;
class equation {
public:
    static map<vector<t_index>, point> initial;

    explicit equation(unsigned int n) : terms(n), length(0) { }
    equation() : equation(0) { }

    vector<pc_pair_ptr> terms;
    pair<point_fast, unsigned short> source;

    unsigned int length{};
};

bool pair_point_COEFF_smaller(const pc_pair &lhs, const pc_pair &rhs);
bool p_is_empty(const point &, sector_count_t fixed_database_sector = 0);
void p_get(const point &p, pc_pair_ptr_vec &terms, sector_count_t fixed_database_sector = 0);
void p_get(const point &p, pc_pair_ptr_lst &terms, sector_count_t fixed_database_sector = 0);

void p_set(const point &p, pc_pair_ptr_vec &&terms, unsigned char level, sector_count_t fixed_database_sector = 0);
inline void p_set(const point &p, pc_pair_vec &&terms, unsigned char level, sector_count_t fixed_database_sector = 0) {
    pc_pair_ptr_vec ptr_terms;
    ptr_terms.reserve(terms.size());
    for(auto & item : terms) ptr_terms.emplace_back(make_pc_ptr(std::move(item)));
    p_set(p, std::move(ptr_terms), level, fixed_database_sector);
}

inline void p_set(const point &p, unsigned int n, pc_pair_ptr_lst::iterator termsB, pc_pair_ptr_lst::iterator termsE, unsigned char level, sector_count_t fixed_database_sector) {
    pc_pair_ptr_vec terms;
    terms.reserve(n);
    for(auto itr=termsB; itr!=termsE; ++itr) terms.emplace_back(*itr);
    p_set(p, std::move(terms), level, fixed_database_sector);
}
inline void p_set(const point &p, pc_pair_ptr_lst &terms, unsigned char level, sector_count_t fixed_database_sector = 0) {
    p_set(p, static_cast<unsigned int>(terms.size()), terms.begin(), terms.end(), level, fixed_database_sector);
}

vector<point> p_get_monoms(const point &p, sector_count_t fixed_database_sector = 0);
point point_reference(const vector<t_index> &v);
point point_reference_fast(const point_fast &v);
bool point_fast_smaller_in_sector(const point_fast & pf1, const point_fast & pf2, SECTOR s);
bool is_lower_in_orbit(const vector<t_index> &lhs, const vector<t_index> &rhs);

struct indirect_more {
    bool operator()(const point &lhs, const point &rhs) const {
        return (rhs) < (lhs);
    }
};
