/** @file functions.cpp
 * @author Alexander Smirnov
 *
 * This file is a part of the FIRE package
*/

#include "functions.h"
#include "common.h"
#include "parser.h"
#include <algorithm>
#include <sys/stat.h>
#include <sys/wait.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <functional>

#if !defined(__APPLE__)
#include <malloc.h>
#endif

inline char digit2char(const int i) {
    if (0 <= i && i <= 9) {
        return static_cast<char>('0' + i);
    }
    cout << "Error in digit2char" << endl;
    abort();
}

inline unsigned int char2digit(const char i) {
    if ('0' <= i && i <= '9') {
        return static_cast<unsigned int>(i - '0');
    }
    cout << "Error in char2digit" << endl;
    abort();
}

void subs2(string &str, const vector<t_index> &v) {
    size_t length = str.size();
    for (unsigned int i = 0; i != length; ++i) {
        if ((i + 2) >= length) { continue; }
        if (!((str[i] == 'a') && (str[i + 1] == 'a'))) { continue; }
        unsigned int ind = 0;
        if ((str[i] == 'a') &&
            (str[i + 1] == 'a') &&
            (str[i + 2] == 'a') &&
            (str[i + 3] == '(') &&
            (str[i + 5] == ')')) {
            ind = char2digit(str[i + 4]);
        } else if ((str[i] == 'a') &&
                   (str[i + 1] == 'a') &&
                   (str[i + 2] == 'a') &&
                   (str[i + 3] == '(') &&
                   (str[i + 6] == ')')) {
            ind = (10 * char2digit(str[i + 4])) + char2digit(str[i + 5]);
            str[i] = ' ';
            i++;
        }
        if (ind == 0) {
            cout << "subs2 error" << endl;
            cout << str << endl;
            cout << i << endl;
            abort();
        }
        ind--;
        str[i] = '(';
        if (v[ind] >= 100) {
            str[i + 1] = ' ';
            str[i + 2] = digit2char(v[ind] / 100);
            str[i + 3] = digit2char((v[ind] % 100) / 10);
            str[i + 4] = digit2char(v[ind] % 10);
        } else if (v[ind] >= 10) {
            str[i + 1] = ' ';
            str[i + 2] = ' ';
            str[i + 3] = digit2char(v[ind] / 10);
            str[i + 4] = digit2char(v[ind] % 10);
        } else if (v[ind] >= 0) {
            str[i + 1] = ' ';
            str[i + 2] = ' ';
            str[i + 3] = ' ';
            str[i + 4] = digit2char(v[ind]);
        } else if (v[ind] <= -100) {
            str[i + 1] = '-';
            str[i + 2] = digit2char(-v[ind] / 100);
            str[i + 3] = digit2char((-v[ind] % 100) / 10);
            str[i + 4] = digit2char(-v[ind] % 10);
        } else if (v[ind] <= -10) {
            str[i + 1] = ' ';
            str[i + 2] = '-';
            str[i + 3] = digit2char(-v[ind] / 10);
            str[i + 4] = digit2char(-v[ind] % 10);
        } else {
            str[i + 1] = ' ';
            str[i + 2] = ' ';
            str[i + 3] = '-';
            str[i + 4] = digit2char(-v[ind]);
        }
        i += 5;
    }
}

bool MONOM_smaller(const pc_pair_ptr & lhs, const pc_pair_ptr  & rhs) {
    return (lhs->first) < (rhs->first);
}


equation apply(const vector<pair<vector<COEFF>, point_fast > >& ibp, point_fast v, const SECTOR ssector_fast) {
    equation result(ibp.size());
    unsigned int length = 0;
    sector_count_t sn = common::sector_numbers_fast[ssector_fast];
    point p_orig(v,ssector_fast);

    bool there_are_symmetries = (common::symmetries.size() > 1);

    for (auto read = ibp.begin(); read != ibp.end(); ++read) {
        point_fast new_v = read->second + v;
        SECTOR s1 = new_v.sector_fast();

        if ((s1 != ssector_fast) && (~(ssector_fast | (~s1)))) { continue; }
        //bitwise not to s1, or ssector_fast is full of 1 if sector_fast is over s1. now again bitwise not, and it should become 0 only if sector_fast is over s1

        point new_p;

        if (!there_are_symmetries && (s1 == ssector_fast))  {
            new_p = point(p_orig,read->second,ssector_fast, point::preferred_fast[sn].find(new_v) == point::preferred_fast[sn].end());
        } else {
            new_p = point_reference_fast(new_v);
        }

        if (new_p.is_zero()) { continue; }
        COEFF o1, o2;
        t_index *pos = v.buf;
        _set_(o1, (read->first)[0]);
        for (unsigned short i = 0; i != common::dimension; ++i, ++pos) {
            _set_(o2, int(*pos));
            _add_mul_(o1, o2, (read->first)[i+1]);
        }
        if (!o1.is_zero()) {
            result.terms[length] = make_pc_ptr(new_p, std::move(o1));
            ++length;
        }
    }

    result.length = length;
    sort(result.terms.begin(), result.terms.begin() + length, MONOM_smaller);

    unsigned int k = 0; // k is where we write to; i is where we read from;
    for (unsigned int i = 0; i != result.length; ++i, ++k) {
        if (k != i) {
            result.terms[k] = result.terms[i]; // just changing the pointer to another portion of data
        }

        while (((i + 1) != result.length) && (result.terms[i + 1]->first == result.terms[k]->first)) {
            // the next term is equal
            COEFF o1;
            _add_(o1, result.terms[k]->second, result.terms[i + 1]->second);
            result.terms[k] = make_pc_ptr(result.terms[k]->first, std::move(o1));
            if (result.terms[k]->second.is_zero()) {
                --k;
                ++i;
                break; // out of the internal cycle. we will anyway increase i and k
            }
            ++i;
        }
    }
    result.length -= (result.length - k);

    return result;
}

bool vector_smaller_in_sector(const point_fast &lhs, const point_fast &rhs, SECTOR s) {
    if (lhs == rhs) {
        return false;
    }
    t_index *ordering_now = common::orderings_fast[s].get();

    vector<t_index> d1;
    vector<t_index> d2;
    SECTOR bit = 1u << (common::dimension - 1);
    for (unsigned int i = 0; i != common::dimension; ++i, bit >>= 1) {
        if (s & bit) {
            d1.push_back(lhs.buf[i]);
            d2.push_back(rhs.buf[i]);
        } else {
            d1.push_back(-lhs.buf[i]);
            d2.push_back(-rhs.buf[i]);
        }
    }

    for (unsigned int i = 0; i != common::dimension; ++i) {
        int pr = 0;
        for (unsigned int j = 0; j != common::dimension; ++j) {
            if (ordering_now[i * common::dimension + j] == 1) pr += (d1[j] - d2[j]);
        }
        if (pr < 0) return true;
        if (pr > 0) return false;
    }
    return false;
}

map<point_fast, COEFF> calculate_product(const list<vector<pair<COEFF, point_fast> > > &product, unsigned int n) {
    // let's calculate the product
    map<point_fast, COEFF> local_coeffs;
    for (const auto &v1 : product) {
        map<point_fast, COEFF> new_local_coeffs;
        if (local_coeffs.empty()) {
            for (const auto &i1 : v1) {
                new_local_coeffs.emplace(i1.second, i1.first);
            }
        } else {
            for (const auto &i1 : v1) {
                for (const auto &local_coeff : local_coeffs) {
                    point_fast new_v;
                    for (unsigned short u = 0; u != n; u++) {
                        new_v.buf[u] = i1.second.buf[u] + local_coeff.first.buf[u];
                    }
                    const auto i3 = new_local_coeffs.find(new_v);
                    if (i3 == new_local_coeffs.end()) {
                        COEFF new_c;
                        _mul_(new_c, i1.first, local_coeff.second);
                        new_local_coeffs.emplace_hint(i3, new_v, std::move(new_c));
                    } else {
                        _add_mul_(i3->second, i1.first, local_coeff.second);
                    }
                }
            }
        }
        local_coeffs.clear();
        for (auto &&new_local_coeff : new_local_coeffs) {
            if (!is_zero(new_local_coeff.second)) {
                local_coeffs.insert(new_local_coeff);
            }
        }
    }
    return local_coeffs;
}

/**
 * Writes Lee symmetries with the given complexity level to the database.
 * @param p_start corner of the sector
 * @param pos number of dots
 * @param neg number of numerators
 * @param sector_level number of positive indices in this sector
 * @return number of symmetries written
 */
int write_symmetries(const point &p_start, const unsigned int pos, const unsigned int neg, int sector_level) {
    int result = 0;
    int n = p_start.get_vector().size();
    auto iitr = point::ibases.find(p_start.s_number());
    if (iitr != point::ibases.end()) {
        point_fast p_fast = point_fast(p_start);
        set<point_fast> s = level_points_fast(p_fast, pos, neg);
        if (neg == 1) {
            set<point_fast> s2 = level_points_fast(p_fast, pos, 0);
            for (const auto &pp : s2) s.insert(pp);
        }
        if (pos == 1) {
            set<point_fast> s2 = level_points_fast(p_fast, 0, neg);
            for (const auto &pp : s2) s.insert(pp);
        }
        if ((neg == 1) && (pos == 1)) {
            set<point_fast> s2 = level_points_fast(p_fast, 0, 0);
            for (const auto &pp : s2) s.insert(pp);
        }
        for (const auto &pp : s) {
            // go through points
            // no need for symmetries here more
            // y is pp.buf
            point p = point_reference_fast(pp);
            if (!p_is_empty(p)) {
                continue;
            }

            for (const auto &ibasis : iitr->second) {
                // go through symmetries

                list<vector<pair<COEFF, point_fast> > > product;
                point_fast new_pp;
                for (int j = 0; j != n; ++j) {  // here we use the permutation
                    if (ibasis.first[j] == 0) {
                        new_pp.buf[j] = 0;
                    } else {
                        new_pp.buf[j] = pp.buf[ibasis.first[j] - 1];
                    }
                }
                vector<pair<COEFF, point_fast> > new_term;
                new_term.emplace_back(CO_1, new_pp);
                product.push_back(new_term);  // starting a product, next parts will be in a cycle

                for (const auto &pitr : ibasis.second) {
                    int power = -pp.buf[pitr.second - 1];
                    if (power < 0) {
                        cout << "Wrong internal symmetry rule: " << p << endl;
                        abort();
                    }
                    for (int j = 0; j != power; j++) {
                        product.push_back(pitr.first);
                    }
                }
                map<point_fast, COEFF> local_coeffs = calculate_product(product, n); // this is a shared function in external and internal symmetry generation

                // product calculated, now it is time to get the rule
                // let's add -p right here
                auto itr = local_coeffs.find(pp);
                if (itr == local_coeffs.end()) {
                    local_coeffs.emplace(pp, CO_1m);
                } else {
                    _sub_(itr->second, itr->second, CO_1);
                }

                // now we will convert here to strings, but...
                pc_pair_vec mon;
                mon.reserve(local_coeffs.size());

                for (const auto &local_coeff : local_coeffs) {
                    if (!is_zero(local_coeff.second))
                    {
                        point new_p = point_reference_fast(local_coeff.first);
                        if (new_p.is_zero()) continue;
                        if (new_p.s_number() == 1) { // it can send to sector 1???
                            cout << "Incorrect internal symmetry for " << p << endl;
                            cout << "Sending to " << new_p << endl;
                            abort();
                        }

                        if ((new_p.s_number() != p.s_number())) { // internal symmetry should not map to another sector of same level?
                            if (new_p.level() >= p.level()) {
                                cout << "Incorrect internal symmetry for " << p << endl;
                                cout << "Sending to " << new_p << endl;
                                abort();
                            }
                        }
                        mon.emplace_back(new_p, local_coeff.second);
                    }
                }

                // product ready, now sorting, joining, evaluating
                sort(mon.begin(), mon.end(), sort_pair_point_coeff_by_point);
                mon = group_equal_in_sorted(mon);
                if ((mon.empty()) || mon[mon.size() - 1].first != p) {
                    continue;
                    // trivial or sending higher symmetries are simply ignored
                }
                p_set(p, std::move(mon), 2 * sector_level);
                ++result;
            }
        }
    }
    return result;
}


point_fast lowest_in_sector_orbit_fast(const point_fast &p, SECTOR s, const vector<vector<vector<t_index> > > &sym) {
    point_fast result = p;
    for (const auto &values : sym) {
        const vector<t_index> &permutation = values[0];
        point_fast p_new;

        for (unsigned int i = 0; i != common::dimension; ++i) { p_new.buf[i] = p.buf[permutation[i] - 1]; }

        // we only use the first part of symmetries, but I do not even know whether the other parts used to work properly
        // part 3 can be added only at parser time and means something related to part 2 = conditional symmetries
        // part 1 i odd symmetries, but they did not work properly even in earlier versions

        // now we need to compare the points and choose the lowest
        if ((p_new.sector_fast() == s) && point_fast_smaller_in_sector(p_new, result, s)) {
            result = p_new;
        }
    }
    return result;
}


/**
 * Comparison functor - to compare ibps applied in a sector according to the current ordering
 * Assumes that every ibp is internally sorted with first being biggest
 */
struct ibp_comparator {
    SECTOR s; ///< sector passes to the comparator
    /**the internal comparator; compares highest term shifts from the corner
    * @param lhs first vector
    * @param rhs first vector
    * @return compare result, true if second is smaller
    */
    bool operator()(const vector<pair<vector<COEFF>, point_fast > >& lhs, const vector<pair<vector<COEFF>, point_fast > >& rhs) {
        auto& v1 = lhs[0].second;
        auto& v2 = rhs[0].second;
        if (vector_smaller_in_sector(v2,v1,s)) return true;
        if (vector_smaller_in_sector(v1,v2,s)) return false;
        return lhs.size()<rhs.size();
    }
};
//thread_local ibp_comparator ibpcompare; // Instance of ibp_comparator.

unsigned int sort_ibps(const point &p, const set<pair<unsigned int, unsigned int> > &current_levels,
                       vector<pair<point, pair<point_fast, unsigned short> > > &ibps_vector,
                       const vector<point_fast> &IBPdegree,
                       const vector<point_fast> &IBPdegreeFull) {

    point_fast p_fast(p);
    SECTOR sector = p_fast.sector_fast();
    set<point_fast> s_fast;
    for (const auto &current_level : current_levels) {
        const unsigned int pos = current_level.first;
        const unsigned int neg = current_level.second;
        set<point_fast> s0 = level_points_fast(p_fast, pos, neg);
        for (const auto &fp : s0) { s_fast.insert(fp); }
        if (neg == 1) {
            set<point_fast> s2 = level_points_fast(p_fast, pos, 0);
            for (const auto &fp : s2) s_fast.insert(fp);
        }
        if (pos == 1) {
            set<point_fast> s2 = level_points_fast(p_fast, 0, neg);
            for (const auto &fp : s2) s_fast.insert(fp);
        }
        if ((neg == 1) && (pos == 1)) {
            set<point_fast> s2 = level_points_fast(p_fast, 0, 0);
            for (const auto &fp : s2) s_fast.insert(fp);
        }
    }

    if (common::symmetries.size() > 1) {  //there are symmetries
        vector<vector<vector<t_index> > > &sym = common::symmetries;
        set<point_fast> s_new_fast;

        // now taking only lowest in orbits
        for (const auto &fp : s_fast) {
            s_new_fast.insert(lowest_in_sector_orbit_fast(fp, sector, sym));
        }
        s_fast = s_new_fast;
    }

    ibps_vector.reserve(s_fast.size() * IBPdegree.size());
    unsigned int counter = 0;

    for (const auto &fp : s_fast) {
        point_fast deg = fp.degree();
        for (unsigned int i = 0; i != IBPdegree.size(); ++i) {
            pc_pair_vec terms;
            point_fast highest_fast;
            for (unsigned short j = 0; j != common::dimension; ++j) {
                highest_fast.buf[j] = fp.buf[j] + IBPdegreeFull[i].buf[j];
            }
            SECTOR bit = 1u << (common::dimension - 1);
            for (unsigned int j = 0; j < common::dimension; ++j, bit >>= 1) {
                if (sector & bit) {
                    if (highest_fast.buf[j] <= 0) highest_fast.buf[j] = 1;
                } else {
                    if (highest_fast.buf[j] > 0) highest_fast.buf[j] = 0;
                }
            }
            ibps_vector.emplace_back(point(highest_fast), make_pair(fp, static_cast<unsigned short>(i)));
            ++counter;
            if (!common::all_ibps) {
                if (over_fast(deg, IBPdegree[i])) break;
            }
        }
    }
    sort(ibps_vector.begin(), ibps_vector.end(), [](const pair<point, pair<point_fast, unsigned short> > &a,
                                                    const pair<point, pair<point_fast, unsigned short> > &b) -> bool {
        return a.first < b.first;
    });
    return counter;
};


vector<pair<unsigned int, unsigned int> > under_levels(const unsigned int p0, const unsigned int m0) {
    unsigned int p = p0;
    unsigned int m = m0;
    vector<pair<unsigned int, unsigned int> > result;
    while ((m > 0) || (p > 0)) {
        while (p > 0) {
            if ((p <= p0) && (m <= m0) && (m > 0)) {
                result.push_back(make_pair(p,m));
            }
            p--;
            m++;
        }
        // p is now 0, m is positive
        p = m - 1;
        m = 0;
    }
    return result;
}

void make_master(const point &p) {
    if (common::split_masters) {
        cout<<"Unspecified master integral in masters mode: "<<p<<endl;
        abort();
    }

    vector<t_index> v2 = p.get_vector();

    point p2(v2, 0, -2);
    pc_pair_ptr_vec t;
    if (!common::split_masters) t.emplace_back(make_pc_ptr(p2,CO_1));
    t.emplace_back(make_pc_ptr(p,CO_1m));
    p_set(p, std::move(t), 127, p.s_number());
}


void mark_master_integrals(const point &Corner, const unsigned int pos, const unsigned int neg) {
    SECTOR s = point_fast(Corner).sector_fast();
    if ((pos > 0) && (neg > 0)) {
        set<point_fast> p_refs_fast = level_points_fast(point_fast(Corner), pos - 1, neg - 1);
        // now taking lower orbit values
        if (common::symmetries.size() > 1) {  //there are symmetries
            vector<vector<vector<t_index> > > &sym = common::symmetries;
            set<point_fast> p_refs_fast_new;
            // now taking only lowest in orbits
            for (const auto &fp : p_refs_fast) {
                point_fast lowest = lowest_in_sector_orbit_fast(fp, s, sym);
                p_refs_fast_new.insert(lowest);
            }
            p_refs_fast = p_refs_fast_new;
        }
        for (const auto &fp : p_refs_fast) {
            point p = point(fp, s);
            if (p_is_empty(p)) make_master(p);
        }
    }
}


/**
 * Functor for comparing levels of points (using dots and numerators).
 */
struct level_smaller {
    /**the internal comparator; levels are compared by their sums (total shift from the corner), in case of equal the priority is to have less dots
    * @param lhs first level
    * @param rhs first level
    * @return compare result, true if first is smaller
    */
    bool operator()(const pair<unsigned int, unsigned int> &lhs, const pair<unsigned int, unsigned int> &rhs) const {
        if (lhs.first + lhs.second < rhs.first + rhs.second) return true;
        if (lhs.first + lhs.second > rhs.first + rhs.second) return false;
        return lhs.first < rhs.first;
    }
};


void add_needed(map<sector_count_t, set<point> > &needed_lower, const point &p) {
    sector_count_t new_sector = p.s_number();
    auto itr = needed_lower.find(new_sector);
    if (itr == needed_lower.end()) {
        set<point> s;
        s.insert(p);
        needed_lower.emplace(new_sector, s);
    } else {
        itr->second.insert(p);
    }
}

/**
 * Finish calculations in sector, save results in database and print statistics.
 * @param needed_lower_temp set in integrals at are needed in lower sectors and should be written to corresponding databases
 * @param test_sector sector we work in
 */
void finish_sector(const set<point> &needed_lower_temp, sector_count_t test_sector) {
    auto & dbn = DBM[test_sector];
    if(dbn.status!=0) {
        cout << "finish_sector: dbn.status!=0" << endl;
        abort();
    }
    dbn.status = 1;
    if (!needed_lower_temp.empty()) {
        vector<point> point_vec(needed_lower_temp.size());
        int i=0;
        for (auto sitr = needed_lower_temp.begin(); sitr != needed_lower_temp.end(); ++sitr, ++i) {
            point_vec[i] = *sitr;
        }
        dbn.lower = point_vec;
    }

    // INT64_MAX means override anyway
    dbn.lower_size = needed_lower_temp.size();
    dbn.need_write = true;
    close_database(test_sector);
}

void add_ibps(const COEFF& first_mul, const COEFF& second_mul, const ibp_type& first, const ibp_type& second, const SECTOR sector_fast, ibp_type& result) {
    result.clear();
    result.reserve(first.size()+second.size());
    unsigned int i=0;
    unsigned int j=0;
    while (i+j != first.size()+second.size()) {
        if ((j == second.size()) || (((i != first.size())) && vector_smaller_in_sector(second[j].second,first[i].second,sector_fast))) {
            // second does not exist or i is bigger
            vector<COEFF> mult = first[i].first;
            for (COEFF & c : mult) _mul_(c, c, first_mul);
            result.emplace_back(mult,first[i].second);
            ++i;
        } else if ((i == first.size()) || (((j != second.size())) && vector_smaller_in_sector(first[i].second,second[j].second,sector_fast))) {
            // first does not exist or j is bigger
            vector<COEFF> mult = second[j].first;
            for (COEFF &c : mult) _mul_neg_(c, c, second_mul);
            result.emplace_back(mult,second[j].second);
            ++j;
        } else {
            // they are equal
            vector<COEFF> added;
            added.reserve(common::dimension + 1);
            bool has_nonzero = false;
            for (int k = 0; k!=common::dimension + 1; ++k) {
                COEFF c;
                _mul_(c, first_mul, first[i].first[k]);
                _sub_mul_(c, second_mul, second[j].first[k]);
                if (!c.is_zero()) has_nonzero = true;
                added.emplace_back(std::move(c));
            }
            if (has_nonzero) {
                result.emplace_back(added,first[i].second);
            }
            ++i;
            ++j;
        }
    }
}

int matching_index(const vector<COEFF>& first, const vector<COEFF>& second) {
    int matching_index = -1;
    for (int k = 0; k!=common::dimension + 1; ++k) {
        if ((!first[k].is_zero()) || (!second[k].is_zero())) {
            if (matching_index >= 0) {
                matching_index = -1;
                break;
            } else {
                matching_index = k;
            }
        }
    }
    return matching_index;
}

void improve_ibps(vector<ibp_type>& ibps,SECTOR sector_fast, ibp_comparator &ibpcompare) {
    for (auto& ibp: ibps) {
        sort(ibp.begin(),ibp.end(), [&sector_fast](const pair<vector<COEFF>, point_fast > &a,
                                                    const pair<vector<COEFF>, point_fast > &b) -> bool {
            return vector_smaller_in_sector(b.second, a.second, sector_fast);
        });
    }
    // sorting each ibp, first are biggest

    if (common::disable_presolve) return;
    
    ibpcompare.s = sector_fast;
    sort(ibps.begin(), ibps.end(), ibpcompare);
    point_fast zero_point{};

    for (unsigned int i = 0; i!=ibps.size(); ++i) {
        //forward pass
        for (unsigned int j = i + 1; j!=ibps.size(); ++j) {
            if (ibps[i][0].second != ibps[j][0].second) {
                break;
            }
            int index = matching_index(ibps[i][0].first,ibps[j][0].first);
            if (index >= 0) {
                // only one index coeff differs
                // need to add i to j
                const COEFF & mul_i = ibps[j][0].first[index];
                const COEFF & mul_j = ibps[i][0].first[index];
                ibp_type res;
                add_ibps(mul_i, mul_j, ibps[i], ibps[j], sector_fast, res);
                // now res is the new relation that should replace j
                if (res.empty()) {
                    ibps.erase(ibps.begin()+j);
                } else if (res[0].second == zero_point) {
                    // we should not result in ibps with zero top shift since in can have cancelling coefficients!
                    ++j;
                } else {
                    ibps[j] = res;
                    sort(ibps.begin() + j, ibps.end(), ibpcompare);
                }
                --j;
            }
        }
    }
}

string now(bool use_date=false) {
    time_t timep;
    time (&timep);
    char tmp[64];
    if(use_date) strftime(tmp, sizeof(tmp), "%Y-%m-%d %H:%M:%S",localtime(&timep) );
    else strftime(tmp, sizeof(tmp), "%H:%M:%S",localtime(&timep) );
    return tmp;
}

bool equations_less(const equation &lhs, const equation &rhs);
/* main worker in a sector
* tries different methods
* such as searching for an sbasis or lbases
* if nothing uses Laporta
*/
void forward_stage(sector_count_t ssector_number) {
    if(common::run_mode>1 && read_status(ssector_number)!=0) return;
    if(common::run_mode==3) { common::sector_tasks.push_back(ssector_number); return; }
    
    atomic<uint64_t> virts_number{1};
    ibp_comparator ibpcompare;
    open_database(ssector_number);
    auto & dbn = DBM[ssector_number];
    
    set<point> set_needed_lower;
    set<sector_count_t> set_needed_lower_sector;

    vector<t_index> ssector = common::ssectors[ssector_number];
    unsigned short sector_level = static_cast<unsigned short>(std::count_if(ssector.begin(), ssector.end(),
            [&](const t_index &elem) {
                return elem == 1;
            }));

    set<point> needed_set;
    auto const & pm = dbn.pmap;
    for(auto const & kv : pm) {
        const point & test = kv.first;
        if ((test.virt()) || (test.s_number() != ssector_number)) continue;
        if (kv.second.first > 2 * sector_level) needed_set.insert(test);
    }

    set<point> *needed = &needed_set;
    point Corner = point_reference(corner(ssector));
    unsigned int n = ssector.size();
    bool first_pass = true;

    while (!needed->empty()) {
        set<point, indirect_more> ivpl;
        set<point>::iterator ivpl_counter;
        {
            for (const auto &read : *needed) {
                ivpl.insert(read);
            }
            bool done = true;
            for (ivpl_counter = ivpl.begin(); ivpl_counter != ivpl.end(); ++ivpl_counter) {
                point p = *ivpl_counter;
                vector<point> monoms = p_get_monoms(p);
                if (!monoms.empty()) {
                    for (const auto &monom : monoms) {
                        if (monom.s_number() == ssector_number) {
                            ivpl.insert(ivpl_counter, monom);
                        }
                    }
                } else {
                    done = false;
                    break;
                }
            }
            if (done) {
                ivpl.clear();
                for (const auto &read : *needed) {
                    ivpl.insert(read);
                }
                for (ivpl_counter = ivpl.begin(); ivpl_counter != ivpl.end(); ++ivpl_counter) {
                    point p = *ivpl_counter;
                    vector<point> monoms = p_get_monoms(p);
                    if (!monoms.empty()) {
                        for (const auto &monom : monoms) {
                            if (monom.s_number() == ssector_number) {
                                ivpl.insert(ivpl_counter, monom);
                            }
                            if (p.s_number() != monom.s_number()) {
                                set_needed_lower.insert(monom);
                                set_needed_lower_sector.insert(monom.s_number());
                            }
                        }
                    }
                }
                clear_sector(ssector_number, &ivpl);
                finish_sector(set_needed_lower, ssector_number);
                return;
            }
        }
        set<pair<unsigned int, unsigned int> > input_levels;
        set<pair<unsigned int, unsigned int> > real_input_levels;
        for (const auto &read : *needed) {
            const vector<t_index> v = read.get_vector();
            auto l = level(v);
            if (l.first == 0) {
                l.first = 1;
            }
            if (l.second == 0) {
                l.second = 1;
            }
            real_input_levels.insert(l);
            l = level(v);
            // using needed_level right here
            if (first_pass) {
                l.first = l.first + 1;
                if (l.second == 0) {
                    l.second = 1;
                }
            } else {
                l.first = l.first + 1;
                l.second = l.second + 1;
            }
            input_levels.insert(l);
            ivpl.insert(read);
        }

        // now checking for a new LBASIS
        auto litr = common::lbases.find(ssector_number);

        if (litr != common::lbases.end()) {
            //if (!common::silent) cout << "L-basis found." << endl;
            vector<pair<vector<pair<vector<t_index>, pair<short, bool> > >, vector<pair<string, vector<pair<vector<t_index>, short> > > > > > lbasis = litr->second;
            for (ivpl_counter = ivpl.begin(); ivpl_counter != ivpl.end(); ++ivpl_counter) {
                point p = *ivpl_counter;
                if (!p_is_empty(p)) {
                    continue;
                }
                vector<t_index> y = p.get_vector();
                bool made_table = false;
                for (const auto &elem : lbasis) {
                    vector<pair<vector<t_index>, pair<short, bool> > > conditions = elem.first;
                    bool satisfies = true;
                    for (const auto &condition : conditions) {
                        int sum = condition.second.first;
                        for (unsigned int i = 0; i != n; ++i) {
                            sum += (y[i] * (condition.first[i]));
                        }
                        if ((sum == 0) && condition.second.second == false) satisfies = false;
                        if ((sum != 0) && condition.second.second == true) satisfies = false;
                        if (!satisfies) {
                            break;
                        }
                    }
                    if (satisfies) { //rule satisfies all conditions
                        vector<pair<point, COEFF> > mon;
                        vector<pair<string, vector<pair<vector<t_index>, short> > > > terms = elem.second;
                        mon.reserve(terms.size());
                        map<vector<t_index>, string> local_coeffs;
                        for (const auto &term : terms) {
                            vector<t_index> new_v;
                            for (unsigned int i = 0; i != n; ++i) {
                                t_index sum = term.second[i].second;
                                for (unsigned int j = 0; j != n; ++j) {
                                    sum += (y[j] * term.second[i].first[j]);
                                }
                                new_v.push_back(sum);
                            }
                            string coeff = term.first;
                            subs2(coeff, y);
                            if (coeff != "0") {
                                auto local_coeffs_itr = local_coeffs.find(new_v);
                                if (local_coeffs_itr == local_coeffs.end()) {
                                    local_coeffs.emplace_hint(local_coeffs_itr, new_v, "(" + coeff + ")");
                                } else {
                                    local_coeffs_itr->second = local_coeffs_itr->second + "+" + "(" + coeff + ")";
                                }
                            }
                        }  // substituted all terms, joining equal
                        for (auto &&local_coeff : local_coeffs) {
                            string &ss = local_coeff.second;
                            ss = calc(ss);   // we are not using thread number here
                            #if defined(PRIME)
                            unsigned long long n = mod(ss);
                            char res[32];
                            sprintf(res, "%llu", n);
                            ss = string(res);
                            #endif
                            if (ss != "0") {
                                point new_p = point_reference(local_coeff.first);
                                if (new_p.is_zero()) {
                                    continue;
                                }
                                if (new_p.s_number() == 1) {
                                    cout << "Incorrect l-basis-rule for " << p << endl;
                                    cout << "Sending to " << new_p << endl;
                                    abort();
                                }
                                if (new_p.s_number() != p.s_number()) {
                                    if (new_p.level() >= p.level()) {
                                        cout << "Incorrect l-basis-rule for " << p << endl;
                                        cout << "Sending to " << new_p << endl;
                                        abort();
                                    }
                                }
                                mon.emplace_back(new_p, ss);
                            }
                        }

                        mon.emplace_back(p, CO_1m);
                        sort(mon.begin(), mon.end(), pair_point_COEFF_smaller);

                        if ((mon.empty()) || mon[mon.size() - 1].first != p) {
                            cout << "LBasis ordering error: " << p << " -> " << endl;
                            for (const auto &it : mon) {
                                cout << it.first << ", " << endl;
                            }
                            abort();
                        }
                        made_table = true;
                        for (unsigned int j = 0; j != mon.size(); ++j) {
                            if (p.s_number() != mon[j].first.s_number()) {
                                set_needed_lower.insert(mon[j].first);
                                set_needed_lower_sector.insert(mon[j].first.s_number());
                            }
                            if ((mon[j].first != p) && (mon[j].first.s_number() == ssector_number)) {
                                ivpl.insert(ivpl_counter, mon[j].first);
                            }
                        }
                        p_set(p, std::move(mon), 2 * sector_level);
                        break;
                    }  // if satisfies
                } // rule cycle
                if (!made_table) {
                    make_master(p);
                }
            }
            clear_sector(ssector_number, nullptr);
            finish_sector(set_needed_lower, ssector_number);
            return;
        }

        // should not we change ivpl to point_fast ???
        auto ditr = point::dbases.find(ssector_number);
        if (ditr != point::dbases.end()) {
            pair<vector<t_index>, vector<pair<vector<pair<COEFF, point_fast> >, t_index> > > dbasis = ditr->second;
            for (ivpl_counter = ivpl.begin(); ivpl_counter != ivpl.end(); ++ivpl_counter) {
                point p = *ivpl_counter;
                if (!p_is_empty(p)) continue;
                point_fast y = point_fast(p);
                list<vector<pair<COEFF, point_fast> > > product;
                point_fast new_y;
                for (unsigned int j = 0; j != n; j++) {
                    if (dbasis.first[j] == 0) {
                        new_y.buf[j] = 0;
                    } else {
                        new_y.buf[j] = y.buf[dbasis.first[j] - 1];
                    }
                }
                vector<pair<COEFF, point_fast> > new_term;
                new_term.emplace_back(CO_1, new_y); // we start building a product of sums. first term is without sum
                product.push_back(new_term);

                for (const auto &pitr : dbasis.second) {
                    int power = -y.buf[pitr.second - 1];
                    if (power < 0) {
                        cout << "Wrong delayed rule: " << p << endl;
                        abort();
                    }
                    for (int j = 0; j != power; j++) {
                        product.push_back(pitr.first);  // product is created from lbases. coefficients are strings
                    }
                }

                map<point_fast, COEFF> local_coeffs = calculate_product(product,
                                                                        n); // this is a shared function in external and internal symmetry generation

                // let's add -y right here
                auto itr = local_coeffs.find(y);
                if (itr == local_coeffs.end()) {
                    local_coeffs.insert(make_pair(y, CO_1m));
                } else {
                    _sub_(itr->second, itr->second, CO_1);
                }

                vector<pair<point, COEFF> > mon;
                mon.reserve(local_coeffs.size());

                for (const auto &local_coeff : local_coeffs) {
                    if (!is_zero(local_coeff.second))
                    {
                        point new_p = point_reference_fast(local_coeff.first);
                        if (new_p.is_zero()) {
                            continue;
                        }
                        if (new_p.s_number() == 1) {
                            cout << "Incorrect l-rule for " << p << endl;
                            cout << "Sending to " << new_p << endl;
                            abort();
                        }
                        if ((new_p.s_number() == p.s_number()) && (new_p != p)) {
                            cout << "Incorrect l-rule for " << p << endl;
                            cout << "Sending to " << new_p << endl;
                            abort();
                        }
                        if (new_p.s_number() != p.s_number()) {
                            if (new_p.level() > p.level()) {
                                cout << "Incorrect l-rule for " << p << endl;
                                cout << "Sending to " << new_p << endl;
                                abort();
                            }
                            if (new_p.level() == p.level()) {
                                if (common::lsectors.find(sector_fast(new_p.get_vector())) == common::lsectors.end()) {
                                    cout << "Incorrect l-rule for " << p << endl;
                                    cout << "Sending to " << new_p << endl;
                                    abort();
                                }
                            }
                        }
                        mon.emplace_back(new_p, local_coeff.second);
                    }
                }

                // product ready, now sorting, joining, evaluating
                sort(mon.begin(), mon.end(), sort_pair_point_coeff_by_point);
                mon = group_equal_in_sorted(mon);
                if ((mon.empty()) || mon[mon.size() - 1].first != p) {
                    cout << "LSymmetry ordering error: " << p << " -> " << endl;
                    for (const auto &it : mon) {
                        cout << it.first << ", " << endl;
                    }
                    abort();
                }

                for (unsigned int j = 0; j != mon.size(); ++j) {
                    if (p.s_number() != mon[j].first.s_number()) {
                        set_needed_lower.insert(mon[j].first);
                        set_needed_lower_sector.insert(mon[j].first.s_number());
                    }
                    if ((mon[j].first != p) && (mon[j].first.s_number() == ssector_number)) {
                        ivpl.insert(ivpl_counter, mon[j].first);
                    }
                }
                p_set(p, std::move(mon), 2 * sector_level + 1);
            } // writing symmetry done
            clear_sector(ssector_number, nullptr);
            finish_sector(set_needed_lower, ssector_number);
            return;
        }

        // ok, back to Laporta
        set<pair<unsigned int, unsigned int>, level_smaller> levels;
        auto const & um = dbn.umap;
        for (const auto &input_level : input_levels) {
            auto here_levels = under_levels(input_level.first, input_level.second);
            for (const auto &current_level : here_levels) {
                if (um.find(current_level)==um.end()) levels.insert(current_level);
            }
        }

        if (levels.empty()) {
            finish_sector(set_needed_lower, ssector_number);
            return;
        }

        point_fast p_fast(Corner);
        SECTOR sector_fast = p_fast.sector_fast();

        auto ibps = point::ibps;
        improve_ibps(ibps,sector_fast,ibpcompare);

        auto itr = levels.begin();

        for (unsigned int current_sum = 2; (itr != levels.end()); ++current_sum) {
            set<pair<unsigned int, unsigned int> > current_levels;
            while ((itr != levels.end()) && (((*itr).first) + ((*itr).second) <= current_sum)) {
                current_levels.insert(*itr);
                ++itr;
            }
            if (current_levels.empty()) {
                continue;
            }

            int symmetries = 0;
            if (common::pos_pref) {
                for (const auto &current_level : current_levels) {
                    if (common::pos_pref > 0) {
                        if ((current_level.first <= static_cast<unsigned int>(abs(common::pos_pref))) && (current_level.second == 1)) {
                            symmetries += write_symmetries(Corner, current_level.first, current_level.second, sector_level);
                        }
                    } else {
                        if ((current_level.first == 1) && (current_level.second <= static_cast<unsigned int>(abs(common::pos_pref)))) {
                            symmetries += write_symmetries(Corner, current_level.first, current_level.second, sector_level);
                        }
                    }
                }
            }
            
            vector<pair<unsigned int, unsigned int> > level_tasks;
            level_tasks.reserve(current_levels.size());
            for (auto level_itr = current_levels.rbegin(); level_itr != current_levels.rend(); ++level_itr) {
                level_tasks.push_back(*level_itr);
            }
            
{ // ======BEGIN reduce_in_level BEGIN======
    std::shared_mutex db_rw_mutex;
    size_t level_tasks_n = level_tasks.size();
    // seems only a few level_tasks is slow
    for(int level_tasks_i=0; level_tasks_i<level_tasks_n; level_tasks_i++) {
        sector_count_t ssector_number = Corner.s_number();

        vector<t_index> v = Corner.get_vector();
        point_fast p_fast(Corner);
        SECTOR sector_fast = p_fast.sector_fast();

        set<pair<unsigned int, unsigned int> > current_levels;
        current_levels.insert(level_tasks[level_tasks_i]); // we take the first task from the queue

        fstream out;
        fstream in;
        char readbuf[128];
        bool hint_exists = false;
        char buf[256];
        // we are expecting length 1 in current_levels now, keeping multiple for easy checks
        if (common::hint) {
            snprintf(buf, 256, "%s/%d-{%u,%u}.m", common::hint_path.c_str(), int(ssector_number),
                    current_levels.begin()->first, current_levels.begin()->second);
            in.open(buf, fstream::in);
            if (in.good()) {
                hint_exists = true;
                in.getline(readbuf, 64);
                if (strcmp(readbuf, "{") != 0) {
                    cout << "incorrect hint" << endl;
                    abort();
                }//}
            } else {
                in.close();
                out.open(buf, fstream::out);
                out << "{" << endl;//}
            }
        }

        vector<pair<point, pair<point_fast, unsigned short> > > ibps_vector; // the first point is just for sorting by sort_ibps; when loaded from file it can be filled with empty
        unsigned int eqs_number;

        if (common::hint && hint_exists) {
            // using the hint file to create lists of ibps
            eqs_number = 0;
            ibps_vector.reserve(std::count(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>(), '\n'));
            in.close();
            in.open(buf, fstream::in);
            in.getline(readbuf, 64);
            while (true) {
                in.getline(readbuf, 64);
                point_fast p;
                const char *pos = readbuf;
                if (*pos != '{') {
                    cout << "Wrong line start in hint;";
                    abort();
                } //}
                ++pos;
                if (*pos != '{') {
                    cout << "Wrong line start in hint;";
                    abort();
                } //}
                ++pos;
                int number;
                int i = 0;
                while (true) {
                    int move = s2i(pos, number);
                    p.buf[i] = number;
                    pos += move; //{
                    if (*pos == '}') {
                        ++pos;
                        break;
                    }
                    ++pos; // passing the comma
                    ++i;
                }
                if (*pos != ',') {
                    cout << "Wrong line middle in hint;";
                    abort();
                }
                ++pos;
                int move = s2i(pos, number);
                pos += move;
                if (*pos != '}') {
                    cout << "Wrong line end in hint;";
                    abort();
                }
                ++pos;
                ibps_vector.emplace_back(make_pair(point(), make_pair(p, number)));
                ++eqs_number;
                if (*pos == '}') break; // final closing bracket
            }
        } else {
            // preparing ibp lists
            vector<point_fast> IBPdegree;
            vector<point_fast> IBPdegreeFull;
            for (const auto &ibp : ibps) {
                // we fill IBPdegree, sscc and mm based of ibps
                point_fast p;
                point_fast ah = ibp[0].second;
                SECTOR bit = 1u << (common::dimension - 1);
                for (unsigned int i = 0; i < common::dimension; ++i, bit >>= 1)
                    if (sector_fast & bit) {
                        p.buf[i] = max(ah.buf[i], char(0));
                    } else {
                        p.buf[i] = max(-ah.buf[i], 0);
                    }
                IBPdegreeFull.emplace_back(ah);
                IBPdegree.emplace_back(p);
            }
            eqs_number = sort_ibps(Corner, current_levels, ibps_vector, IBPdegree, IBPdegreeFull);
        }

        unsigned int used_number = 0;
        vector<equation> eqs;

        if (hint_exists) eqs.resize(1);
        
        struct worker_item_type {
            int i1;
            int i2;
            const COEFF & o;
            uint64_t i3;
            worker_item_type(int _i1, int _i2, const COEFF & _o) : i1(_i1), i2(_i2), o(_o) { }
            worker_item_type(int _i1, int _i2, uint64_t _i3) : i1(_i1), i2(_i2), i3(_i3), o(CO_0) { }
        };
        
        vector<list<pc_pair_ptr_lst>::iterator> worker_itr;
        list<worker_item_type> worker_items;
        map<int,int> item_left;
        int worker_left = 0;
        mutex worker_mutex;
        condition_variable worker_cond;
        condition_variable worker_done_cond;
        bool worker_stop = false;
        int total_threads = common::lt1;
        if(total_threads<1) total_threads = 1;
        thread worker[total_threads];
        
        std::function<void(void)> worker_main;
        std::function<void(list<pc_pair_ptr_lst> &to_substitute)> worker_to_substitue;
        map<int, list<pc_pair_ptr_lst>> cache_mul_eqs;
        map<int,bool> item_in_use;
        
if(common::ltm) { // using cached equation

        worker_main = [&]() {
            int i1, i2;
            uint64_t i3;
            while(true) {
                auto itr = worker_items.begin();
                unique_lock<mutex> guard(worker_mutex);
                worker_cond.wait(guard, [&]() {
                    if(worker_stop) return true;
                    if(!worker_left) return false;
                    itr = worker_items.begin();
                    while(itr!=worker_items.end()) {
                        auto & ii1 = itr->i1;
                        auto & ii2 = itr->i2;
                        if(ii1==ii2 && item_left[ii1]==1) break;
                        if(ii1<ii2 && item_left[ii1]==0) break;
                        itr++;
                    }
                    if(itr!=worker_items.end()) return true;
                    else return false;
                });
                
                if (worker_stop) return;

                i1 = itr->i1;
                i2 = itr->i2;
                i3 = itr->i3;
                const COEFF & o = itr->o;
                worker_items.erase(itr);
                guard.unlock();
                    
                auto itrFrom = worker_itr[i1];
                auto itrTo = worker_itr[i2];
                
                if(i1==i2) {
                    auto const & eqs = cache_mul_eqs[i2];
                    for(auto eq : eqs) add_to(*itrTo, eq, false);
                    pc_pair_ptr_lst::iterator new_start;
                    { // write locker block
                        std::lock_guard<std::shared_mutex> write_locker(db_rw_mutex);
                        new_start = split(*itrTo, ssector_number, i3);
                    }
                    while (itrTo->begin() != new_start) itrTo->pop_front();
                    {
                        lock_guard<mutex> guard(worker_mutex);
                        worker_left--;
                        item_left[i2]--; // now item_left[i2] should be 0
                        if(!worker_left) cache_mul_eqs.clear();
                    }
                    if(!worker_left) worker_done_cond.notify_one();
                    else worker_cond.notify_all();
                } else {
                    COEFF c;
                    _div_neg_(c, o, (*(itrFrom->back())).second);
                    pc_pair_ptr_lst eq;
                    for(auto item : *itrFrom) {
                        COEFF co;
                        _mul_(co, c, item->second);
                        eq.emplace_back(make_pc_ptr(item->first, std::move(co)));
                    }
                    {
                        lock_guard<mutex> guard(worker_mutex);
                        item_left[i2]--;
                        cache_mul_eqs[i2].emplace_back(std::move(eq));
                    }
                }
            }
        };
        
        worker_to_substitue = [&](list<pc_pair_ptr_lst> &to_substitute) {
            //if(common::run_sector && !common::silent) cout << to_substitute.size() << flush;
            bool use_parallel = false;
            int eval_n = 0, eval_p;
            {
                lock_guard<mutex> guard(worker_mutex);
                int worker_total = to_substitute.size();
                worker_itr.clear();
                worker_items.clear();
                item_left.clear();
                cache_mul_eqs.clear();
                worker_itr.reserve(worker_total);
                for(auto itr = to_substitute.begin(); itr != to_substitute.end(); ++itr) {
                    worker_itr.push_back(itr);
                }
                vector<pc_pair_ptr_lst::iterator> term_vec(worker_total);
                vector<bool> changed_vec(worker_total);
                for(int i=0; i<worker_total; i++) {
                    term_vec[i] = worker_itr[i]->begin();
                    item_left[i] = 0;
                    changed_vec[i] = false;
                    cache_mul_eqs[i];
                }
                eval_n = 0;
                for(int i=0; i<worker_total; i++) {
                    if (changed_vec[i]) {
                        item_left[i]++;
                        worker_left++;
                        worker_items.emplace_back(i,i,++virts_number);
                    }
                    auto itrFrom = worker_itr[i];
                    const point &p = (*(itrFrom->rbegin()))->first;
                    for (int j=i+1; j<worker_total; j++) {
                        auto itrTo = worker_itr[j];
                        auto itrTerm = term_vec[j];
                        for (; itrTerm != itrTo->end(); ++itrTerm) {
                            if (p == (*itrTerm)->first) {
                                worker_items.emplace_back(i,j,(*itrTerm)->second);
                                changed_vec[j] = true;
                                item_left[j]++;
                                eval_n++;
                                ++itrTerm; // we move it before, or it will be invalidated
                                break;
                            } else if (p < (*itrTerm)->first) {
                                break; // this relation does not go there
                            }
                        }
                        term_vec[j] = itrTerm;
                    }
                }
                if(total_threads>1 && worker_left>common::lmt1) {
                    use_parallel = true;
                    eval_p = worker_left;
                } else {
                    worker_left = 0;
                }
            }
            
            if(use_parallel) {
                //if(common::run_sector && !common::silent) cout << "|" << eval_n << flush;
                worker_cond.notify_all();
                {
                    unique_lock<mutex> guard(worker_mutex);
                    worker_done_cond.wait(guard,[&worker_left](){ return !worker_left; });
                }
                //if(common::run_sector && !common::silent) cout << "|" << eval_p << flush;
            } else {
                for(auto & itr : worker_items) {
                    int i1 = itr.i1;
                    int i2 = itr.i2;
                    int i3 = itr.i3;
                    COEFF c(itr.o);
                    auto itrFrom = worker_itr[i1];
                    auto itrTo = worker_itr[i2];
                    if(i1==i2) {
                        pc_pair_ptr_lst::iterator new_start;
                        { // write locker block
                            std::lock_guard<std::shared_mutex> write_locker(db_rw_mutex);
                            new_start = split(*itrTo, ssector_number, i3);
                        }
                        while (itrTo->begin() != new_start) itrTo->pop_front();
                    } else {
                        _div_neg_(c, c, (*(itrFrom->back())).second);
                        mul_add_to(*itrTo, *itrFrom, c, false);
                    }
                }
                //if(common::run_sector && !common::silent) cout << "|" << eval_n << flush;
            }
        };
    
} else {
        
        worker_main = [&]() {
            while(true) {
                list<worker_item_type> wi_items;
                {
                    unique_lock<mutex> guard(worker_mutex);
                    worker_cond.wait(guard, [&]() {
                        if(worker_stop) return true;
                        if(!worker_left) return false;
                        bool first = true;
                        int last_i2;
                        auto itr = worker_items.begin();
                        while(itr!=worker_items.end()) {
                            auto & ii1 = itr->i1;
                            auto & ii2 = itr->i2;
                            if(first) {
                                if(ii1<ii2 && item_left[ii1]==0 && !item_in_use[ii2]) {
                                    wi_items.emplace_back(itr->i1, itr->i2, itr->o);
                                    itr = worker_items.erase(itr);
                                    first = false;
                                    last_i2 = ii2;
                                } else if(ii1==ii2 && item_left[ii1]==1 && !item_in_use[ii2]) {
                                    wi_items.emplace_back(itr->i1, itr->i2, itr->i3);
                                    itr = worker_items.erase(itr);
                                    last_i2 = ii2;
                                    break;
                                } else itr++;
                            } else if(ii2==last_i2 && ii1<ii2 && item_left[ii1]==0) {
                                wi_items.emplace_back(itr->i1, itr->i2, itr->o);
                                itr = worker_items.erase(itr);
                            } else if(ii2==last_i2 && ii1==ii2 && item_left[ii1]==wi_items.size()+1) {
                                wi_items.emplace_back(itr->i1, itr->i2, itr->i3);
                                itr = worker_items.erase(itr);
                                break;
                            } else itr++;
                        }
                        if(wi_items.size()>0) {
                            item_in_use[last_i2] = true;
                            return true;
                        } else return false;
                    });
                }
                if (worker_stop) return;
                
                bool hasii = false;
                int i2;
                for(auto item : wi_items) {
                    auto const & i1 = item.i1;
                    i2 = item.i2;
                    auto const & i3 = item.i3;
                    const COEFF & o = item.o;

                    auto itrFrom = worker_itr[i1];
                    auto itrTo = worker_itr[i2];

                    if(i1==i2) {
                        pc_pair_ptr_lst::iterator new_start;
                        { // write locker block
                            std::lock_guard<std::shared_mutex> write_locker(db_rw_mutex);
                            new_start = split(*itrTo, ssector_number, i3);
                        }
                        while (itrTo->begin() != new_start) itrTo->pop_front();
                        hasii = true;
                    } else {
                        COEFF c;
                        _div_neg_(c, o, (*(itrFrom->back())).second);
                        mul_add_to(*itrTo, *itrFrom, c, false);
                    }
                }
                {
                    lock_guard<mutex> guard(worker_mutex);
                    item_left[i2] -= wi_items.size();
                    item_in_use[i2] = false;
                    if(hasii) worker_left--;
                }
                if(hasii) {
                    if(!worker_left) worker_done_cond.notify_one();
                    else worker_cond.notify_all();
                }
            }
        };

        worker_to_substitue = [&](list<pc_pair_ptr_lst> &to_substitute) {
            if(total_threads>1) {
                //if(common::run_sector && !common::silent) cout << to_substitute.size() << flush;
                bool use_parallel = false;
                int eval_n = 0, eval_p;
                {
                    lock_guard<mutex> guard(worker_mutex);
                    int worker_total = to_substitute.size();
                    worker_itr.clear();
                    worker_items.clear();
                    item_in_use.clear();
                    item_left.clear();
                    worker_itr.reserve(worker_total);
                    for(auto itr = to_substitute.begin(); itr != to_substitute.end(); ++itr) {
                        worker_itr.push_back(itr);
                    }
                    vector<pc_pair_ptr_lst::iterator> term_vec(worker_total);
                    vector<bool> changed_vec(worker_total);
                    for(int i=0; i<worker_total; i++) {
                        term_vec[i] = worker_itr[i]->begin();
                        item_left[i] = 0;
                        changed_vec[i] = false;
                        item_in_use[i] = false;
                    }
                    eval_n = 0;
                    for(int i=0; i<worker_total; i++) {
                        if (changed_vec[i]) {
                            item_left[i]++;
                            worker_left++;
                            worker_items.emplace_back(i,i,++virts_number);
                        }
                        auto itrFrom = worker_itr[i];
                        const point &p = (*(itrFrom->rbegin()))->first;
                        for (int j=i+1; j<worker_total; j++) {
                            auto itrTo = worker_itr[j];
                            auto itrTerm = term_vec[j];
                            for (; itrTerm != itrTo->end(); ++itrTerm) {
                                if (p == (*itrTerm)->first) {
                                    worker_items.emplace_back(i,j,(*itrTerm)->second);
                                    changed_vec[j] = true;
                                    item_left[j]++;
                                    eval_n++;
                                    ++itrTerm; // we move it before, or it will be invalidated
                                    break;
                                } else if (p < (*itrTerm)->first) {
                                    break; // this relation does not go there
                                }
                            }
                            term_vec[j] = itrTerm;
                        }
                    }
                    if(worker_left>common::lmt1) {
                        use_parallel = true;
                        eval_p = worker_left;
                    } else {
                        worker_left = 0;
                    }
                }

                if(use_parallel) {
                    //if(common::run_sector && !common::silent) cout << "|" << eval_n << flush;
                    worker_cond.notify_all();
                    {
                        unique_lock<mutex> guard(worker_mutex);
                        worker_done_cond.wait(guard,[&worker_left](){ return !worker_left; });
                    }
                    //if(common::run_sector && !common::silent) cout << "|" << eval_p << flush;
                } else {
                    for(auto & itr : worker_items) {
                        int i1 = itr.i1;
                        int i2 = itr.i2;
                        int i3 = itr.i3;
                        COEFF c(itr.o);
                        auto itrFrom = worker_itr[i1];
                        auto itrTo = worker_itr[i2];
                        if(i1==i2) {
                            pc_pair_ptr_lst::iterator new_start;
                            { // write locker block
                                std::lock_guard<std::shared_mutex> write_locker(db_rw_mutex);
                                new_start = split(*itrTo, ssector_number, i3);
                            }
                            while (itrTo->begin() != new_start) itrTo->pop_front();
                        } else {
                            _div_neg_(c, c, (*(itrFrom->back())).second);
                            mul_add_to(*itrTo, *itrFrom, c, false);
                        }
                    }
                    //if(common::run_sector && !common::silent) cout << "|" << eval_n << flush;
                }
            } else { // original version
                int eval_n = 0;
                COEFF o;
                //if(common::run_sector && !common::silent) cout << to_substitute.size() << flush;
                for (auto itrTo = to_substitute.begin(); itrTo != to_substitute.end(); ++itrTo) {
                    bool changed = false;
                    auto startTerm = itrTo->begin();
                    for (auto itrFrom = to_substitute.cbegin(); itrFrom != itrTo; ++itrFrom) {
                        const point& p = (*(itrFrom->rbegin()))->first;
                        auto itrTerm = startTerm;
                        for (; itrTerm != itrTo->end(); ++itrTerm) {
                            if (p == (*itrTerm)->first) {
                                eval_n++;
                                _div_neg_(o, (*itrTerm)->second, (*(itrFrom->back())).second);
                                mul_add_to(*itrTo, *itrFrom, o, false);
                                changed = true;
                                ++itrTerm; // we move it before, or it will be invalidated
                                break;
                            } else if (p < (*itrTerm) -> first) {
                                break; // this relation does not go there
                            }
                        }
                        startTerm = itrTerm;
                    }
                    // top point cannot be changed, so it's safe to split
                    if (changed) {
                        pc_pair_ptr_lst::iterator new_start;
                        { // write locker block
                            std::lock_guard<std::shared_mutex> write_locker(db_rw_mutex);
                            new_start = split(*itrTo, ssector_number, ++virts_number);
                        }
                        while (itrTo->begin() != new_start) itrTo->pop_front();
                    }
                }
                //if(common::run_sector && !common::silent) cout << "|" << eval_n << flush;
            }
        };
        
}

    
        if(total_threads>1) for(int i=0; i<total_threads; ++i) worker[i] = thread(worker_main);
        
        auto worker_join = [&]() {
            if(total_threads<2) return;
            {
                lock_guard<mutex> guard(worker_mutex);
                worker_stop = true;
            }
            worker_cond.notify_all();
            for(int i=0; i<total_threads; ++i) worker[i].join();
        };
        
        /*
        if(common::run_sector && !common::silent) {
            cout << "\r";
            for(int j=0; j<100; j++) cout << " ";
            cout << "\r";
            cout << "Task: " << (level_tasks_i+1) << "/" << level_tasks_n << ", Equation: " << ibps_vector.size() << " @ " << now(true) << flush;
        }
        */
        int run_sector_i = 1;
        for (vector<pair<point, pair<point_fast, unsigned short> > >::const_iterator ibps_itr = ibps_vector.begin(); ibps_itr != ibps_vector.end(); ++ibps_itr) {
            if(common::run_sector && !common::silent) {
                cout << "\r";
                for(int j=0; j<100; j++) cout << " ";
                cout << "\r";
                cout << "Task: " << (level_tasks_i+1) << "/" << level_tasks_n << ", Equation: U" << used_number << "/E" << run_sector_i << "/T" << ibps_vector.size() << " @ " << now(true) << flush;
            }
            
            auto itr2 = ibps_itr;
            int k;
            int write;

            if (hint_exists) {
                // just generating the first equation from hint
                unsigned short i = (ibps_itr)->second.second;
                const point_fast &p = (ibps_itr)->second.first;
                eqs[0] = apply(ibps[i], p, sector_fast);
                write = 1; // there is just one equation
            } else { // moving until the highest member is changed
                int dist = 0;
                write = 0;
                while ((itr2 != ibps_vector.end()) && (itr2->first == ibps_itr->first)) {
                    ++itr2;
                    ++dist;
                }
                eqs.resize(dist); // with same highest member
                for (k = 0; k != dist; ++k) {
                    unsigned short i = (ibps_itr + k)->second.second;
                    const point_fast &p = (ibps_itr + k)->second.first;
                    // apply will call new
                    eqs[write] = apply(ibps[i], p, sector_fast);
                    // generating
                    if (eqs[write].length!=0) {
                        eqs[write].source = (ibps_itr + k)->second;
                        ++write;
                    }
                }
                sort(eqs.begin(), eqs.begin() + write, equations_less); // sorting those with same highest member
            }

            for (k = 0; k != write; ++k) {  //cycle of same starting point
                pc_pair_ptr_lst result;
                for (unsigned int i = 0; i != eqs[k].length; ++i) {
                    result.emplace_back(std::move(eqs[k].terms[i]));
                }

                list<pc_pair_ptr_lst> to_substitute;
                for (auto itr = result.rbegin(); itr != result.rend(); ++itr) {
                    const point &p = (*itr)->first;
                    if ((p.s_number() < ssector_number) ||
                        (p.virt()) ||
                        (p.s_number() == 1)) {
                        break; // no need to touch lower
                    } else {
                        pc_pair_ptr_lst terms2l;
                        { // read locker block
                            std::shared_lock<std::shared_mutex> read_locker(db_rw_mutex);
                            p_get(p, terms2l, ssector_number);
                        }
                        if (terms2l.empty()) {
                            continue;
                        } else {
                            to_substitute.emplace_front(terms2l);
                            COEFF o;
                            _div_neg_(o, (*itr)->second, terms2l.back()->second);
                            mul_add_to(result, terms2l, o, true); // last term still stays here
                            result.erase(next(itr--).base());
                            if (result.empty()) break;
                        }
                    }
                }
                
                // let's write the table for the current equation if needed
                //if(common::run_sector && !common::silent) cout << " " << flush;
                if (!result.empty()) {
                    const point &p = (*(result.rbegin()))->first;
                    if ((p.s_number() == ssector_number) && (!p.virt())) {
                        { // write locker block
                            std::lock_guard<std::shared_mutex> write_locker(db_rw_mutex);
                            split(result, ssector_number, ++virts_number);
                        }
                        if (common::hint && !hint_exists) {
                            unsigned short i = eqs[k].source.second;
                            point_fast &p = eqs[k].source.first;
                            if (used_number) out << "," << endl;
                            out << "{{";
                            for (unsigned j = 0; j + 1 != common::dimension; ++j) {
                                out << int(p.buf[j]) << ",";
                            }
                            out << int(p.buf[common::dimension - 1]) << "}" << "," << i << "}";
                        }
                        ++used_number;
                    }
                }
                
                // now substituting back
                worker_to_substitue(to_substitute);
            } // cycle of same starting point
            
            if (!hint_exists) {
                vector<equation> null_vec;
                std::swap(eqs, null_vec);
                ibps_itr = itr2;
                --ibps_itr;
                if(!common::silent) run_sector_i += k;
            }
        } // equation cycle
        if(common::run_sector && !common::silent) cout << endl;
        if (hint_exists) {
            vector<equation> null_vec;
            std::swap(eqs, null_vec);
        }
        if (common::hint) {
            if (hint_exists) in.close();
            else { out << "}" << endl; out.close(); }
        }
        worker_join();
    }
    //if(common::run_sector) cout << endl;
} // ======END reduce_in_level END ======

            auto & um = dbn.umap;
            for (const auto &current_level : current_levels) {
                um[current_level] = true;
                dbn.need_write = true;
                if (common::pos_pref > 1 && current_level.first == 1 && current_level.second == 1 && levels.size() > 1) {
                    // do not mark corner
                } else if (common::pos_pref > 1 && current_level.first + current_level.second == 3) {
                    mark_master_integrals(Corner, 1, 1);
                    mark_master_integrals(Corner, current_level.first, current_level.second);
                } else {
                    mark_master_integrals(Corner, current_level.first, current_level.second);
                }
            }

            ivpl.clear();
            for (const auto &read : *needed) {
                ivpl.insert(read);
            }
            bool done = true;
            for (ivpl_counter = ivpl.begin(); ivpl_counter != ivpl.end(); ++ivpl_counter) {
                point p = *ivpl_counter;
                vector<point> monoms = p_get_monoms(p);
                if (!monoms.empty()) {
                    for (const auto &monom : monoms) {
                        if (monom.s_number() == ssector_number) {
                            ivpl.insert(ivpl_counter, monom);
                        }
                    }
                } else {
                    done = false;
                    break;
                }
            }
            
            if (done) {

                ivpl.clear();
                for (const auto &read : *needed) {
                    ivpl.insert(read);
                }
                for (ivpl_counter = ivpl.begin(); ivpl_counter != ivpl.end(); ++ivpl_counter) {
                    point p = *ivpl_counter;
                    vector<point> monoms = p_get_monoms(p);
                    if (!monoms.empty()) {
                        for (const auto &monom : monoms) {
                            if (monom.s_number() == ssector_number) {
                                ivpl.insert(ivpl_counter, monom);
                            }
                            if (p.s_number() != monom.s_number()) {
                                set_needed_lower.insert(monom);
                                set_needed_lower_sector.insert(monom.s_number());
                            }
                        }
                    } else {
                        make_master(p);
                    }
                }
                clear_sector(ssector_number, &ivpl);
                finish_sector(set_needed_lower, ssector_number);
                return;
            }
        }

        ivpl.clear();
        for (const auto &read : *needed) {
            ivpl.insert(read);
        }

        for (ivpl_counter = ivpl.begin(); ivpl_counter != ivpl.end(); ++ivpl_counter) {
            point p = *ivpl_counter;
            vector<point> monoms = p_get_monoms(p);
            if (!monoms.empty()) {
                for (const auto &monom : monoms) {
                    if (monom.s_number() == ssector_number) {
                        ivpl.insert(ivpl_counter, monom);
                    }
                }
            } else {
                needed->insert(p);
            }
        }

        first_pass = false;
    }
    finish_sector(set_needed_lower, ssector_number);
}

/**
 * Function that compares two equations.
 * @param lhs first equation
 * @param rhs second equation
 * @return true if lhs is less than rhs, false otherwise.
 */
bool equations_less(const equation &lhs, const equation &rhs) {
    int i = lhs.length;
    int j = rhs.length;
    while ((i != 0) && (j != 0)) {
        const point &p1 = lhs.terms[i - 1]->first;
        const point &p2 = rhs.terms[j - 1]->first;
        if ((p1) < (p2)) return true;
        if ((p2) < (p1)) return false;
        --i;
        --j;
    }
    return (j!=0);
}

void clear_sector(sector_count_t sn, set<point, indirect_more> *ivpl) {
    // if we do not clear, there are too many extra entries in the database for relations we started solving
    // and they hurt at the substitution stage... unless we rewrite it with finding a list of needed for higher and then building the tree
    // and since we clear, the rules points are marked as absolutely needed
    
    if (common::keep_all) return;
    if (ivpl != nullptr) {
        auto & dbn = DBM[sn];
        dbn.need_write = true;
        auto & pm = dbn.pmap;
        for (auto it = pm.begin(); it != pm.end();) {
            const point & test = it->first;
            bool rm = (ivpl->find(test) == ivpl->end()) && (!common::only_masters || it->second.first != 127);
            if(rm) it = pm.erase(it);
            else ++it;
        }

    }
}

void perform_substitution(sector_count_t ssector_number) {
    if(common::run_mode>1 && read_status(ssector_number)==2) return;
    if(common::run_mode==3) { common::sector_tasks.push_back(-ssector_number); return; }
    
    open_database(ssector_number);
    auto & dbn = DBM[ssector_number];
    unsigned short v_level;
    if (common::keep_all) {
        v_level = (1 - in_lsectors(ssector_number)) + 2 * positive_index(common::ssectors[ssector_number]);
    } else {
        v_level = 0;
        // if we do not keep all entries, we just can pick all
    }
    
    set<point, indirect_more> v_s;
    auto & pm = dbn.pmap;
    for(auto const & kv : pm) {
        if(kv.second.first <= v_level) continue;
        const point & test = kv.first;
        if(test.s_number() == ssector_number) v_s.insert(test);
    }
    
    set<point, indirect_more> &needed = v_s;
    set<point, indirect_more>::reverse_iterator ritr;
    if (common::keep_all) {
        ritr = expressed_by(needed, ssector_number);  // if all entries are kept we have to pick the remaining
    } else {
        ritr = needed.rbegin();
    }
    pass_back(needed, ritr, ssector_number);

    if (common::keep_all) {
        v_level = 0;
    } else {
        v_level = (1 - in_lsectors(ssector_number)) + 2 * positive_index(common::ssectors[ssector_number]);
        // if we do not keep all entries, we just can pick all
    }
    
    for (auto it = pm.begin(); it != pm.end();) {
        const point & test = it->first;
        bool rm = it->second.first <= v_level;
        if(rm) it = pm.erase(it);
        else ++it;
    }
    if(dbn.status!=1) {
        cout << "finish_sector: dbn.status!=0" << endl;
        abort();
    }
    dbn.need_write = true;
    dbn.status = 2;
    close_database(ssector_number);
}

void run_child(int sector) {
    
    if (sector > 0) forward_stage(sector);
    else if (sector < 0) perform_substitution(-sector);
    else {
        cout << "Sector not specified" << endl;
        abort();
    }
    
    flint_cleanup();
    
}

void Evaluate() {
    int last_level = common::abs_max_level;

    // forward reduction level by level
    while (last_level >= common::abs_min_level) {
        int inlsectors = 1;
        while (inlsectors != -1) {
            bool inlsectorsbool = (inlsectors == 0);
            inlsectors--;

            // list of sectors of this level
            set<sector_count_t> sector_set_this_level;
            set<sector_count_t>::reverse_iterator sector_set_this_level_itr;
            for (sector_count_t test_sector = 2; test_sector <= common::abs_max_sector; ++test_sector) {
	      
                vector<t_index> ssector = common::ssectors[test_sector];

                unsigned short current_level = static_cast<unsigned short>(std::count_if(ssector.begin(), ssector.end(),
                        [&](const t_index &elem) {
                            return elem == 1;
                        }));

                if (!common::sector_numbers_fast[sector_fast(ssector)]) continue;

                if ((inlsectorsbool == in_lsectors(test_sector)) && (current_level == last_level) &&
                    database_exists(test_sector)) {
                    sector_set_this_level.insert(test_sector);
                }
            }
            
            if (!common::silent) cout << "> LEVEL " << last_level << "." << inlsectors + 1 << " @ " << now(true) << endl;
            vector<int> worker_tasks;
            for (sector_set_this_level_itr = sector_set_this_level.rbegin();
                sector_set_this_level_itr != sector_set_this_level.rend(); ++sector_set_this_level_itr) {
                worker_tasks.push_back(*sector_set_this_level_itr);  // we created the list of jobs
            }
        
            std::atomic_flag lock = ATOMIC_FLAG_INIT;
            std::atomic<int> n{0};
            auto worker_tasks_n = worker_tasks.size();
            if (!common::silent && worker_tasks_n) cout << "  \\--Start @ " << now(true) << endl;
            int nts = common::t1;
            if(nts>worker_tasks_n) nts = worker_tasks_n;
            #pragma omp parallel for schedule(dynamic,1) num_threads(nts) if(common::t1>1 && worker_tasks_n>1 && common::run_mode!=3)
            for(int ti=0; ti<worker_tasks_n; ++ti) {
                while (lock.test_and_set(std::memory_order_acquire)) ;
                if (!common::silent) cout << "\r                       \r  \\--Evaluating " << (++n) << "/" << worker_tasks_n << " @ " << now() << flush;
                lock.clear();
                run_child(worker_tasks[ti]);
            }
            if (!common::silent && worker_tasks_n) cout << endl << "  \\--Done  @ " << now(true) << endl;
            if(common::run_mode==3 && common::sector_tasks.size()>0) {
                for(auto item : common::sector_tasks) cout << item << endl;
                exit(0);
            }

            // here we will be writing the integral requests into lower sector databases
            map<sector_count_t, set<point> > needed_lower;

            for (sector_set_this_level_itr = sector_set_this_level.rbegin(); sector_set_this_level_itr != sector_set_this_level.rend(); ++sector_set_this_level_itr) {
                sector_count_t test_sector = (*sector_set_this_level_itr);

                open_database(test_sector, 1);
                auto & dbn = DBM[test_sector];
                int64_t new_size = dbn.lower_size;
                if (new_size > 0) {
                    auto const & point_vec = dbn.lower;
                    for (int i = 0; i != new_size; ++i) {
                        add_needed(needed_lower, point_vec[i]);
                    }
                }
                close_database(test_sector);
            }
            
            for (const auto &sector : needed_lower) {
                if(common::run_mode>1 && read_status(sector.first)!=0) continue;
                open_database(sector.first);
                for (const auto &pnt : sector.second) {
                    pc_pair_ptr_vec t;
                    p_set(pnt, std::move(t), 2 * last_level + inlsectors + 1);
                }
                close_database(sector.first);
            }

            if (!common::silent) {
                cout << "  LEVEL " << last_level << "." << inlsectors + 1 << " @ " << now(true) << endl;
            }
        }
        last_level--;
    }

    // backward substitutions level by level
    if (not common::only_masters) {
        if (!common::silent) cout << "----------------------------------------" << endl;
        while (last_level <= common::abs_max_level) {
            if (!common::silent) {
                cout << "< LEVEL " << last_level << " @ " << now(true) << endl;
            }
            int inlsectors = 0;
            while (inlsectors != 2) {
                bool inlsectorsbool = (inlsectors == 0);

                // list of sectors of this level
                set<sector_count_t> sector_set_this_level;
                set<sector_count_t>::iterator sector_set_this_level_itr;

                for (sector_count_t i = 2; i != common::abs_max_sector + 1; i++) {
                    if (!common::sector_numbers_fast[sector_fast(common::ssectors[i])]) continue;
                    if (database_exists(i) && (positive_index(common::ssectors[i]) == last_level) && (in_lsectors(i) == inlsectorsbool)) {
                        sector_set_this_level.insert(i);
                    }
                }

                inlsectors++;
                if (sector_set_this_level.empty()) {
                    continue;
                }

                for (sector_set_this_level_itr = sector_set_this_level.begin();
                        sector_set_this_level_itr != sector_set_this_level.end();
                        ++sector_set_this_level_itr) {
                    sector_count_t test_sector = *sector_set_this_level_itr;
                    if(common::run_mode>1 && read_status(test_sector)==2) continue;
                    open_database(test_sector);
                    auto & dbn = DBM[test_sector];
                    map<sector_count_t, set<point> > needed_lower;
                    
                    // here we read and fill needed_for for the current sector
                    int64_t new_size = dbn.lower_size;
                    if (new_size > 0) {
                        auto const & point_vec = dbn.lower;
                        for (int i = 0; i != new_size; ++i) {
                            //needed_for[test_sector].insert(point_vec[i].s_number());
                            add_needed(needed_lower, point_vec[i]);
                        }
                    }

                    // we are moving entries from lower database to higher
                    for (const auto &sector : needed_lower) {
                        open_database(sector.first, 0, &(sector.second));
                        auto const & pm = DBM[sector.first].pmap;
                        for (const auto &p : sector.second) {
                            const auto & itr = pm.find(p);
                            if (itr == pm.end()) {
                                cout << "Error on point transfer: " << p.s_number() << p << endl;
                                abort();
                            }
                            if(common::run_mode) dbn.pmap[p] = std::move(itr->second); // maybe some improvement
                            else dbn.pmap[p] = itr->second;
                        }
                        close_database(sector.first);
                    }
                    dbn.need_write = true;
                    close_database(test_sector);
                }

                vector<int> worker_tasks;
                for (sector_set_this_level_itr = sector_set_this_level.begin(); sector_set_this_level_itr != sector_set_this_level.end(); ++sector_set_this_level_itr) {
                    worker_tasks.push_back(-(*sector_set_this_level_itr));  // we created the list of jobs
                }
                auto worker_tasks_n = worker_tasks.size();
                std::atomic<int> n{0};
                std::atomic_flag lock = ATOMIC_FLAG_INIT;
                if (!common::silent && worker_tasks_n) cout << "  \\--Start @ " << now(true) << endl;
                int nts = common::t2;
                if(nts>worker_tasks_n) nts = worker_tasks_n;
                #pragma omp parallel for schedule(dynamic,1) num_threads(nts) if(common::t2>1 && worker_tasks_n>1 && common::run_mode!=3)
                for(int ti=0; ti<worker_tasks_n; ++ti) {
                    while (lock.test_and_set(std::memory_order_acquire)) ;
                    if (!common::silent) cout << "\r                       \r  \\--Evaluating " << (++n) << "/" << worker_tasks_n << " @ " << now() << flush;
                    lock.clear();
                    run_child(worker_tasks[ti]);
                }
                if (!common::silent && worker_tasks_n) cout << endl << "  \\--Done  @ " << now(true) << endl;
                if(common::run_mode==3 && common::sector_tasks.size()>0) {
                    for(auto item : common::sector_tasks) cout << item << endl;
                    exit(0);
                }
            }
            
            if (!common::silent) {
                cout << "  LEVEL " << last_level << " @ " << now(true) << endl;
            }
            last_level++;
        }

    }  // if making substitutions

}

set<point, indirect_more>::reverse_iterator expressed_by(set<point, indirect_more> &to_test, sector_count_t sector_number) {
    for (auto itr = to_test.begin(); itr != to_test.end(); ) {
        point p = *itr;
        if (p.s_number() != sector_number) {
            return set<point, indirect_more>::reverse_iterator(itr);
        }
        vector<point> monoms = p_get_monoms(p, sector_number);
        if (!monoms.empty()) {
            auto last = monoms.end();
            last--;
            for (auto mitr = monoms.begin(); mitr != last; ++mitr) {
                if (mitr->s_number() == sector_number) {
                    to_test.insert(itr, (*mitr));
                }
            }
            ++itr;
        } else {
            to_test.erase(itr++);
        }
    }
    return to_test.rbegin();
}

void mul_add_to(pc_pair_ptr_lst &terms1, pc_pair_ptr_lst::const_iterator termsB, pc_pair_ptr_lst::const_iterator termsE, const COEFF &coeff, bool skip_last) {
    auto eq2end = termsE;
    if (skip_last) { eq2end--; }
    auto itr1 = terms1.begin();
    auto itr2 = termsB;
    while (itr2 != eq2end) {
        if ((itr1 != terms1.end()) && ((*itr1)->first < (*itr2)->first)) { // first equation only. just adding to result
            ++itr1;
        } else if ((itr1 == terms1.end()) || ((*itr2)->first < (*itr1)->first)) { // second equation only. have to multiply
            if(!is_zero((*itr2)->second)) {
                COEFF o;
                _mul_(o, (*itr2)->second,  coeff);
                itr1 = terms1.emplace(itr1, make_pc_ptr((*itr2)->first, std::move(o)));
                ++itr1;
            }
            ++itr2;
        } else { // both equations. have to multiply and add
            COEFF o;
            _mul_(o, (*itr2)->second, coeff);
            _add_(o, o, (*itr1)->second);
            if(!o.is_zero()) {
                (*itr1) = make_pc_ptr((*itr1)->first, std::move(o));
                ++itr1;
            } else {
                itr1 = terms1.erase(itr1);
            }
            ++itr2;
        }
    }
}
void mul_add_to(pc_pair_ptr_lst &terms1, const pc_pair_ptr_lst &terms2, const COEFF &coeff, bool skip_last) {
    mul_add_to(terms1, terms2.begin(), terms2.end(), coeff, skip_last);
}

void add_to(pc_pair_ptr_lst &terms1, pc_pair_ptr_lst::const_iterator termsB, pc_pair_ptr_lst::const_iterator termsE, bool skip_last) {
    auto eq2end = termsE;
    if (skip_last) { eq2end--; }
    auto itr1 = terms1.begin();
    auto itr2 = termsB;
    while (itr2 != eq2end) {
        if ((itr1 != terms1.end()) && ((*itr1)->first < (*itr2)->first)) { // first equation only. just adding to result
            ++itr1;
        } else if ((itr1 == terms1.end()) || ((*itr2)->first < (*itr1)->first)) { // second equation only. have to multiply
            if(!is_zero((*itr2)->second)) {
                itr1 = terms1.emplace(itr1, *itr2);
                ++itr1;
            }
            ++itr2;
        } else { // both equations. have to multiply and add
            COEFF o;
            _add_(o, (*itr1)->second, (*itr2)->second);
            if(!o.is_zero()) {
                (*itr1) = make_pc_ptr((*itr1)->first, std::move(o));
                ++itr1;
            } else {
                itr1 = terms1.erase(itr1);
            }
            ++itr2;
        }
    }
}
void add_to(pc_pair_ptr_lst &terms1, const pc_pair_ptr_lst &terms2, bool skip_last) {
    add_to(terms1, terms2.begin(), terms2.end(), skip_last);
}

void apply_table(const pc_pair_ptr_vec &terms,
                 sector_count_t fixed_database_sector, unsigned short sector_level) {
    bool changed = false;
    auto end = terms.cend();
    --end;
    pc_pair_ptr_lst result;
    COEFF o;
    for (auto itr = terms.cbegin(); itr != end; ++itr) {
        const point &p = (*itr)->first;
        if (p.s_number() == 1) {
            result.emplace_back(*itr);
        } else {
            pc_pair_ptr_lst terms2;
            p_get(p, terms2, fixed_database_sector);
            if (terms2.empty()) {
                result.emplace_back(*itr);  // we just put the monomial at the end with no substitution
            } else {
                changed = true;
                _div_neg_(o, (*itr)->second, terms2.back()->second);
                mul_add_to(result, terms2, o, true); // we add the second expression with the last term killed
            }
        }
    }
    if (!changed) return;
    result.emplace_back(terms.back());
    p_set(result.back()->first, result, 2 * sector_level, fixed_database_sector);
}

void pass_back(const set<point, indirect_more> &cur_set, set<point, indirect_more>::const_reverse_iterator ritr, sector_count_t fixed_database_sector) {
// -- original version
//    for (auto itr = ritr; itr != cur_set.rend(); ++itr) {
//        point p = *itr;
//        pc_pair_ptr_vec terms;
//        p_get(p, terms, fixed_database_sector);
//        if (!terms.empty()) {
//            apply_table(terms, fixed_database_sector, 0);  // backward reduction with 0 as sector and 0 as thread_number
//        }
//    }

    struct worker_item_type {
        int i;
        set<int> iset;
        worker_item_type(int _i, set<int> && _iset) : i(_i), iset(std::move(_iset)) { }
    };
    
    vector<point> worker_point;
    list<worker_item_type> worker_items;
    int worker_left = 0;
    mutex worker_mutex;
    condition_variable worker_cond;
    int total_threads = common::lt2;
    if(total_threads<1) total_threads = 1;
    thread worker[total_threads];
    map<int,bool> item_ready;

    auto worker_main = [&]() {
        int i;
        while(true) {
            auto itr = worker_items.begin();
            unique_lock<mutex> guard(worker_mutex);
            worker_cond.wait(guard, [&]() {
                if(!worker_left) return true;
                itr = worker_items.begin();
                while(itr!=worker_items.end()) {
                    auto & iset = itr->iset;
                    bool ok = true;
                    for(auto ii : iset) if(!item_ready[ii]) {
                        ok = false;
                        break;
                    }
                    if(ok) break;
                    itr++;
                }
                if(itr!=worker_items.end()) return true;
                else return false;
            });

            if (!worker_left) return;
            i = itr->i;
            worker_items.erase(itr);
            guard.unlock();

            const point & p = worker_point[i];
            pc_pair_ptr_vec terms;
            p_get(p, terms, fixed_database_sector);
            if (!terms.empty()) apply_table(terms, fixed_database_sector, 0);
            {
                lock_guard<mutex> guard(worker_mutex);
                worker_left--;
                item_ready[i] = true;
            }
            worker_cond.notify_all();
            if(!worker_left) return;
        }
    };

    auto worker_pass_back = [&](const set<point, indirect_more> &cur_set) {
        if(total_threads>1) {
            if(common::run_sector && !common::silent) cout << cur_set.size() << flush;
            {
                lock_guard<mutex> guard(worker_mutex);
                worker_left = 0;
                worker_point.clear();
                worker_items.clear();
                worker_point.reserve(cur_set.size());
                map<point, int> p2i;
                for (auto citr = ritr; citr != cur_set.rend(); ++citr) {
                    const point & p = *citr;
                    worker_point.emplace_back(p);
                    
                    auto dsector = (fixed_database_sector == 0) ? p.s_number() : fixed_database_sector;
                    auto const & pm = DBM[dsector].pmap;
                    auto itr = pm.find(p);
                    if(itr==pm.end()) continue;
                    
                    auto const & terms = itr->second.second;
                    if(terms.empty()) continue;
                    
                    set<int> iset;
                    bool changed = false;
                    auto end = terms.cend();
                    --end;
                    for (auto itr2 = terms.cbegin(); itr2 != end; ++itr2) {
                        const point & p2 = (*itr2)->first;
                        if (p2.s_number() == 1) continue;
                        else {
                            auto itr3 = p2i.find(p2);
                            if(itr3 != p2i.end()) iset.insert(itr3->second);
                            if (!changed && pm.find(p2) != pm.end()) changed = true;
                        }
                    }
                    p2i[p] = worker_left;
                    item_ready[worker_left] = !changed;
                    worker_items.emplace_back(worker_left, std::move(iset));
                    worker_left++;
                }
            }

            if(worker_left>common::lmt2) {
                for(int i=0; i<total_threads; ++i) worker[i] = thread(worker_main);
                for(int i=0; i<total_threads; ++i) worker[i].join();
                if(common::run_sector && !common::silent) cout << ".." << flush;
            } else {
                for(auto & item : worker_items) {
                    const point & p = worker_point[item.i];
                    pc_pair_ptr_vec terms;
                    p_get(p, terms, fixed_database_sector);
                    if (!terms.empty()) {
                        apply_table(terms, fixed_database_sector, 0);  // backward reduction with 0 as sector and 0 as thread_number
                    }
                }
                if(common::run_sector && !common::silent) cout << "." << flush;
            }
        } else { // original version
            if(common::run_sector && !common::silent) cout << cur_set.size() << flush;
            for (auto itr = ritr; itr != cur_set.rend(); ++itr) {
                point p = *itr;
                pc_pair_ptr_vec terms;
                p_get(p, terms, fixed_database_sector);
                if (!terms.empty()) {
                    apply_table(terms, fixed_database_sector, 0);  // backward reduction with 0 as sector and 0 as thread_number
                }
            }
            if(common::run_sector && !common::silent) cout << "." << flush;
        }
    };
    
    worker_pass_back(cur_set);
}

pc_pair_ptr_lst::iterator split(pc_pair_ptr_lst &terms, sector_count_t sector_number, uint64_t virts_temp) {
    if (terms.empty()) return terms.begin();
    auto itr = terms.begin();
    size_t size = 0;
    while (((*itr)->first.s_number() < sector_number) || (*itr)->first.virt()) {
        ++itr;
        ++size;
    }
    auto pidx = 2 * positive_index(common::ssectors[sector_number]);
    if ((itr != terms.begin()) && (itr != terms.end())) {
        //uint64_t virts_temp = ++virts_number; // atomic
        point p(common::ssectors[sector_number], virts_temp); // virtual
        pc_pair_ptr save = *itr;
        (*itr) = make_pc_ptr(p, CO_1m);
        
        ++itr;
        p_set(p, size + 1, terms.begin(), itr, pidx, 0);  // split won't happen in symmetries sector
        --itr;
        (*itr) = save;
        
        --itr;
        (*itr) = make_pc_ptr(p, CO_1);
        
        p_set(terms.back()->first, terms.size() - size + 1, itr, terms.end(), pidx, 0);
        return itr;
    } else p_set(terms.back()->first, terms, pidx, 0);  // split won't happen in symmetries sector
    return terms.begin();
}

bool sort_pair_point_coeff_by_point(const pair<point, COEFF> &lhs, const pair<point, COEFF> &rhs) {
    return lhs.first < rhs.first;
}

pc_pair_vec group_equal_in_sorted(const vector<pair<point, COEFF> > &mon) {
    pc_pair_vec terms;
    terms.reserve(mon.size());
    for (auto read = mon.begin(); read != mon.end(); ++read) {
        // there is already a check here for equal points due to symmetries
        COEFF o;
        _set_(o, read->second);
        auto read2 = read;
        read2++;
        while ((read2 != mon.end()) && (read2->first == read->first)) {
            _add_(o, o, read2->second);
            read2++;
        }
        read2--;
        read = read2;
        if (!o.is_zero()) {
            terms.emplace_back(read->first, std::move(o));
        }
    }
    return terms;
}
