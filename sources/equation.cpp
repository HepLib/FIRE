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

mutex equation::f_submit_mutex[MAX_THREADS];
mutex equation::f_receive_mutex[MAX_THREADS];
condition_variable equation::f_submit_cond[MAX_THREADS];
condition_variable equation::f_receive_cond[MAX_THREADS];
list<pair<int, pair<point, string> > > equation::f_jobs[MAX_THREADS];
bool equation::f_stop = false;
list<pair<point, string> > equation::f_result[MAX_THREADS];
thread equation::f_threads[MAX_THREADS];


// get monoms from a point (Feynman integral). database access used
vector<point> p_get_monoms(const point &p, unsigned short fixed_database_sector) {
    vector<point> result;
    unsigned short dsector = (fixed_database_sector == 0) ? p.s_number() : fixed_database_sector;

    size_t size;
    char *res = common::points[dsector]->get(p.ww, sizeof(point), &size);
    if (res == nullptr) {
        if (common::points[dsector]->error().code() != kyotocabinet::BasicDB::Error::Code::NOREC) {
            cout << p.s_number() << p << string(common::points[dsector]->error().message()) << endl;
            abort();
        }
        return result;
    }

    unsigned short len = reinterpret_cast<unsigned short *>(res)[0];
    result.reserve(len);
    for (unsigned int i = 0; i != len; ++i) {
        result.emplace_back((reinterpret_cast<point *>(res + 3))[i]);
    }

    delete[] res;
    if ((len != 0) && (p != result.back())) {
        cout << "get_monoms error" << endl;
        cout << p << endl;
        cout << result.size() << endl;
        for (const auto & pnt : result) {
            cout << pnt << endl;
        }
        abort();
    }
    return result;
}

bool p_is_empty(const point &p, unsigned short fixed_database_sector) {
    unsigned short dsector = (fixed_database_sector == 0) ? p.s_number() : fixed_database_sector;
    size_t size;
    char *res = common::points[dsector]->get(p.ww, sizeof(point), &size);
    if (res == nullptr) {
        if (common::points[dsector]->error().code() != kyotocabinet::BasicDB::Error::Code::NOREC) {
            cout << p.s_number() << p << string(common::points[dsector]->error().message()) << endl;
            abort();
        }
        return true;
    }
    unsigned short len = reinterpret_cast<unsigned short *>(res)[0];
    delete[] res;
    return len == 0;
}

void p_get(const point &p, vector<pair<point, COEFF> > &terms, unsigned short fixed_database_sector) {
    unsigned short dsector = (fixed_database_sector == 0) ? p.s_number() : fixed_database_sector;
    size_t size;
    char *res = common::points[dsector]->get(p.ww, sizeof(point), &size);
    if (res == nullptr) {
        if (common::points[dsector]->error().code() != kyotocabinet::BasicDB::Error::Code::NOREC) {
            cout << p.s_number() << p << string(common::points[dsector]->error().message()) << endl;
            abort();
        }
        terms.clear();
        return;
    }

    unsigned short len = reinterpret_cast<unsigned short *>(res)[0];
    if (len == 0) {
        terms.clear();
        delete[] res;
        return;
    } // empty point

    terms.reserve(len);  // this is the only difference! I could not make this with templates in c++11 (only 17 I think)

    p_get_internal(p, terms, dsector, len, res);
}

void p_get(const point &p, list<pair<point, COEFF>, ALLOCATOR2> &terms, unsigned short fixed_database_sector) {
    unsigned short dsector = (fixed_database_sector == 0) ? p.s_number() : fixed_database_sector;
    size_t size;
    char *res = common::points[dsector]->get(p.ww, sizeof(point), &size);
    if (res == nullptr) {
        if (common::points[dsector]->error().code() != kyotocabinet::BasicDB::Error::Code::NOREC) {
            cout << p.s_number() << p << string(common::points[dsector]->error().message()) << endl;
            abort();
        }
        terms.clear();
        return;
    }

    unsigned short len = reinterpret_cast<unsigned short *>(res)[0];
    if (len == 0) {
        terms.clear();
        delete[] res;
        return;
    } // empty point

    p_get_internal(p, terms, dsector, len, res);
}

#ifdef MPQBin
size_t mpq_size(const mpq_class & q) {
    if(sizeof(char)!=1) throw std::runtime_error("char size is NOT 1!");
    size_t s1 = (mpz_sizeinbase(q.get_num().get_mpz_t(), 2) + CHAR_BIT-1) / CHAR_BIT;
    size_t s2 = (mpz_sizeinbase(q.get_den().get_mpz_t(), 2) + CHAR_BIT-1) / CHAR_BIT;
    return sizeof(bool)+2*sizeof(size_t)+s1+s2;
}
void mpq_write(const mpq_class & q, char* &pos) {
    if(sizeof(char)!=1) throw std::runtime_error("char size is NOT 1!");
    bool neg = q.get_num() < 0;
    memcpy(pos, &neg, sizeof(neg));
    pos += sizeof(bool);
    size_t c;
    
    c = (mpz_sizeinbase(q.get_num().get_mpz_t(), 2) + CHAR_BIT-1) / CHAR_BIT;
    memcpy(pos, &c, sizeof(size_t));
    pos += sizeof(size_t);
    mpz_export(pos, &c, 1, 1, 0, 0, q.get_num().get_mpz_t());
    pos += c;
    
    c = (mpz_sizeinbase(q.get_den().get_mpz_t(), 2) + CHAR_BIT-1) / CHAR_BIT;
    memcpy(pos, &c, sizeof(size_t));
    pos += sizeof(size_t);
    mpz_export(pos, &c, 1, 1, 0, 0, q.get_den().get_mpz_t());
    pos += c;
}
mpq_class mpq_read(char* &pos) {
    if(sizeof(char)!=1) throw std::runtime_error("char size is NOT 1!");
    bool neg;
    mpz_class num, den;
    memcpy(&neg, pos, sizeof(neg));
    pos += sizeof(bool);
    size_t c;
    memcpy(&c, pos, sizeof(size_t));
    pos += sizeof(size_t);
    mpz_import(num.get_mpz_t(), c, 1, 1, 0, 0, pos);
    pos += c;
    memcpy(&c, pos, sizeof(size_t));
    pos += sizeof(size_t);
    mpz_import(den.get_mpz_t(), c, 1, 1, 0, 0, pos);
    pos += c;
    mpq_class res = mpq_class(num)/mpq_class(den);
    if(neg) res = -res;
    return res;
}
#endif

#if defined(FlintC3)
inline slong str_len3(const fmpz_poly_t p) {
    static auto l62 = log((double)62);
    slong len = p->length;
    if(len==0) return 1;
    fmpz* poly = p->coeffs;
    slong bound;
    bound = (slong) (ceil(log((double) (len + 1))/l62));
    for (slong i = 0; i < len; i++) bound += fmpz_sizeinbase(poly + i, 62) + 1;
    bound += len + 2;
    return bound;
}
slong str_len3(const fmpz_poly_q_t p) {
    slong len = str_len3(fmpz_poly_q_numref(p));
    if(!fmpz_poly_is_one(fmpz_poly_q_denref(p))) len += 1+str_len3(fmpz_poly_q_denref(p));
    return len;
}

inline void _fmpz_poly_get_str_(const fmpz_poly_t p, char* &str) {
    slong len = p->length;
    if (len == 0) {
        str[0] = '0';
        str += 1;
    }
    fmpz * poly = p->coeffs;
    str += flint_sprintf(str, "%wd ", len);
    do {
        *str++ = ' ';
        if (!COEFF_IS_MPZ(*poly)) {
            mpz_t z;
            flint_mpz_init_set_si(z, *poly);
            mpz_get_str(str, 62, z);
            mpz_clear(z);
        } else mpz_get_str(str, 62, COEFF_TO_PTR(*poly));
        while(*str!='\0') str++;
    } while (poly++, --len);
}

void write_str3(const fmpz_poly_q_struct* p, char* &pos) {
    _fmpz_poly_get_str_(fmpz_poly_q_numref(p), pos);
    if(!fmpz_poly_is_one(fmpz_poly_q_denref(p))) {
        *pos = '/';
        pos++;
        _fmpz_poly_get_str_(fmpz_poly_q_denref(p), pos);
    }
}

inline int _fmpz_poly_set_str_(fmpz * poly, char *str) {
    slong i, len;
    if (!isdigit((unsigned char) str[0])) return -1;
    len = atol(str);
    if (len < 0) return -1;
    if (len == 0) return 0;
    while (*str!=' ') str++; 
    char *v;
    for (i = 0; i < len; i++) {
        while (*str==' ') str++;
        v = str;
        while(*str != ' ' && *str != '\0') str++;
        *str = '\0';
        int ans = fmpz_set_str(poly++, v, 62);
        if (ans) return -1;
        str++;
    }
    return 0;
}

inline int fmpz_poly_set_str_(fmpz_poly_t poly, char * str) {
    if (!isdigit((unsigned char) str[0])) return -1;
    slong len = atol(str);
    if (len < 0) return -1;
    if (len == 0) {
        fmpz_poly_zero(poly);
        return 0;
    }
    fmpz_poly_fit_length(poly, len);
    int ans = _fmpz_poly_set_str_(poly->coeffs, str);
    if (ans == 0) {
        _fmpz_poly_set_length(poly, len);
        _fmpz_poly_normalise(poly);
    } else {
        _fmpz_vec_zero(poly->coeffs, len);
        _fmpz_poly_set_length(poly, 0);
    }
    return ans;
}

inline void read_str3(fmpz_poly_q_struct* p, char* pos) {
    char *end = pos;
    while (*end != '/' && *end != '\0') ++end;
    if('/'==*end) { /* numerator / denominator */
        *end = '\0';
        fmpz_poly_set_str_(fmpz_poly_q_numref(p),pos);
        pos = end+1;
        fmpz_poly_set_str_(fmpz_poly_q_denref(p),pos);
    } else {
        fmpz_poly_set_str_(fmpz_poly_q_numref(p),pos);
        fmpz_poly_set_si(fmpz_poly_q_denref(p),1);
    }
}
#endif

// get monoms and coefficients from a point (Feynman integral). database access used
template<class I>
void p_get_internal(const point &p, I &terms, unsigned short dsector, unsigned short len, char* res) {

    char *pos = res + 3 + len * sizeof(point);

    for (unsigned int j = 0; j != len; ++j) {
        COEFF c;
#ifdef PRIME
        c.n = *reinterpret_cast<unsigned long long *>(pos);
        pos += sizeof(unsigned long long);
#elif defined(MPQStr)
        char *end = pos;
        while (*end != '|') ++end;
        *end = '\0';
        c.s.set_str(string(pos),62); // note the 62 base
        pos = end;
        ++pos;
#elif defined(MPQBin)
        c.s = mpq_read(pos);
#elif defined(FlintX)
        char *end = pos;
        while (*end != '|') ++end;
        *end = '\0';
        c.s = pos;
        pos = end;
        ++pos;
#elif defined(FlintC3)
        char *end = pos;
        while (*end != '|') ++end;
        *end = '\0';
        read_str3(c.s[0],pos);
        pos = end;
        ++pos;
#elif defined(FMPQ)
        char *end = pos;
        while (*end != '|') ++end;
        *end = '\0';
        fmpq_set_str(c.s[0],pos,62); // note the 62 base
        pos = end;
        ++pos;
#elif defined(FlintC)
        char *end = pos;
        while (*end != '|') ++end;
        *end = '\0';
        fmpz_poly_q_set_str2(c.s[0],pos);
        pos = end;
        ++pos;
#elif defined(FlintM)
        char *end = pos;
        while (*end != '|') ++end;
        *end = '\0';
        fmpz_mpoly_q_set_str(c.s[0],pos);
        pos = end;
        ++pos;
#else
        char *end = pos;
        while (*end != '|') ++end;
        *end = '\0';
        c.s = string(pos);
        pos = end;
        ++pos;
#endif
        terms.emplace_back(reinterpret_cast<point *>(res + 3)[j], c);
    }

    delete[] res;
    if ((len != 0) && (p != terms.back().first)) {
        cout << "p_get error" << endl;
        cout << p << endl;
        cout << terms.size() << endl;
        for (const auto & term : terms) {
            cout << term.first << endl;
        }
        abort();
    }
}

void p_set(const point &p, const vector<pair<point, COEFF> > &terms, unsigned char level, unsigned short fixed_database_sector) {
    p_set(p, static_cast<unsigned int>(terms.size()), terms.begin(), terms.end(), level, fixed_database_sector);
}

void p_set(const point &p, const list<pair<point, COEFF>, ALLOCATOR2 > &terms, unsigned char level, unsigned short fixed_database_sector) {
    p_set(p, static_cast<unsigned int>(terms.size()), terms.begin(), terms.end(), level, fixed_database_sector);
}

bool is_lower_in_orbit(const vector<t_index> &lhs, const vector<t_index> &rhs) {
    if (lhs == rhs) return false;
    vector<t_index> s1 = sector(lhs);
    vector<t_index> s2 = sector(rhs);
    if (s1 != s2) {
        unsigned short sn1 = common::sector_numbers_fast[sector_fast(s1)];
        unsigned short sn2 = common::sector_numbers_fast[sector_fast(s2)];
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
bool pair_point_COEFF_smaller(const pair<point, COEFF> &lhs, const pair<point, COEFF> &rhs) {
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
    unsigned short sn = common::sector_numbers_fast[ssector];
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
        unsigned short best_sn = sn;
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
            unsigned short new_sn = common::sector_numbers_fast[new_sector];
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
    unsigned short sn = common::sector_numbers_fast[ssector];
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

/* time to normalize an equation - to use the GCD everywhere and throw away zero members
* calls to fermat come from here
*/

#if defined(PolyMode)

void normalize(vector<pair<point, COEFF> > &terms, unsigned short thread_number) {
    vector<pair<point, COEFF> > mon;
    mon.reserve(terms.size());

    if (common::send_to_parent) {
        auto start_time = chrono::steady_clock::now();

        char lbuf[10];
        sprintf(lbuf, "%d\n", int(terms.size()));
        fputs(lbuf, common::child_stream_from_child);
        fflush(common::child_stream_from_child);
        char buf[sizeof(point) * 2 + 1];
        buf[sizeof(point) * 2] = '\0';
        int submitted = 0;


        auto itr = terms.begin();
        while (itr != terms.end()) {
            // we put the point
            itr->first.safe_string(buf);
            fputs(buf, common::child_stream_from_child);
            // we put the coefficient
            fputs(itr->second.s.c_str(), common::child_stream_from_child);
            fputs("\n", common::child_stream_from_child);
            fflush(common::child_stream_from_child);
            submitted++;
            ++itr;
        }

        while (submitted != 0) {
            int buf_size = 3;
            char *bbuf = static_cast<char *>(malloc(static_cast<size_t>(buf_size)));
            if (!bbuf) {
                cout<<"Cannot malloc in normalize"<<endl;
                abort();
            }

            read_from_stream(&bbuf, &buf_size, common::child_stream_to_child);

            bbuf[strlen(bbuf) - 1] = '\0';
            string ss = (bbuf + sizeof(point) * 2);
            if (ss != "0") {
                COEFF c;
                c.s = ss;
                mon.emplace_back(point(bbuf), c);
            }
            submitted--;
            free(bbuf);
        }
        auto stop_time = chrono::steady_clock::now();

        common::fermat_time += chrono::duration_cast<chrono::microseconds>(stop_time - start_time).count();

    } else {
        list<pair<int, pair<point, string> > > to_submit;
        int submitted = 0;
        auto itr = terms.begin();
        while (itr != terms.end()) {
            to_submit.emplace_back(thread_number, pair<point, string>(itr->first, itr->second.s));
            ++itr;
            submitted++;
        }

        {
            lock_guard<mutex> guard(equation::f_submit_mutex[thread_number % common::f_queues]);
            for (auto &f_job : to_submit) {
                equation::f_jobs[(thread_number % common::f_queues)].push_back(f_job);
            }
        }

        equation::f_submit_cond[(thread_number % common::f_queues)].notify_all();

        to_submit.clear();

        while (submitted != 0) {
            auto start_timeA = chrono::steady_clock::now();
            unique_lock<mutex> guard(equation::f_receive_mutex[thread_number]); // mutex is locked
            equation::f_receive_cond[thread_number].wait(guard,[thread_number](){return !equation::f_result[thread_number].empty();});
            // mutex is unlocked while waiting, then locked again
            pair<point, string> res = *(equation::f_result[thread_number].begin()); // take the result
            equation::f_result[thread_number].pop_front(); // remove result from queue
            guard.unlock(); // unlock guard, proceed to work;
            auto stop_timeA = chrono::steady_clock::now();
            common::thread_time -= chrono::duration_cast<chrono::microseconds>(stop_timeA - start_timeA).count();

            if (res.second != "0") {
                COEFF c;
                c.s = res.second;
                mon.emplace_back(res.first, c);
            }
            submitted--;
        }
    }
    sort(mon.begin(), mon.end(), pair_point_COEFF_smaller);
    terms = mon;
}


void normalize_eq(equation *eq, unsigned short thread_number) {
    if (common::send_to_parent) {
        auto start_time = chrono::steady_clock::now();
        char lbuf[10];
        sprintf(lbuf, "%d\n", int(eq->length));
        fputs(lbuf, common::child_stream_from_child);
        fflush(common::child_stream_from_child);
        char buf[sizeof(point) * 2 + 1];
        buf[sizeof(point) * 2] = '\0';
        int submitted = 0;

        MONOM **itr = eq->terms;
        while (itr != eq->terms + eq->length) {
            // we put the point
            (*itr)->p.safe_string(buf);
            fputs(buf, common::child_stream_from_child);
            // we put the coefficient
            fputs((*itr)->c.s.c_str(), common::child_stream_from_child);
            fputs("\n", common::child_stream_from_child);
            fflush(common::child_stream_from_child);
            submitted++;
            ++itr;
        }

        int new_length = 0;
        while (submitted != 0) {
            int buf_size = 3;
            char *bbuf = static_cast<char *>(malloc(static_cast<size_t>(buf_size)));
            if (!bbuf) {
                cout<<"Cannot malloc in normalize_eq"<<endl;
                abort();
            }

            read_from_stream(&bbuf, &buf_size, common::child_stream_to_child);
            bbuf[strlen(bbuf) - 1] = '\0';
            string ss = (bbuf + sizeof(point) * 2);
            if (ss != "0") {
                eq->terms[new_length][0].p = point(bbuf);
                eq->terms[new_length][0].c.s = ss;
                ++new_length;
            }
            submitted--;
            free(bbuf);
        }
        eq->length = new_length;
        auto stop_time = chrono::steady_clock::now();
        common::fermat_time += chrono::duration_cast<chrono::microseconds>(stop_time - start_time).count();

    } else {
        list<pair<int, pair<point, string> > > to_submit;
        int submitted = 0;
        MONOM **itr = eq->terms;
        while (itr != eq->terms + eq->length) {
            to_submit.emplace_back(thread_number, pair<point, string>((*itr)->p, (*itr)->c.s));
            ++itr;
            submitted++;
        }

        {
            lock_guard<mutex> guard(equation::f_submit_mutex[thread_number % common::f_queues]);
            for (const auto &f_job : to_submit) {
                equation::f_jobs[(thread_number % common::f_queues)].push_back(f_job);
            }
        }

        equation::f_submit_cond[(thread_number % common::f_queues)].notify_all();

        to_submit.clear();

        int new_length = 0;
        while (submitted != 0) {

            auto start_timeA = chrono::steady_clock::now();
            unique_lock<mutex> guard(equation::f_receive_mutex[thread_number]); // mutex is locked
            equation::f_receive_cond[thread_number].wait(guard,[thread_number](){return !equation::f_result[thread_number].empty();});
            pair<point, string> res = *(equation::f_result[thread_number].begin());
            equation::f_result[thread_number].pop_front();
            guard.unlock();
            auto stop_timeA = chrono::steady_clock::now();
            common::thread_time -= chrono::duration_cast<chrono::microseconds>(stop_timeA - start_timeA).count();

            if (res.second != "0") {
                eq->terms[new_length][0].p = res.first;
                eq->terms[new_length][0].c.s = res.second;
                ++new_length;
            }
            submitted--;
        }
        eq->length = new_length;
    }
}
#endif


// submit to fermat evaluation queue and wait for the result
void calc_wrapper(string &s, unsigned short thread_number) {
    if (s=="") return;
    if (common::send_to_parent) {
        fputs("1\n", common::child_stream_from_child);
        fflush(common::child_stream_from_child);
        char buf[sizeof(point) * 2 + 1];
        buf[sizeof(point) * 2] = '\0';
        point p;
        p.safe_string(buf);

        fputs(buf, common::child_stream_from_child);
        fputs(s.c_str(), common::child_stream_from_child);
        fputs("\n", common::child_stream_from_child);
        fflush(common::child_stream_from_child);

        int buf_size = 3;
        char *bbuf = static_cast<char *>(malloc(static_cast<size_t>(buf_size)));
        if (!bbuf) {
            cout<<"Cannot malloc in calc_wrapper"<<endl;
            abort();
        }

        read_from_stream(&bbuf, &buf_size, common::child_stream_to_child);
        bbuf[strlen(bbuf) - 1] = '\0';
        s = (bbuf + sizeof(point) * 2);
        free(bbuf);
    } else {
        {
            lock_guard<mutex> guard(equation::f_submit_mutex[thread_number % common::f_queues]);
            equation::f_jobs[(thread_number % common::f_queues)].emplace_back(thread_number, pair<point, string>(point(), s));
        }

        equation::f_submit_cond[(thread_number % common::f_queues)].notify_one(); // one because there is only one evaluation task



        unique_lock<mutex> guard(equation::f_receive_mutex[thread_number]);
        equation::f_receive_cond[thread_number].wait(guard,[thread_number](){return !equation::f_result[thread_number].empty();});
        pair<point, string> res = *(equation::f_result[thread_number].begin());
        equation::f_result[thread_number].pop_front();
        guard.unlock();
        s = res.second;
    }
}


void f_worker(unsigned short fnum, unsigned short qnum) {

    while (true) {
        unique_lock<mutex> guard(equation::f_submit_mutex[qnum]); // it's locked
        equation::f_submit_cond[qnum].wait(guard,[qnum](){return (equation::f_stop || !equation::f_jobs[qnum].empty());});
        // predicate is checked on locked mutex. if we are out, is is still locked, but while we are waiting, it's not
        if (equation::f_stop) {
            break; // time to stop, and the mutex get's unlocked on contex end
        }
        int t_number = equation::f_jobs[qnum].begin()->first;
        pair<point, string> f_submit = equation::f_jobs[qnum].begin()->second;
        equation::f_jobs[qnum].pop_front(); // we take the data out of the queue
        guard.unlock(); // mutex is unlocked for other threads, and we start working

        auto start_time = chrono::steady_clock::now();
        f_submit.second.erase(std::remove(f_submit.second.begin(), f_submit.second.end(), ' '), f_submit.second.end());
        f_submit.second = calc(f_submit.second, fnum);

        if ((f_submit.second[0] == '0') && (f_submit.second[1] == '/')) {
            f_submit.second = "0";
        }

        auto stop_time = chrono::steady_clock::now();
        common::fermat_time += chrono::duration_cast<chrono::microseconds>(stop_time - start_time).count();

        {
            lock_guard<mutex> receive_guard(equation::f_receive_mutex[t_number]);
            equation::f_result[t_number].push_back(f_submit);
        }
        equation::f_receive_cond[t_number].notify_one(); // only one should be waiting
    }
}
