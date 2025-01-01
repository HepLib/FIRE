#pragma once
#include <gmp.h>

#include "flint/fmpq.h"
#include "flint/fmpz_mpoly_q.h"
#include "flint/mpoly.h"
#include "flint/gr.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <map>
#include <list>
#include <set>
#include <vector>
#include <unistd.h>
#include <algorithm>
#include <cstdint>
#include <fcntl.h>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <shared_mutex>
#include "taskspool.h"

using namespace std;

typedef unsigned long long prime_type;
string calc(const string &buf); // from gateToFermat

extern prime_type primes[128];
class COEFF;
extern COEFF CO_0;
extern COEFF CO_1;
extern COEFF CO_1m;
typedef __int128 int128_t;
typedef unsigned __int128 uint128_t;

#if defined(FlintM)
#define TO_CTX(ring_ctx) ((fmpz_mpoly_ctx_struct *)(GR_CTX_DATA_AS_PTR(ring_ctx)))
inline string get_str(const fmpz_mpoly_q_struct* mp, const char** xs, fmpz_mpoly_ctx_t ctx) {
    auto cs = fmpz_mpoly_get_str_pretty(fmpz_mpoly_q_numref(mp), xs, ctx);
    string s(cs);
    flint_free(cs);
    if(!fmpz_mpoly_is_one(fmpz_mpoly_q_denref(mp), ctx)) {
        cs = fmpz_mpoly_get_str_pretty(fmpz_mpoly_q_denref(mp), xs, ctx);
        string ds(cs);
        flint_free(cs);
        bool is_nc = fmpz_mpoly_is_fmpz(fmpz_mpoly_q_numref(mp), ctx);
        bool is_dc = fmpz_mpoly_is_fmpz(fmpz_mpoly_q_denref(mp), ctx);
        if(!is_nc) s = "("+s+")";
        if(!is_dc) ds = "("+ds+")";
        s = s+"/"+ds;
    }
    return s;
}
inline void set_str(fmpz_mpoly_q_struct* mp, char* pos, const char** xs, fmpz_mpoly_ctx_t ctx) {
    char *end = pos;
    while (*end != '/' && *end != '\0' && *end != '|') ++end;
    if('/'==*end || '|'==*end) { /* numerator/denominator & numerator|denominator*/ 
        *end = '\0';
        int nok = fmpz_mpoly_set_str_pretty(fmpz_mpoly_q_numref(mp), pos, xs, ctx);
        if(nok) {
            cout << "not supported polynormial: " << pos << endl;
            abort();
        }
        pos = end+1;
        nok = fmpz_mpoly_set_str_pretty(fmpz_mpoly_q_denref(mp), pos, xs, ctx);
        if(nok) {
            cout << "not supported polynormial: " << pos << endl;
            abort();
        }
    } else {
        int nok = fmpz_mpoly_set_str_pretty(fmpz_mpoly_q_numref(mp), pos, xs, ctx);
        if(nok) {
            cout << "not supported polynormial: " << pos << endl;
            abort();
        }
        fmpz_mpoly_set_si(fmpz_mpoly_q_denref(mp), 1, ctx);
    }
}
inline void set_str(fmpz_mpoly_q_struct* mp, const char* pos, const char** xs, fmpz_mpoly_ctx_t ctx) {
    int n = strlen(pos);
    char *buff = new char[n+1];
    strcpy(buff, pos);
    set_str(mp, buff, xs, ctx);
    delete []buff;
}
#endif

/** Used for storing indices of integrals */
typedef char t_index;

#ifndef SMALL_POINT
constexpr size_t MAX_IND = {22};
#else
constexpr size_t MAX_IND = {18};
constexpr size_t BITS_PER_INDEX = {6};
#endif
// 22 indices should be enough for most problems;
// 5-loop propagator = 20;
// 6-loop bubble = 21;


typedef uint16_t sector_count_t;
// Typedef for variables that involve the number of sectors.
// point::ww and database functions have compile-time constants defined
// based on the size of the type
constexpr sector_count_t MAX_SECTORS  = {1<<15};
// Number of sectors should fit into 2^15.
// It is possible to increase the allowed MAX_SECTORS by changing
// the associated sector_count_t to uint32_t.  However, because
// there are a number of statically-allocated data structures
// sized by MAX_SECTORS, there are potential linker/assembler problems
// when using values larger than ~2^21.  Any possible need to use
// values greater than 2^31 will also require some type promotions
// in the signatures of functions::run/watch_child (which take negative
// sector numbers during back substitution), and for a few variables
// in functions::Evaluate

typedef uint32_t SECTOR;
//static set<unsigned short> needed_for[MAX_SECTORS + 1];

typedef map<sector_count_t, vector<pair<vector<pair<vector<t_index>, pair<short, bool>>>, vector<pair<string, vector<pair< vector<t_index>, short>>>>>>> common_lbases_t;

int positive_index(vector<t_index> v);

vector<t_index> _sector_(const vector<t_index> &v);
vector<t_index> corner(const vector<t_index> &v);
vector<t_index> degree(const vector<t_index> &v);
pair<unsigned int, unsigned int> level(const vector<t_index> &v);
void print_vector(const vector<t_index> &v);
void print_sector_fast(const SECTOR &sf);
void Orbit(const vector<t_index> &v, vector<vector<t_index> > &orbit, vector<vector<vector<t_index> > > &sym);
void make_ordering(t_index *mat, const vector<t_index> &sector);
vector<vector<t_index> > all_sectors(unsigned int d, int positive, int positive_start);
SECTOR sector_fast(const vector<t_index> &v);

class common {
public:
    static map<string, bool> opt_set; // check the option has been set in command line argument
    static map<string, string> prt_replace;
    static int prt_rule_counter;
    const static string version;
    static int run_sector; // int not SECTOR
    static int verb;
    static vector<int> sector_tasks; // not sector_count_t, due to negative sector
    static bool skip_if_exist;
    static int run_mode;
    static bool o_output;
    static string dep_file;
    static int pos_pref;
    static unique_ptr<sector_count_t[]> sector_numbers_fast;
    static unique_ptr<unique_ptr<t_index[]>[]> orderings_fast;
    static unsigned short dimension;
    static bool disable_presolve;
    static string suffix;
    static std::map<string, string> var_values_from_arv;
    static SECTOR virtual_sector;
    static int abs_min_level;
    static int abs_max_level;
    static sector_count_t abs_max_sector;
    static bool split_masters;
    static unsigned int master_number_min;
    static unsigned int master_number_max;
    static bool only_masters;
    static bool silent;
    static bool keep_all;
    static string fermat;
    static string variables;
    static string config_file;
    static string path; // Path to databases.
    static string hint_path;
    static bool all_ibps;
    static vector<vector<vector<t_index> > > iorderings;
    static vector<vector<t_index> > ssectors;
    static vector<vector<vector<t_index> > >  symmetries;
    static set<SECTOR> lsectors;

    static common_lbases_t lbases;
    static unsigned int global_pn;

    static unsigned int t1; // threads @ Forward
    static unsigned int t2; // threads @ Backward
    static unsigned int lt1; // level threads
    static unsigned int lt2; // level threads
    static unsigned int tp; // pool size
    static unsigned int lmt1; // Level Limit @ Forward
    static unsigned int lmt2;
    static unsigned int len; // len limit in Tasks Pools
    static TasksPool TPool;
    static int ifm; // 0-auto, 1-prime, 2-poly
    static map<string, string> variable_replacements;
    static bool small;
    static bool hint;
     
};

bool in_lsectors(sector_count_t test_sector);

int128_t mul_inv(int128_t a, int128_t b);
prime_type mod(const string &s, bool positive_number_in_range = false);

string ReplaceAll(string str, const string &from, const string &to);
string ReplaceAllVariables(string str);

class COEFF {
public:
    gr_ptr r = NULL;
    prime_type p = 0;
    
    inline static bool ctx_init = false;
    inline static bool is_float = false;
    inline static gr_ctx_t ctx;
    inline static uint64_t prime = 0;
    inline static unsigned short prime_number = 0;
    inline static const char ** xs = NULL;
    inline static vector<string> vs;
    inline static unsigned int fp = 500;
    
    void static init_ctx_m(int n) {
        if(ctx_init) return;
        ctx_init = true;
        gr_ctx_init_fmpz_mpoly_q(ctx, n, ORD_LEX);
        CO_0.r = gr_heap_init(ctx);
        gr_zero(CO_0.r, ctx);
        CO_1.r = gr_heap_init(ctx);
        gr_one(CO_1.r, ctx);
        CO_1m.r = gr_heap_init(ctx);
        gr_neg(CO_1m.r, CO_1.r, ctx);
    }
    
    void static init_ctx_q() {
        if(ctx_init) return;
        ctx_init = true;
        gr_ctx_init_fmpq(ctx);
        CO_0.r = gr_heap_init(ctx);
        gr_zero(CO_0.r, ctx);
        CO_1.r = gr_heap_init(ctx);
        gr_one(CO_1.r, ctx);
        CO_1m.r = gr_heap_init(ctx);
        gr_neg(CO_1m.r, CO_1.r, ctx);
    }
    
    void static init_ctx_f() {
        if(ctx_init) return;
        ctx_init = true;
        is_float = true;
        gr_ctx_init_real_float_arf(ctx, fp);
        CO_0.r = gr_heap_init(ctx);
        gr_zero(CO_0.r, ctx);
        CO_1.r = gr_heap_init(ctx);
        gr_one(CO_1.r, ctx);
        CO_1m.r = gr_heap_init(ctx);
        gr_neg(CO_1m.r, CO_1.r, ctx);
    }
    
    void static init_prime(int pn) {
        prime_number = pn;
        prime = primes[prime_number];
        CO_0.p = 0;
        CO_1.p = 1;
        CO_1m.p = prime-1;
    }
        
    inline void init() {
        r = gr_heap_init(ctx);
    }

    explicit COEFF(slong n=0) {
        if(prime) {
            int128_t nn = n;
            nn = nn % prime;
            if(nn<0) nn += prime;
            p = nn;
        }
        if(ctx_init) {
            init();
            #if defined(FMPQ)
            fmpq_set_si((fmpq*)r, n, 1);
            #elif defined(FlintM)
            fmpz_mpoly_q_set_si((fmpz_mpoly_q_struct*)r, n, TO_CTX(ctx));
            #else
            gr_set_si(r, n, ctx);
            #endif
        }
    }
    
    size_t length() const {
        #if defined(PRIME) || defined(FMPQ) || defined(FloatR)
        return 1;
        #elif defined(FlintM)
        fmpz_mpoly_q_struct* mp = (fmpz_mpoly_q_struct*)r;
        size_t num_len = fmpz_mpoly_length(fmpz_mpoly_q_numref(mp), TO_CTX(ctx));
        size_t den_len = fmpz_mpoly_length(fmpz_mpoly_q_denref(mp), TO_CTX(ctx));
        return num_len>den_len ? num_len : den_len;
        #endif
    }
    
    void init_p(const string str, int base=10) {
        string pstr = str;
        for (int i=0; i<vs.size(); i++) {
            pstr = ReplaceAll(pstr, vs[i], "("+to_string(11+i*100)  + "/349374932742332831)");
        }
        if(vs.size()>0 && base==10) {
            pstr = calc(pstr); // note we need calc below, since pstr may be overflow in llu read
            pstr = calc("(Numer("+pstr+")|"+to_string(prime)+")/(Denom("+pstr+")|"+to_string(prime)+")");
            p = mod(pstr, true); // note true here
        } else p = mod(pstr, base==62); // base==62 -> positive_in_range
    }
    void init_r(const string str, int base=10) {
        #if defined(FMPQ)
        int nok = fmpq_set_str((fmpq*)r, str.c_str(), base);
        if(nok) throw runtime_error("fmpq_set_str error, str='" + str + "'");
        #elif defined(FlintM)
        set_str((fmpz_mpoly_q_struct*)r, str.c_str(), xs, TO_CTX(ctx));
        #elif defined(FloatR)
        fmpq_t q;
        fmpq_init(q);
        int nok = fmpq_set_str(q, str.c_str(), base);
        if(nok) throw runtime_error("fmpq_set_str error, str='" + str + "'");
        gr_set_fmpq(r, q, ctx);
        fmpq_clear(q);
        #endif
    }
    explicit COEFF(const string str, int base=10) {
        if(prime) init_p(str, base);
        if(ctx_init) { init(); init_r(str, base); }
    }
    
    explicit COEFF(const COEFF & o) {
        p = o.p;
        if(ctx_init) { init(); gr_set(r, o.r, ctx); }
    }
    
    COEFF & operator=(const COEFF & o) {
        p = o.p;
        if(ctx_init) gr_set(r, o.r, ctx);
        return *this;
    }
    
    COEFF(COEFF && o) {
        p = o.p;
        if(ctx_init) { init(); gr_swap(r, o.r, ctx); }
    }

    COEFF & operator=(COEFF && o) {
        p = o.p;
        if(ctx_init) gr_swap(r, o.r, ctx);
        return *this;
    }

    ~COEFF() {
        if(ctx_init) if(r!=NULL) gr_heap_clear(r, ctx);
    }
    
    inline bool is_zero_p() const { return !p; }
    inline bool is_zero_r() const { return (gr_is_zero(r, ctx)==truth_t::T_TRUE); }
    bool is_zero() const {
        if(is_float) {
            bool z = is_zero_p();
            if(!z) return false;
            // TODO: add a check to test numerical value
            return true;
        }
        if(prime && ctx_init) return is_zero_p() && is_zero_r();
        else if(prime) return is_zero_p();
        else return is_zero_r();
    }
    
    #if defined(FloatR)
    inline static int io_lines = 3;
    #else
    inline static int io_lines = 1;
    #endif
    void export_to(ostream & oss, int base) const { // total lines: io_lines
        #if defined(PRIME)
        oss << static_cast<uint64_t>(p) << '\n';
        #elif defined(FMPQ)
        char * cs = fmpq_get_str(NULL,base,(fmpq*)r);
        oss << cs << '\n';
        flint_free(cs);
        #elif defined(FlintM)
        fmpz_mpoly_q_struct* mp = (fmpz_mpoly_q_struct*)r;
        auto cs = fmpz_mpoly_get_str_pretty(fmpz_mpoly_q_numref(mp), xs, TO_CTX(ctx));
        oss << cs;
        flint_free(cs);
        if(!fmpz_mpoly_is_one(fmpz_mpoly_q_denref(mp), TO_CTX(ctx))) {
            cs = fmpz_mpoly_get_str_pretty(fmpz_mpoly_q_denref(mp), xs, TO_CTX(ctx));
            oss << "|" << cs;
            string ds(cs);
            flint_free(cs);
        }
        oss << '\n';
        #elif defined(FloatR)
        fmpz_t z1, z2;
        fmpz_init(z1);
        fmpz_init(z2);
        gr_get_fmpz_2exp_fmpz(z1, z2, r, ctx);
        auto s = fmpz_get_str(NULL, base, z1);
        oss << s << '\n';
        flint_free(s);
        s = fmpz_get_str(NULL, base, z2);
        oss << s << '\n';
        flint_free(s);
        fmpz_clear(z1);
        fmpz_clear(z2);
        oss << static_cast<uint64_t>(p) << '\n'; // export prime
        #endif
    }
    void import_from(fstream & db, int base) {
        string s;
        db >> s;
        #if defined(FloatR)
        fmpz_t z1, z2;
        fmpz_init(z1);
        fmpz_init(z2);
        fmpz_set_str(z1, s.c_str(), base);
        db >> s;
        fmpz_set_str(z2, s.c_str(), base);
        gr_set_fmpz_2exp_fmpz(r, z1, z2, ctx);
        fmpz_clear(z1);
        fmpz_clear(z2);
        db >> s; // import prime
        p = mod(s, base==62);
        #elif defined(PRIME)
        p = mod(s, base==62);
        #else
        init_r(s, base);
        #endif
    }
};

inline bool is_zero_p(const COEFF & c) { return c.is_zero_p(); }
inline bool is_zero_r(const COEFF & c) { return c.is_zero_r(); }
inline bool is_zero(const COEFF & c) { return c.is_zero(); }

inline void _add_(prime_type & pr, const prime_type & p1, const prime_type & p2) {
    uint128_t n1 = p1;
    uint128_t n2 = p2;
    n1 += n2;
    n1 = n1 % COEFF::prime;
    pr = n1;
}
inline void _add_(COEFF &cr, const COEFF &c1, const COEFF &c2) {
    if(COEFF::prime) _add_(cr.p, c1.p, c2.p);
    if(COEFF::ctx_init) gr_add(cr.r, c1.r, c2.r, COEFF::ctx);
}

inline void _sub_(prime_type & pr, const prime_type & p1, const prime_type & p2) {
    uint128_t n1 = p1;
    uint128_t n2 = p2;
    n2 = COEFF::prime - n2;
    n1 += n2;
    n1 = n1 % COEFF::prime;
    pr = n1;
}
inline void _sub_(COEFF &cr, const COEFF &c1, const COEFF &c2) {
    if(COEFF::prime) _sub_(cr.p, c1.p, c2.p);
    if(COEFF::ctx_init) gr_sub(cr.r, c1.r, c2.r, COEFF::ctx);
}

inline void _mul_(prime_type & pr, const prime_type & p1, const prime_type & p2) {
    uint128_t n1 = p1;
    uint128_t n2 = p2;
    n1 *= n2;
    n1 = n1 % COEFF::prime;
    pr = n1;
}
inline void _mul_(COEFF &cr, const COEFF &c1, const COEFF &c2) {
    if(COEFF::prime) _mul_(cr.p, c1.p, c2.p);
    if(COEFF::ctx_init) gr_mul(cr.r, c1.r, c2.r, COEFF::ctx);
}

inline void _mul_(prime_type &pr, const prime_type &p1, slong i) {
    uint128_t n1 = p1;
    int128_t n2 = i;
    n2 = n2 % COEFF::prime;
    if(n2<0) n2 += COEFF::prime;
    n1 *= n2;
    n1 = n1 % COEFF::prime;
    pr = n1;
}
inline void _mul_(COEFF &cr, const COEFF &c1, slong i) {
    if(COEFF::prime) _mul_(cr.p, c1.p, i);
    if(COEFF::ctx_init) gr_mul_si(cr.r, c1.r, i, COEFF::ctx);
}

inline void _div_(prime_type & pr, const prime_type & p1, const prime_type & p2) {
    uint128_t nn = p1;
    uint128_t nd = p2;
    nd = mul_inv(nd, COEFF::prime);
    nn *= nd;
    nn = nn % COEFF::prime;
    pr = nn;
}
inline void _div_(COEFF &cr, const COEFF &c1, const COEFF &c2) {
    if(COEFF::prime) _div_(cr.p, c1.p, c2.p);
    if(COEFF::ctx_init) gr_div(cr.r, c1.r, c2.r, COEFF::ctx);
}

inline void _neg_(prime_type & p) {
    if(p) p = COEFF::prime - p;
}
inline void _neg_(COEFF & o) {
    if(COEFF::prime) if(o.p) o.p = COEFF::prime - o.p;
    if(COEFF::ctx_init) gr_neg(o.r, o.r, COEFF::ctx);
}

inline void _swap_(prime_type & p1, prime_type & p2) {
    std::swap(p1, p2);
}
inline void _swap_(COEFF & o1, COEFF & o2) {
    if(COEFF::prime) std::swap(o1.p, o2.p);
    if(COEFF::ctx_init) gr_swap(o1.r, o2.r, COEFF::ctx);
}

inline void _set_(COEFF & o1, const COEFF & o2) {
    o1.p = o2.p;
    if(COEFF::ctx_init) gr_set(o1.r, o2.r, COEFF::ctx);
}
inline void _set_(COEFF & o1, COEFF && o2) {
    o1.p = o2.p;
    if(COEFF::ctx_init) gr_swap(o1.r, o2.r, COEFF::ctx);
}

inline void _set_(prime_type & p, slong i) {
    int128_t nn = i;
    nn = nn % COEFF::prime;
    if(nn<0) nn += COEFF::prime;
    p = nn;
}
inline void _set_(COEFF & o, slong i) {
    if(COEFF::prime) _set_(o.p, i);
    if(COEFF::ctx_init) gr_set_si(o.r, i, COEFF::ctx);
}

inline void _add_mul_(prime_type &pr, const prime_type &p1, const prime_type &p2) {
    prime_type p;
    _mul_(p, p1, p2);
    _add_(pr, pr, p);
}
inline void _add_mul_(COEFF &cr, const COEFF &c1, const COEFF &c2) {
    if(COEFF::prime) _add_mul_(cr.p, c1.p, c2.p);
    if(COEFF::ctx_init) gr_addmul(cr.r, c1.r, c2.r, COEFF::ctx);
}

inline void _sub_mul_(prime_type &pr, const prime_type &p1, const prime_type &p2) {
    prime_type p;
    _mul_(p, p1, p2);
    _sub_(pr, pr, p);
}
inline void _sub_mul_(COEFF &cr, const COEFF &c1, const COEFF &c2) {
    if(COEFF::prime) _sub_mul_(cr.p, c1.p, c2.p);
    if(COEFF::ctx_init) gr_submul(cr.r, c1.r, c2.r, COEFF::ctx);
}

inline void _mul_neg_(prime_type & pr, const prime_type &p1, const prime_type &p2) {
    _mul_(pr, p1, p2);
    _neg_(pr);
}
inline void _mul_neg_(COEFF & cr, const COEFF &c1, const COEFF &c2) {
    _mul_(cr, c1, c2);
    _neg_(cr);
}

inline void _div_neg_(prime_type & pr, const prime_type &p1, const prime_type &p2) {
    _div_(pr, p1, p2);
    _neg_(pr);
}
inline void _div_neg_(COEFF & cr, const COEFF &c1, const COEFF &c2) {
    _div_(cr, c1, c2);
    _neg_(cr);
}

inline ostream & operator<<(ostream & os, const COEFF &o) {
    #if defined(PRIME)
    os << static_cast<uint64_t>(o.p);
    #elif defined(FMPQ)
    char * cs;
    cs = fmpq_get_str(NULL, 10, (fmpq*)o.r);
    os << cs;
    flint_free(cs);
    #elif defined(FlintM)
    string res = get_str((fmpz_mpoly_q_struct*)o.r, COEFF::xs, TO_CTX(COEFF::ctx));
    os << res;
    #else
    char *cs;
    gr_get_str(&cs, o.r, COEFF::ctx);
    os << cs;
    flint_free(cs);
    #endif
    return os;
}
