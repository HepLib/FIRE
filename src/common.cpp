/** @file common.cpp
 *  @author Alexander Smirnov
*
*  This file is a part of the FIRE package.
*  It contains the initializations of static variables in the common class
*  and a number of general functions used in other parts of the program.
*/


#include "common.h"
#include <sys/stat.h>
#include <sys/mman.h>
#include <omp.h>
/* initializations of static variables
* refer to common.h for details
*/

#include "gateToFermat.h"

unsigned short common::dimension;
unique_ptr<sector_count_t[]> common::sector_numbers_fast;
unique_ptr<unique_ptr<t_index[]>[]> common::orderings_fast;

const string common::version = "2024/05/13 (based on v6.5.2)";
map<string, string> common::prt_replace;
int common::prt_rule_counter = 0;
bool common::skip_if_exist = false;
int common::run_mode = 0; // 0: not use db_file, 1: use db_file 2: restart from db_file
bool common::o_output = false; // use offical FIRE output
string common::dep_file{};
int common::pos_pref = 1;
bool common::silent = false;
bool common::disable_presolve = false;
bool common::hint;
map<string,bool> common::opt_set;
string common::hint_path;
string common::fermat;
string common::variables;
string common::config_file;

map<string, string> common::variable_replacements; // while parsing

bool common::all_ibps = false;
bool common::split_masters = false;
unsigned int common::master_number_min = 0;
unsigned int common::master_number_max = 0;

SECTOR common::virtual_sector = 0;

string common::path = "";

string common::suffix = "";
std::map<string, string> common::var_values_from_arv;

vector<int> common::sector_tasks; // not sector_count_t, due to negative sector
int common::run_sector = 0; // int not SECTOR
bool common::only_masters = false;
bool common::keep_all = false;

sector_count_t common::abs_max_sector = 3;

int common::abs_max_level = 0;

int common::abs_min_level = 100;

TasksPool common::TPool;
unsigned int common::global_pn = 0;
int common::ifm = 0;
unsigned int common::t1 = 0; // 0 for all CPU cores
unsigned int common::t2 = 0; // 0 for all CPU cores
unsigned int common::lt1 = 1;
unsigned int common::lt2 = 1;
#if defined(FMPQ) || defined(FloatR)
unsigned int common::lmt1 = 1000; // Level Limit
unsigned int common::lmt2 = 1000; // Level Limit
unsigned int common::len = 10; // no effect
unsigned int common::tp = 4; // pool size
#else
unsigned int common::lmt1 = 100; // Level Limit
unsigned int common::lmt2 = 100; // Level Limit
unsigned int common::len = 50;
unsigned int common::tp = omp_get_num_procs()/2; // pool size
#endif

vector<vector<vector<t_index> > > common::iorderings;
vector<vector<vector<t_index> > > common::symmetries;

vector<vector<t_index> > common::ssectors;
set<SECTOR> common::lsectors;

common_lbases_t common::lbases;
//      sector    condition+result        coefficients  free term   eq                 coeff     indices:   coefficients and free term

vector<t_index> _sector_(const vector<t_index> &v) {
    vector<t_index> result;
    for (const t_index i : v)
        if (i > 0) {
            result.push_back(1);
        } else {
            result.push_back(-1);
        }
    return result;
}

SECTOR sector_fast(const vector<t_index> &v) {
    SECTOR result = 0;
    for (const t_index i : v)
        if (i > 0) {
            result = ((result << 1) ^ 1);
        } else {
            result = result << 1;
        }
    return result;
}

vector<t_index> corner(const vector<t_index> &v) {
    vector<t_index> result;
    for (const t_index i : v)
        if (i > 0) {
            result.push_back(1);
        } else {
            result.push_back(0);
        }
    return result;
}


vector<t_index> degree(const vector<t_index> &v) {
    vector<t_index> result;
    for (const t_index i : v)
        if (i > 0) {
            result.push_back(i - 1);
        } else {
            result.push_back(-i);
        }
    return result;
}


pair<unsigned int, unsigned int> level(const vector<t_index> &v) {
    unsigned int p = 0;
    unsigned int m = 0;
    for (const t_index i : v)
        if (i > 0) {
            p += (i - 1);
        } else {
            m += (-i);
        }
    return make_pair(p,m);
}


// output routines for vectors and double vectors
void print_vector(const vector<t_index> &v) {
    if (common::silent) return;
    auto it = v.begin();
    cout << "{" << int(*it);
    for (it++; it != v.end(); it++) {
        cout << "," << int(*it);
    }
    cout << "}";
}

// output routines for vectors and double vectors
void print_sector_fast(const SECTOR &sf) {
    if (common::silent) return;
    SECTOR i = 1;
    i <<= (common::dimension -1);
    cout << "{" << (i&sf ? 1 : -1);
    for (i>>=1; i; i>>=1) {
        cout << "," << (i&sf ? 1 : -1);
    }
    cout << "}";
}


void Orbit(const vector<t_index> &v, vector<vector<t_index> > &orbit, vector<vector<vector<t_index> > > &sym) {
    for (const auto &values : sym) {
        if (values.size() == 4) {
            const vector<t_index> &restrictions = values[2];
            for (unsigned int i = 0; i != v.size(); ++i) {
                if ((restrictions[i] == 1) && (v[i] <= 0)) continue;
                if ((restrictions[i] == -1) && (v[i] > 0)) continue;
            }
        }
        const vector<t_index> &permutation = values[0];
        vector<t_index> result;
        for (unsigned int i = 0; i != v.size(); ++i) {
            result.push_back(v[permutation[i] - 1]);
        }
        bool to_add = true;
        for (unsigned int i = (*(values.rbegin()))[0]; i != result.size(); ++i) {
            if (result[i] > 0) {
                to_add = false;
            }
        }
        if (to_add) orbit.push_back(result);
    }
}

void make_ordering(t_index *mat, const vector<t_index> &sector) {
    vector<vector<t_index> > result;
    result.reserve(sector.size());
    vector<t_index> neg;
    vector<t_index> one_vector;
    vector<t_index> zero_vector;
    zero_vector.reserve(sector.size());

    bool all_plus = true;
    bool all_minus = true;
    for (const t_index i : sector) {
        if (i != 1) all_plus = false;
        if (i != -1) all_minus = false;
    }

    auto l_pos = static_cast<unsigned int>(-1);
    auto l_neg = static_cast<unsigned int>(-1);

    for (unsigned int i = 0; i < sector.size(); ++i) {
        one_vector.push_back(1);
        zero_vector.push_back(0);
        if (sector[i] == 1) {
            l_pos = i;
            neg.push_back(0);
        } else {
            l_neg = i;
            neg.push_back(1);
        }
    }

    if ((all_plus) || all_minus) {
        result.push_back(one_vector);
        for (unsigned int i = 0; i + 1 < sector.size(); ++i) {
            vector<t_index> v;
            for (unsigned int j = 0; j < sector.size(); ++j) {
                if (i == j) {
                    v.push_back(1);
                } else {
                    v.push_back(0);
                }
            }
            result.push_back(v);
        }
    } else {
        if (l_pos == static_cast<unsigned int>(-1) || l_neg == static_cast<unsigned int>(-1)) {
            cout << "Something wrong with make_ordering!" << endl;
            abort();
        }
        result.push_back(one_vector);
        result.push_back(neg);
        for (unsigned int i = 0; i < sector.size() - 1; ++i) {
            if (sector[i] == 1) {
                if (i != l_pos) {
                    vector<t_index> v = zero_vector;
                    v[i] = 1;
                    result.push_back(v);
                }
            } else {
                if (i != l_neg) {
                    neg[i] = 0;
                    result.push_back(neg);
                }
            }
        }
    }

    int row = 0;
    for (auto itr_row = result.begin(); itr_row != result.end(); ++itr_row, ++row) {
        // now we stopped using vectors, moving to arrays, it's much faster
        // and array of vectors leads to crashes
        int column = 0;
        for (auto itr_column = itr_row->begin(); itr_column != itr_row->end(); ++itr_column, ++column) {
            mat[row * common::dimension + column] = *itr_column;
        }
    }
}


vector<vector<t_index> > all_sectors(unsigned int d, int positive, int positive_start) {
    vector<t_index> vv3;
    for (unsigned int i = 0; i != d; ++i) vv3.push_back(-1);
    vector<vector<t_index> > vv;
    vv.push_back(vv3);
    for (int i = positive_start - 1; i < positive; ++i) {
        vector<vector<t_index> > ww = vv;
        vv.reserve(2 * vv.size());
        for (auto &current : ww) {
            current[i] = 1;
            vv.push_back(current);
        };
    };
    return (vv);
};

bool in_lsectors(sector_count_t test_sector) {
    if (common::lsectors.empty()) return true;
    return (common::lsectors.find(sector_fast(common::ssectors[test_sector]))!= common::lsectors.end());
}

int positive_index(vector<t_index> v) {
    int res = 0;
    for (unsigned int i = 0; i != v.size(); ++i) {
        if (v[i] == 1) res++;
    }
    return res;
}

int128_t mul_inv(int128_t a, int128_t b) {
    int128_t b0 = b, t, q;
    int128_t x0 = 0, x1 = 1;
    if (b == 1) {
        return 1;
    }
    while (a > 1) {
        q = a / b;
        t = b, b = a % b, a = t;
        t = x0, x0 = x1 - q * x0, x1 = t;
    }
    if (x1 < 0) x1 += b0;
    return x1;
}

unsigned long long mod(const string &s, bool positive_number_in_range) {
    uint128_t num, denum;
    int128_t temp128;
    long long int temp;
    unsigned long long int utemp;
    if (positive_number_in_range) {
        sscanf(s.c_str(), "%llu", &utemp);
        num = utemp;
    } else {
        sscanf(s.c_str(), "%lld", &temp); // there should be no overflow here, just reading
        temp128 = temp;
        temp128 = temp128 % (int128_t(COEFF::prime));
        if (temp128 < 0) {
            temp128 += (int128_t(COEFF::prime));
        }
        num = temp128;
    }

    unsigned long pos = s.find('/');
    if (pos != string::npos) {
        if (positive_number_in_range) {
            sscanf(s.c_str() + pos + 1, "%llu", &utemp);
            denum = utemp;
        } else {
            sscanf(s.c_str() + pos + 1, "%lld", &temp);
            if (temp < 0) {
                cout << "Bad denominator: " << s << endl;
                abort();
                // this cannot happen. we are reading a small number and the denominator cannot be big or negative
            }
            temp128 = temp;
            temp128 = temp128 % int128_t(COEFF::prime);
            denum = temp128;
        }

        denum = mul_inv(denum, COEFF::prime);
        num *= denum;
        num %= COEFF::prime;
    }
    return num;
}

string ReplaceAll(string str, const string &from, const string &to) {
    size_t start_pos = 0;
    while ((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
    }
    return str;
}

string ReplaceAllVariables(string str) {
    for (const auto &variable_replacement : common::variable_replacements) {
        str = ReplaceAll(str, variable_replacement.first, "(" + variable_replacement.second + ")");
    }
    return str;
}

COEFF CO_0(0);
COEFF CO_1(1);
COEFF CO_1m(-1); // depend on common::prime
