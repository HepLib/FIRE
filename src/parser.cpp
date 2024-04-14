/** @file parser.cpp
 *  @author Alexander Smirnov
 *
 *  This file is a part of the FIRE package.
*/


#include "parser.h"
#include <math.h>
#include <filesystem>

inline bool file_remove(string fn) {
    return remove(fn.c_str())==0;
}

constexpr size_t LOAD_STR_SIZE = {1024}; ///<size of the string loaded from files
constexpr size_t COEFF_BUF_SIZE = {1024}; ///<size of buffer for parsing coefficients

extern unsigned long long primes[128];

/**
 * String to store temporary relative folder paths.
 */
string folder;

/**
 * Matrix of sector needed by sector relations
 */
static vector<vector<bool> > dependencies{};
/**
 * Sectors that contain non-zero master integrals in master split mode
 */
static vector<bool> needed_sectors{};
/**
 * Contains the list of non-zero masters in parser. Used in \#masters and adding rules
 */
static set<point> masters;
/**
 * Indication whether lbases have been loaded
 */
static bool lbases_loaded = false;

/**
 * In-place inversion of ordering matrix, called only during parsing.
 * @tparam M allocated dimension.
 * @param m pointer to square matrix of size M*M.
 * @param N real dimension of matrix.
 */
template<int M>
inline static void invers(double m[M][M], const int N) {
    int gaus_ostatok;
    long gaus_deter;
    double gaus_minved;
    int rn[N], cn[N];
    int j, k;

    for (j = N; j--;) {
        rn[j] = cn[j] = j;
    }
    gaus_minved = 1e99;
    gaus_deter = 1;
    for (gaus_ostatok = N; gaus_ostatok; gaus_ostatok--) {
        int jved{}, kved{};
        double vved = -1, t;

        // search of leading element
        for (j = N; j--;) {
            if (~rn[j]) {
                for (k = N; k--;) {
                    if ((~cn[k]) && (vved < fabs(m[j][k]))) {
                        vved = fabs(m[j][k]), jved = j, kved = k;
                    }
                }
            }
        }

        if (gaus_minved > vved) { gaus_minved = vved; }
        gaus_deter *= static_cast<long int>(m[jved][kved]);
        if (abs(gaus_deter) != 1) {
            cout << "error in invers: " << abs(gaus_deter) << endl;
            abort();
        }
        if (vved < 1e-99) {
            cout << "too small leading term" <<endl;
            abort();
        }

        int jt = rn[jved], kt = cn[kved];

        // rearrangement
        for (j = N; j--;)
            t = m[kt][j], m[kt][j] = m[jved][j], m[jved][j] = t;
        for (j = N; j--;)
            t = m[j][jt], m[j][jt] = m[j][kved], m[j][kved] = t;

        rn[jved] = rn[kt];
        cn[kved] = cn[jt];
        rn[kt] = cn[jt] = -1;

        vved = m[kt][jt];
        m[kt][jt] = 1;
        for (j = N; j--;) {
            if (j == kt) {
                continue;
            }
            double mul = m[j][jt] / vved;

            m[j][jt] = 0;
            for (k = N; k--;) {
                m[j][k] -= m[kt][k] * mul;
            }
        }
        for (k = N; k--;) {
            m[kt][k] /= vved;
        }
    }
}

/** Compare sector priority:
* 1) total number of positive entries
* 2) existence of an sbasis
* 3) lsectors - make sectors without Lee symmetries lower
* 4) compare by indices
* @param v1 first sector vector
* @param v2 second second vector
* @return true if first is smaller
*/
bool sector_sort_function(const vector<t_index> &v1, const vector<t_index> &v2) {
    unsigned long l = v1.size();
    int c = 0;

    for (unsigned int i = 0; i != l; ++i) {
        c = c + v1[i] - v2[i];
    }
    if (c < 0) return true;
    if (c > 0) return false;

    auto it1 = common::lsectors.find(sector_fast(v1));
    auto it2 = common::lsectors.find(sector_fast(v2));

    if ((it1 == common::lsectors.end()) && (it2 != common::lsectors.end())) return false;
    if ((it2 == common::lsectors.end()) && (it1 != common::lsectors.end())) return true;


    for (int i = l - 1; i != -1; --i) {
        c = v1[i] - v2[i];
        if (c < 0) return true;
        if (c > 0) return false;
    }
    return false;
}


int s2i(const char *digit, int &result) {
    int sign = 1;
    result = 0;
    int move = 0;
    //--- Convert each digit char and add into result.
    while ((*(digit + move) >= '0' && *(digit + move) <= '9') || *(digit + move) == '-') {
        if (*(digit + move) == '-') {
            sign = -1;
        } else {
            result = (result * 10) + (*(digit + move) - '0');
        }
        move++;
    }
    result = result * sign;
    return move;
}

int s2u(const char *digit, unsigned int &result) {
    result = 0;
    int move = 0;
    //--- Convert each digit char and add into result.
    while ((*(digit + move) >= '0' && *(digit + move) <= '9') || *(digit + move) == '-') {
        if (*(digit + move) == '-') {
            if (!move) {
                cout << "Unexpected -"<<endl;
                abort();
            } else {
                break;
            }
        } else {
            result = (result * 10) + (*(digit + move) - '0');
        }
        move++;
    }
    return move;
}

int s2v(const char *digit, vector<t_index> &result) {
    result.clear();
    result.reserve(common::dimension);
    int move = 0;
    while (*(digit + move) == ' ') { move++; }
    if (*(digit + move) != '{') {
        cout << "!error in s2v" << string(digit) << "!" << endl;
        cout << endl;
        abort();
    }
    move++;
    while (*(digit + move) == ' ') { move++; }
    while (*digit != '}') {
        int number;
        int move_now = s2i(digit + move, number);
        result.push_back(number);
        move += move_now;

        while (*(digit + move) == ' ') { move++; }
        if (*(digit + move) == '}') {
            move++;
            break;
        }
        if (*(digit + move) != ',') {
            cout << "!error in s2v" << string(digit) << "!" << endl;
            cout << endl;
            abort();
        }
        move++;
        while (*(digit + move) == ' ') { move++; }
    }
    return move;
}

int s2sf(const char *digit, SECTOR &result) {
    result = 0;
    int move = 0;
    while (*(digit + move) == ' ') { move++; }
    if (*(digit + move) != '{') {
        cout << "!error in s2v" << string(digit) << "!" << endl;
        cout << endl;
        abort();
    }
    move++;
    while (*(digit + move) == ' ') { move++; }
    while (*digit != '}') {
        int number;
        int move_now = s2i(digit + move, number);
        if (number > 0) {
            result = ((result << 1) ^ 1);
        } else {
            result = result << 1;
        }
        move += move_now;

        while (*(digit + move) == ' ') { move++; }
        if (*(digit + move) == '}') {
            move++;
            break;
        }
        if (*(digit + move) != ',') {
            cout << "!error in s2sf" << string(digit) << "!" << endl;
            cout << endl;
            abort();
        }
        move++;
        while (*(digit + move) == ' ') { move++; }
    }
    return move;
}

int s2vv(const char *digit, vector<vector<t_index> > &result) {
    result.clear();
    int move = 0;
    while (*(digit + move) == ' ') { move++; }
    if (*(digit + move) != '{') {
        cout << "error in s2vv";
        abort();
    }
    move++;
    if (*(digit + move) == '}') {
        return 2;
    }
    while (*digit != '}') {
        vector<t_index> small;
        int move_now = s2v(digit + move, small);
        result.push_back(small);
        { move += move_now; }

        while (*(digit + move) == ' ') { move++; }
        if (*(digit + move) == '}') {
            move++;
            break;
        }
        if (*(digit + move) != ',') {
            cout << "error in s2vv" << *(digit + move) << endl;
            abort();
        }
        { move++; }
        while (*(digit + move) == ' ') { move++; }
    }
    return move;
}

int s2vvv(const char *digit, vector<vector<vector<t_index> > > &result) {
    result.clear();
    int move = 0;
    while (*(digit + move) == ' ') { move++; }
    if (*digit != '{') {
        cout << "error in s2vvv" <<endl;
        abort();
    }
    move++;
    while (*digit != '}') {
        vector<vector<t_index> > small;
        int move_now = s2vv(digit + move, small);
        result.push_back(small);
        { move += move_now; }

        while (*(digit + move) == ' ') { move++; }
        if (*(digit + move) == '}') {
            move++;
            break;
        }
        if (*(digit + move) != ',') {
            cout << "error in s2vvv" << *(digit + move) << endl;
            abort();
        }
        { move++; }
        while (*(digit + move) == ' ') { move++; }
    }
    return move;
}


/**
 * @brief Generate coefficient from coefficient string in LiteRed
 * @param cc coefficient string.
 * The coefficient is free meaning it does not contain any index dependencies as most things in LiteRed do
 * Variable substitutions are performed and fermat is called in the prime version to simplify it to a number.
 * @return coefficient.
 */
COEFF generate_free_coeff_from_string(string cc) {
    for (const char &symb : cc) {
        if (symb == '[' || symb == ']') {
            cout << "strange symbol in internal symmetries coefficient\n";
            abort();
        }
    }
    cc.erase(remove(cc.begin(), cc.end(), ' '), cc.end());
    cc = ReplaceAllVariables(cc);
    cc = calc(cc);
    return COEFF(cc);
}

/**
 * Analyse internal or external symmetries to write out sector dependencies
 * @param sn of sector we are working in.
 * @param right the right-hand side of the symmetry that is a permutation and list of sums in different powers
 */
void store_symmetry_dependencies(const unsigned short sn,const pair<vector<t_index>, vector<pair<vector<pair<COEFF, point_fast> >, t_index> > >& right) {
    vector<t_index>& v = common::ssectors[sn];
    vector<t_index> vres;
    vres.reserve(common::dimension);
    for (unsigned int j = 0; j != common::dimension; j++) {
        if (right.first[j] == 0) {
            vres.push_back(-1); // no index goes there and sums can have only negative
        } else {
            vres.push_back(v[right.first[j] - 1]); // here we get the current index
        }
    }
    for (const auto& sum: right.second) {
        for (const auto& term: sum.first) {
            for (unsigned int i = 0; i!=common::dimension; ++i) {
                if ((term.second.buf[i] < 0) && (vres[i]==1)) {
                    vres[i]=0; // meaning that this index can be positive and negative
                }
            }
        }
    }
    // we have this vres, now we need to build all sectors where 0 goes to -1 or 1;
    int zeros = count(vres.begin(),vres.end(), 0);
    vector<vector<t_index> >lower_secs;
    lower_secs.reserve(1<<zeros);
    lower_secs.emplace_back(vres);
    while (zeros) {
        unsigned int pos = 0;
        while (lower_secs[0][pos]) ++pos;
        for (auto& vec: lower_secs) {
            if (!vec[pos]) { // this check is needed because we are adding new
                vector<t_index> vec2 = vec;
                vec[pos] = -1;
                vec2[pos] = 1;
                lower_secs.emplace_back(vec2);
            }
        }
        --zeros;
    }
    for (auto& vec: lower_secs) {
        unsigned short sn2 = common::sector_numbers_fast[sector_fast(vec)];
        if (sn2 && sn!=sn2) {
            dependencies[sn][sn2] = true;
            //cout<<sn<<"->"<<sn2<<endl;
        }
    }
}


/**
 * Parse lbases (Lee rules and symmetries) file.
 * @param c lbases file filename.
 * @return true if lbases were parsed successfully, false otherwise.
 */
bool add_lbases(const char *c) {
    FILE *file;
    char load_string[LOAD_STR_SIZE] = "none";
    string str;
    vector<t_index> v;
    int move;
    int n;
    unsigned int pn;
    lbases_loaded = true;

    file = fopen(c, "r");
    if (file == nullptr) {
        cerr << string(c) << endl << "Invalid lbases filename, exiting" << endl;
        return false;
    }

    while (fgets(load_string, sizeof(load_string), file)) {
        str += load_string;
    }

    for (auto & symb : str) {
        if ((symb == '\r') || (symb == '\n') || (symb == '\\')) {
            symb = ' ';
        }
    }

    const char *it = str.c_str();
    if (*it != '{') {
        cout << "wrong start (should be a pair)" << endl;
        abort();
    }
    it++;
    if (*it != '{') {
        cout << "wrong start of normal rules (should be a list)" << endl;
        abort();
    }
    it++;
    while (*it != '}') {
        if (*it != '{') {
            cout << "wrong start of rule (should be a triple)" << endl;
            abort();
        }
        it++;
        if (*it != '{') {
            cout << "wrong start of point (should be a pair)" << endl;
            abort();
        }
        it++;
        move = s2u(it, pn);
        if (pn != common::global_pn) {
            cout << "wrong problem number" << endl;
            abort();
        }
        it += move;
        if (*it != ',') {
            cout << "no comma after problem number" << endl;
            it -= 3;
            for (int i = 1; i <= 20; i++) {
                cout << *it;
                it++;
            }
            cout << endl;
            abort();
        }
        it++;
        while (*it == ' ') it++;
        if (*it != '{') {
            cout << "wrong start of sector" << endl;
            abort();
        }
        SECTOR sf;
        move = s2sf(it, sf);
        it += move;
        unsigned short sn = common::sector_numbers_fast[sf];
        if (sn == 0) {
            cout << "Undefined sector: {" << pn << ", ";
            print_sector_fast(sf);
            cout << "}" << endl;
        }
        if (*it != '}') {
            cout << "closing bracket after sector missing" << endl;
            abort();
        }
        it++;
        if (*it != ',') {
            cout << "no comma after sector" << endl;
            abort();
        }
        it++;
        while (*it == ' ') it++;
        if (*it != '{') {
            cout << "wrong start of conditions" << endl;
            abort();
        }
        it++;
        vector<pair<vector<t_index>, pair<short, bool> > > conditions;
        while (*it != '}') {
            pair<short, bool> p2; //free term and condition
            pair<vector<t_index>, pair<short, bool> > p1; //complete condition
            if (*it != '{') {
                cout << "wrong start of condition (should be a triple)" << endl;
                abort();
            }
            it++;
            if (*it != '{') {
                cout << "wrong start of indices coefficients list" << endl;
                abort();
            }
            move = s2v(it, v);
            p1.first = v;
            it += move;
            if (*it != ',') {
                cout << "no comma after indices coefficients" << endl;
                abort();
            }
            it++;
            while (*it == ' ') it++;
            move = s2i(it, n);
            p2.first = n;
            it += move;
            if (*it != ',') {
                cout << "no comma after free term" << endl;
                abort();
            }
            it++;
            while (*it == ' ') it++;
            move = s2i(it, n);
            if (n == 1) {
                p2.second = true;
            } else if (n == 0) {
                p2.second = false;
            } else {
                cout << "incorrect condition" << endl;
                abort();
            }
            p1.second = p2;
            conditions.push_back(p1);
            it += move;
            if (*it != '}') {
                cout << "no closing bracket after bool check" << endl;
                abort();
            }
            it++;
            if (*it == ',') it++;
            while (*it == ' ') it++;
        }
        it++;
        if (*it != ',') {
            cout << "comma missing after conditions" << *it << endl;
            abort();
        }
        it++;
        while (*it == ' ') it++;
        if (*it != '{') {
            cout << "wrong start of right side" << endl;
            abort();
        }
        it++;
        vector<pair<string, vector<pair<vector<t_index>, short> > > > right;
        while (*it != '}') {
            pair<string, vector<pair<vector<t_index>, short> > > term;
            while (*it == ' ') { it++; }
            if (*it != '{') {
                cout << "wrong start of a term (should be a pair)" << endl;
                abort();
            }
            it++;
            while (*it == ' ') it++;
            if (*it != '"') {
                cout << "starting quote missing for a coefficient" << endl;
                abort();
            }
            it++;
            move = 0;
            while (*(it + move) != '"') move++;
            term.first = string(it, move);
            for (char &sit : term.first) {
                if (sit == '[') sit = '(';
                if (sit == ']') sit = ')';
            }
            term.first.erase(remove(term.first.begin(), term.first.end(), ' '), term.first.end());
            term.first = ReplaceAllVariables(term.first);
            term.first = ReplaceAll(term.first, "xxx", "aaa");
            // cout << term.first << endl;
            // if you ever decide to rewrite lbases with COEFF::prime, uncomment it and have a look at the v2l example
            // the coefficients here are rational functions of multiple a[i]
            it += move;
            it++;
            while (*it == ' ') { it++; }
            if (*it != ',') {
                cout << "no comma after coefficient" << endl;
                abort();
            }
            it++;
            while (*it == ' ') { it++; }
            if (*it != '{') {
                cout << "wrong start of target indices (should be a list)" << endl;
                abort();
            }
            it++;
            vector<pair<vector<t_index>, short> > indices;

            while (*it != '}') {
                pair<vector<t_index>, short> index;
                if (*it != '{') {
                    cout << "wrong start of index definition (should be a pair)" << endl;
                    abort();
                }
                it++;
                if (*it != '{') {
                    cout << "wrong start of index row" << endl;
                    abort();
                }
                move = s2v(it, v);
                index.first = v;
                it += move;
                if (*it != ',') {
                    cout << "no comma after index row" << endl;
                    abort();
                }
                it++;
                while (*it == ' ') { it++; }
                move = s2i(it, n);
                index.second = n;
                it += move;
                while (*it == ' ') it++;
                if (*it != '}') {
                    cout << "missing closing bracket after index definition" << endl;
                    abort();
                }
                it++;
                if (*it == ',') it++;
                while (*it == ' ') it++;
                indices.push_back(index);
            }
            term.second = indices;
            it++;
            if (*it != '}') {
                cout << "missing closing bracket after a term" << endl;
                abort();
            }
            it++;
            if (*it == ',') it++;
            while (*it == ' ') it++;
            right.push_back(term);
        }
        it++;
        if (*it != '}') {
            cout << "missing closing bracket after a rule" << endl;
            abort();
        }
        it++;
        if (*it == ',') it++;
        while (*it == ' ') it++;
        if (sn != 0) {   // loading it only for our sector
            auto it0 = common::lbases.find(sn);
            vector<pair<vector<pair<vector<t_index>, pair<short, bool> > >, vector<pair<string, vector<pair<vector<t_index>, short> > > > > > lbasis;
            if (it0 != common::lbases.end()) {
                lbasis = it0->second;
            }
            lbasis.emplace_back(conditions, right);
            if (it0 != common::lbases.end()) {
                it0->second = lbasis;
            } else {
                common::lbases.emplace(sn, lbasis);
            }
        }
    }

    it++;
    while (*it == ' ') { it++; }
    if (*it != ',') {
        cout << "missing comma after normal rules" << endl;
        abort();
    }
    it++;

    // here we started the second part
    while (*it == ' ') { it++; }
    if (*it != '{') {
        cout << "wrong start of delayed rules (should be a list)" << endl;
        abort();
    }
    it++;
    while (*it != '}') {
        while (*it == ' ') { it++; }
        if (*it != '{') {
            cout << "wrong start of a delayed term (should be a triple)" << endl;
            abort();
        }
        it++;
        while (*it == ' ') { it++; }

        unsigned int here_pn;
        if (*it != '{') {
            cout << "wrong start of a point in delayed rule" << *it << endl;
            abort();
        }
        it++;
        move = s2u(it, here_pn);
        if (here_pn != common::global_pn) {
            cout << "wrong problem number" << endl;
            abort();
        }
        it += move;
        while (*it == ' ') it++;
        if (*it != ',') {
            cout << "missing comma after a problem number in delayed rule" << endl;
            abort();
        }
        it++;
        while (*it == ' ') it++;
        if (*it != '{') {
            cout << "wrong start of a sector" << *it << endl;
            abort();
        }
        SECTOR sf;
        move = s2sf(it, sf);
        it += move;
        while (*it == ' ') it++;
        if (*it != '}') {
            cout << "missing closing bracket after a sector" << endl;
            abort();
        }
        it++;
        while (*it == ' ') it++;
        if (*it != ',') {
            cout << "missing comma after a point" << endl;
            abort();
        }
        it++;
        unsigned short sn = common::sector_numbers_fast[sf];
        if ((sn == 0) || (sn == 1) || (sn == common::virtual_sector))  {  // printing only in main thread
            cout << "Ignoring rule for sector: {" << common::global_pn << ", ";
            print_sector_fast(sf);
            cout << "}, reason: ";
            if (sn == 0) { cout << "zero sector"; }
            if (sn == 1) { cout << "sector does not exist"; }
            if (sn == common::virtual_sector) { cout << "global symmetry used"; }
            cout << endl;
        }
        while (*it == ' ') it++;
        vector<t_index> permutation;
        if (*it != '{') {
            cout << "wrong start of a permutation" << *it << endl;
            abort();
        }
        move = s2v(it, permutation);
        it += move;
        while (*it == ' ') it++;
        if (*it != ',') {
            cout << "missing comma after a permutation" << endl;
            abort();
        }
        it++;
        while (*it == ' ') it++;
        if (*it != '{') {
            cout << "wrong start of product (should be a list)" << *it << endl;
            abort();
        }
        it++;
        pair<vector<t_index>, vector<pair<vector<pair<COEFF, point_fast> >, t_index> > > right;
        right.first = permutation;
        while (*it != '}') {  //product
            pair<vector<pair<COEFF, point_fast> >, t_index> pterm;
            while (*it == ' ') it++;
            if (*it != '{') {
                cout << "wrong start of a product term (should be a pair)" << endl;
                abort();
            }
            it++;
            while (*it == ' ') it++;
            if (*it != '{') {
                cout << "wrong start of a sum (should be a list)" << endl;
                abort();
            }
            it++;
            while (*it != '}') {  //sum
                pair<COEFF, point_fast> term;
                while (*it == ' ') it++;
                if (*it != '{') {
                    cout << "wrong start of a sum term (should be a pair)" << endl;
                    abort();
                }
                it++;
                while (*it == ' ') it++;
                if (*it != '"') {
                    cout << "starting quote missing for a coefficient" << *it << endl;
                    abort();
                }
                it++;
                move = 0;
                while (*(it + move) != '"') move++;

                string cc = string(it, move);
                term.first = generate_free_coeff_from_string(cc);

                it += move;
                it++;
                while (*it == ' ') it++;
                if (*it != ',') {
                    cout << "no comma after coefficient" << endl;
                    abort();
                }
                it++;
                while (*it == ' ') it++;
                if (*it != '{') {
                    cout << "wrong start of a sum term indices (should be a list)" << endl;
                    abort();
                }
                vector<t_index> tempv;
                move = s2v(it, tempv);
                term.second = point_fast(tempv);

                it += move;
                while (*it == ' ') it++;
                if (*it != '}') {
                    cout << "wrong end of a sum term (should be a pair)" << endl;
                    abort();
                }
                it++;
                while (*it == ' ') it++;
                if (*it == ',') it++;
                while (*it == ' ') it++;
                pterm.first.push_back(term);
            }
            it++;
            while (*it == ' ') it++;
            if (*it != ',') {
                cout << "comma missing after sum" << endl;
                abort();
            }
            it++;
            while (*it == ' ') it++;
            move = s2i(it, n);
            pterm.second = n;
            it += move;
            while (*it == ' ') it++;
            if (*it != '}') {
                cout << "missing closing bracket after product term" << endl;
                abort();
            }
            it++;
            while (*it == ' ') it++;
            if (*it == ',') it++;
            while (*it == ' ') it++;
            right.second.push_back(pterm);
        }
        if ((sn != 0) && (sn != 1) && (sn != common::virtual_sector)) {  // loading only our sector
            auto it0 = point::dbases.find(sn);
            if (it0 != point::dbases.end()) {
                cout << "Repeated delayed rule!!!";
                abort();
            }
            point::dbases.emplace(sn, right);
        }
        if ((sn != 0) && (sn != 1) && (sn != common::virtual_sector)) {
            if (common::split_masters) {
                store_symmetry_dependencies(sn,right);
            }
        }
        it++;
        while (*it == ' ') it++;
        if (*it != '}') {
            cout << "missing closing bracket after delayed rule" << endl;
            abort();
        }
        it++;
        while (*it == ' ') it++;
        if (*it == ',') it++;
        while (*it == ' ') it++;
    }
    it++;

    if (*it != '}') { // old files had only 2 members
        int internal_count = 0;
        if (*it != ',') {
            cout << "wrong end" << endl;
            abort();
        }
        it++;
        while (*it == ' ') it++;
        if (*it != '{') {
            cout << "wrong start of symmetry rules (should be a list)" << endl;
            abort();
        }
        it++;
        while (*it != '}') {
            while (*it == ' ') it++;
            if (*it != '{') {
                cout << "wrong start of a symmetry term (should be a triple)" << endl;
                abort();
            }
            it++;
            while (*it == ' ') it++;
            unsigned int here_pn;
            if (*it != '{') {
                cout << "wrong start of a point in symmetry rule" << *it << endl;
                abort();
            }
            it++;
            move = s2u(it, here_pn);
            if (here_pn != common::global_pn) {
                cout << "wrong problem number" << endl;
                abort();
            }
            it += move;
            while (*it == ' ') it++;
            if (*it != ',') {
                cout << "missing comma after a problem number in symmetry rule" << endl;
                abort();
            }
            it++;
            while (*it == ' ') it++;
            if (*it != '{') {
                cout << "wrong start of a sector" << *it << endl;
                abort();
            }
            SECTOR sf;
            move = s2sf(it, sf);
            it += move;
            while (*it == ' ') it++;
            if (*it != '}') {
                cout << "missing closing bracket after a sector" << endl;
                abort();
            }
            it++;
            while (*it == ' ') it++;
            if (*it != ',') {
                cout << "missing comma after a point" << endl;
                abort();
            }
            it++;
            unsigned short sn = common::sector_numbers_fast[sf];
            if ((sn == 0) || (sn == 1) || (sn == common::virtual_sector)) {  // printing only in main thread
                cout << "Ignoring rule for sector: {" << common::global_pn << ", ";
                print_sector_fast(sf);
                cout << "}, reason: ";
                if (sn == 0) cout << "zero sector";
                if (sn == 1) cout << "sector does not exist";
                if (sn == common::virtual_sector) cout << "global symmetry used";
                cout << endl;
            }
            while (*it == ' ') it++;
            vector<t_index> permutation;
            if (*it != '{') {
                cout << "wrong start of a permutation" << *it << endl;
                abort();
            }
            move = s2v(it, permutation);
            it += move;
            while (*it == ' ') it++;
            if (*it != ',') {
                cout << "missing comma after a permutation" << endl;
                abort();
            }
            it++;
            while (*it == ' ') it++;
            if (*it != '{') {
                cout << "wrong start of product (should be a list)" << *it << endl;
                abort();
            }
            it++;
            pair<vector<t_index>, vector<pair<vector<pair<COEFF, point_fast> >, t_index> > > right;
            right.first = permutation;
            while (*it != '}') { //product
                pair<vector<pair<COEFF, point_fast> >, t_index> pterm;
                while (*it == ' ') it++;
                if (*it != '{') {
                    cout << "wrong start of a product term (should be a pair)" << endl;
                    abort();
                }
                it++;
                while (*it == ' ') it++;
                if (*it != '{') {
                    cout << "wrong start of a sum (should be a list)" << endl;
                    abort();
                }
                it++;
                while (*it != '}') {  //sum
                    pair<COEFF, point_fast> term;
                    while (*it == ' ') it++;
                    if (*it != '{') {
                        cout << "wrong start of a sum term (should be a pair)" << endl;
                        abort();
                    }
                    it++;
                    while (*it == ' ') it++;
                    if (*it != '"') {
                        cout << "starting quote missing for a coefficient" << *it << endl;
                        abort();
                    }
                    it++;
                    move = 0;
                    while (*(it + move) != '"') move++;

                    string cc = string(it, move);
                    term.first = generate_free_coeff_from_string(cc);

                    it += move;
                    it++;
                    while (*it == ' ') it++;
                    if (*it != ',') {
                        cout << "no comma after coefficient" << endl;
                        abort();
                    }
                    it++;
                    while (*it == ' ') it++;
                    if (*it != '{') {
                        cout << "wrong start of a sum term indices (should be a list)" << endl;
                        abort();
                    }
                    vector<t_index> tempv;
                    move = s2v(it, tempv);
                    term.second = point_fast(tempv);
                    it += move;
                    while (*it == ' ') it++;
                    if (*it != '}') {
                        cout << "wrong end of a sum term (should be a pair)" << endl;
                        abort();
                    }
                    it++;
                    while (*it == ' ') it++;
                    if (*it == ',') it++;
                    while (*it == ' ') it++;
                    pterm.first.push_back(term);
                }
                it++;
                while (*it == ' ') it++;
                if (*it != ',') {
                    cout << "comma missing after sum" << endl;
                    abort();
                }
                it++;
                while (*it == ' ') it++;
                move = s2i(it, n);
                pterm.second = n;
                it += move;
                while (*it == ' ') it++;
                if (*it != '}') {
                    cout << "missing closing bracket after product term" << endl;
                    abort();
                }
                it++;
                while (*it == ' ') it++;
                if (*it == ',') it++;
                while (*it == ' ') it++;
                right.second.push_back(pterm);
            }

            // loading only our sector
            if ((sn != 0) && (sn != 1) && (sn != common::virtual_sector)) {
                auto it0 = point::ibases.find(sn);
                vector<pair<vector<t_index>, vector<pair<vector<pair<COEFF, point_fast> >, t_index> > > > irules;
                if (it0 != point::ibases.end()) {
                    irules = it0->second;
                }
                irules.push_back(right);
                ++internal_count;
                if (it0 != point::ibases.end()) {
                    it0->second = irules;
                } else {
                    point::ibases.emplace(sn, irules);
                }
            }

            if ((sn != 0) && (sn != 1) && (sn != common::virtual_sector)) {
                if (common::split_masters) {
                    store_symmetry_dependencies(sn,right);
                }
            }

            it++;
            while (*it == ' ') it++;
            if (*it != '}') {
                cout << "missing closing bracket after delayed rule" << endl;
                abort();
            }
            it++;
            while (*it == ' ') it++;
            if (*it == ',') it++;
            while (*it == ' ') it++;
        }
        it++;
    }

    if (*it != '}') {
        cout << "missing final closing bracket" << endl;
        abort();
    }

    fclose(file);
    return true;
}

/**
 * Add a point to list of preferred points in point's sector.
 * @param v vector of point's simple coefficients.
 * @return True if it is a new integral
 */
bool add_single_preferred(const vector<t_index> &v) {
    SECTOR s = sector_fast(v);
    unsigned short sn = common::sector_numbers_fast[s];
    bool result = (point::preferred[sn].find(v)==point::preferred[sn].end());
    point::preferred[sn].insert(v);
    point::preferred_fast[sn].insert(point_fast(v));
    return result;
}

/**
 * Read the list of prefered integrals (or masters in case of split master mode) and write them to memory
 * @param c preferred file filename.
 * @return the number of integrals in the file
 */
int load_preferred(const char* c) {
    FILE *preferred_int_file = fopen(c,"r");
    char load_string[LOAD_STR_SIZE] = "none";
    string contents;
    if (preferred_int_file == nullptr) {
        cout << "File with preferred integrals could not be opened, exiting" << endl;
        return -1;
    }
    while (fgets(load_string, sizeof(load_string), preferred_int_file)) {
        contents += load_string;
    }
    contents = contents.substr(contents.find('{') + 1); //"}"
    for (unsigned int i = 0; i != contents.size(); ++i) {
        if (contents[i] == '\n') contents[i] = ' ';
    }
    for (unsigned int i = 0; i != contents.size(); ++i) {
        if (contents[i] == '\r') contents[i] = ' ';
    }
    const char *all = contents.c_str();
    int move = 0;
    int pref_count = 0;
    while (all[move] != '}') {
        while (all[move] == ' ') move++;
        if (all[move] != '{') {
            cout << "Error in pref:" << all[move] << endl;
            abort();
        }
        move++;
        while (all[move] == ' ') move++;
        unsigned int i;
        move += s2u(all + move, i);
        if (i != common::global_pn) {
            cerr << "Unknown problem number " << i << endl;
            fclose(preferred_int_file);
            return -1;
        }
        while (all[move] == ' ') move++;
        if (all[move] != ',') {
            cout << "Error in pref:";
            abort();
        }
        move++;
        while (all[move] == ' ') move++;
        vector<t_index> v;
        move += s2v(all + move, v);

        // the set of manually set preferred
        bool new_pref = add_single_preferred(v);
        if (common::split_masters && new_pref) {
            cout << "WARNING: PREFERED INTEGRAL THAT IS NOT IN THE LIST OF MASTERS: ";
            print_vector(v);
            cout<<endl;
        }
        ++pref_count;

        while (all[move] == ' ') move++;
        if (all[move] != '}') {
            cout << "Error in pref:";
            abort();
        }
        move++;
        while (all[move] == ' ') move++;
        if (all[move] == '}') break;
        if (all[move] != ',') {
            cout << "Error in pref:";
            abort();
        }
        move++;
    }
    fclose(preferred_int_file);
    return pref_count;
}

/**
 * Read the list of masters and set most of them to zero
 * @param c masters file filename.
 * @param zero_sector indication whether it is the main process and we should write masters to databases
 */
void set_masters_to_zero(const char *c, bool zero_sector =false) {
    FILE *file;
    char load_string[LOAD_STR_SIZE] = "none";
    string str;
    string str_p;
    //size_t found;
    unsigned int master_counter = 1;

    file = fopen(c, "r");
    if (file == nullptr) {
        cout << string(c) << endl << "Invalid masters filename, exiting" << endl;
        abort();
    }

    string contents{};
    while (fgets(load_string, sizeof(load_string), file)) {
        contents += load_string;
    }
    contents = contents.substr(contents.find('{') + 1); //"}"
    for (unsigned int i = 0; i != contents.size(); ++i) {
        if (contents[i] == '\n') contents[i] = ' ';
    }
    for (unsigned int i = 0; i != contents.size(); ++i) {
        if (contents[i] == '\r') contents[i] = ' ';
    }
    const char *all = contents.c_str();
    int move = 0;
    while (all[move] != '}') {
        while (all[move] == ' ') move++;
        if (all[move] != '{') {
            cout << "Error in masters: " << all[move] << endl;
            abort();
        }
        move++;
        while (all[move] == ' ') move++;
        unsigned int i;
        move += s2u(all + move, i);
        if (i != common::global_pn) {
            cout << "Unknown problem number in master" << i << endl;
            abort();
        }

        while (all[move] == ' ') move++;
        if (all[move] != ',') {
            cout << "Error in integrals, missing comma";
            abort();
        }
        move++;
        while (all[move] == ' ') move++;
        vector<t_index> v;
        move += s2v(all + move, v);

        point p = point_reference(v);
        if (!p.is_zero()) {
            add_single_preferred(v);
            pc_pair_ptr_vec terms;

            p = point_reference(v); //once again, now preferred
            if ((master_counter >= common::master_number_min) && (master_counter <= common::master_number_max)) {
                needed_sectors[p.s_number()] = true;
                masters.insert(p);
                point p1(p.get_vector(), 0, -2); // sending to sector 1
                terms.emplace_back(make_pc_ptr(p1,CO_1));
            } else {
                if (!common::silent && zero_sector) cout<<"Setting master integral to zero: "<<p<<endl;
            }

            if (zero_sector) {
                terms.emplace_back(make_pc_ptr(p, CO_1m));
                
                if(!(common::run_mode>1 && read_status(p.s_number())!=0) && !common::run_sector) {
                    open_database(p.s_number());
                    p_set(p, std::move(terms), 127);
                    close_database(p.s_number());
                }
            }
            ++master_counter;
        }

        while (all[move] == ' ') move++;
        if (all[move] != '}') {
            cout << "Error in masters, missing closing bracket ";
            abort();
        }
        move++;
        while (all[move] == ' ') move++;
        if (all[move] == '}') break;
        if (all[move] != ',') {
            cout << "Error in masters, missing comma between integrals";
            abort();
        }
        move++;
    }
    if (master_counter <= common::master_number_max) {
        cout << "Specified maximal master number ("<<common::master_number_max<<") is bigger than the number of masters in file ("<<master_counter-1<<")"<<endl;
        abort();
    }
    fclose(file);
}


/**
 * Parse rules (master integral mappings) file.
 * @param c rules file filename.
 * @return true if lbases were parsed successfully, false otherwise.
 */
bool add_rules(const char *c) {
    FILE *file;
    char load_string[LOAD_STR_SIZE] = "none";
    string str;
    string str_p;
    size_t found;

    file = fopen(c, "r");
    if (file == nullptr) {
        cerr << string(c) << endl << "Invalid rules filename, exiting" << endl;
        return false;
    }

    int rule_counter = 0;
    while (fgets(load_string, sizeof(load_string), file)) {
        for (int j = 0; j != sizeof(load_string); ++j) {
            if (load_string[j] == '\r') { load_string[j] = '\n'; }
        }

        if (load_string[1] == '\n' || load_string[0] == '\n') {
            for (auto & symb : str) {
                if (symb == '\n') symb = ' ';
            }

            int n, n2;
            rule_counter++;
            found = str.find('G');
            str = str.substr(found + 2);
            int move = s2i(str.c_str(), n);

            str = str.substr(move + 2);
            vector<t_index> v;
            s2v(str.c_str(), v);
            found = str.find("->");
            str = str.substr(found + 2);

            // left-hand side of the rule
            add_single_preferred(v);

            point p = point_reference(v);
            if ((p.is_zero())) {
                cout << "Attempt to create a rule for {" << n << ", ";
                print_vector(v);
                cout << "} ignored" << endl;
                rule_counter--;
                continue;
            }

            if ((p.s_number() == 1)) {
                cout << "Attempt to create a rule for {" << n << ", ";
                print_vector(v);
                cout << "} ignored" << endl;
                rule_counter--;
                continue;
            }

            pc_pair_vec mon;
            if (str.find('{') < str.find('G')) { //"}"
                str = str.substr(str.find('{') + 1); //"}"
            }
            while (str.find('G') != string::npos) {   // here we build it just with strings
                while ((str[0] == ' ') || (str[0] == '{') || (str[0] == '}') || (str[0] == ',')) {
                    str = str.substr(1);
                }
                found = str.find('G');
                string current_coeff = found ? str.substr(0, found) : "1";
                if (current_coeff.find(',') != string::npos) {
                    current_coeff = current_coeff.substr(0, current_coeff.find(','));
                }
                current_coeff = ReplaceAllVariables(current_coeff);
                
                current_coeff.erase(std::remove(current_coeff.begin(), current_coeff.end(), ' '), current_coeff.end());
                current_coeff = calc(current_coeff);
                
                #if defined(PRIME)
                char buf[COEFF_BUF_SIZE];
                sprintf(buf, "%s|%lu", current_coeff.c_str(), static_cast<unsigned long int>(COEFF::prime)); // or it can become big!
                current_coeff = calc(buf);
                #endif
                
                COEFF current_COEFF(current_coeff);
                str = str.substr(found + 2);
                move = s2i(str.c_str(), n2);
                str = str.substr(move + 2);
                vector<t_index> v2;
                move = s2v(str.c_str(), v2);
                str = str.substr(move + 2);

                // right-hand side
                add_single_preferred(v2);
                point p2 = point_reference(v2);

                if (common::split_masters) {
                    dependencies[p.s_number()][p2.s_number()] = true;
                    //if (!sector) cout<<p.s_number()<<"->"<<p2.s_number()<<endl;
                }

                // reading is done, preferred is set, that's all we need for other sectors
                {

                    point p3(v2, 0, -2); // sending to sector 1

                    if (masters.find(p2) == masters.end()) {
                        // it's not listed in masters
                        // either it is split mode, and we have to check if it is zero
                        // or it is normal mode and we need to write a relation
                        if(!(common::run_mode>1 && read_status(p2.s_number())!=0) && !common::run_sector) {
                            open_database(p2.s_number());

                            if (common::split_masters) {
                                if (p_get_monoms(p2,p2.s_number()).size() != 1) {
                                    cout << "Master integral in right-hand side of rules not listed in masters file"<<endl;
                                    abort();
                                }
                            } else {
                                pc_pair_ptr_vec t;
                                t.emplace_back(make_pc_ptr(p3, CO_1));
                                t.emplace_back(make_pc_ptr(p2, CO_1m));
                                p_set(p2, std::move(t), 127);  // the rule right-hand sides are masters
                            }
                            close_database(p2.s_number());
                        }
                        masters.insert(p2);
                        mon.emplace_back(p3, std::move(current_COEFF));
                    } else {
                        mon.emplace_back(p3, std::move(current_COEFF));
                    }
                }
            }

            {
                // no need to write entry for the rule for other sectors
                sort(mon.begin(), mon.end(), pair_point_COEFF_smaller);

                pc_pair_ptr_vec terms;

                for (auto read = mon.begin(); read != mon.end(); ++read) {
                    COEFF coeff2;
                    _set_(coeff2, read->second);
                    auto read2 = read;
                    read2++;
                    while ((read2 != mon.end()) && (read2->first == read->first)) {
                        _add_(coeff2, coeff2, read2->second);
                        read2++;
                    }
                    read2--;
                    read = read2;

                    if (!coeff2.is_zero()) {
                        terms.emplace_back(make_pc_ptr(read->first, std::move(coeff2)));
                    }
                }

                terms.emplace_back(make_pc_ptr(p, CO_1m));

                if(!(common::run_mode>1 && read_status(p.s_number())!=0) && !common::run_sector) {
                    open_database(p.s_number());
                    p_set(p, std::move(terms), 126);  // the left rule side is not really needed... but should not hurt
                    close_database(p.s_number());
                }
            }
            str = "";
        } else {
            str += load_string;
            if (str == string("Null\n")) {
                str = "";
            }
        }
    }
    common::prt_rule_counter = rule_counter;
    fclose(file);
    return true;
}

/**
 * Add corresponding problem to global data, parsing of start file or sbases file.
 * @param problem_number number of problem.
 * @param cc problem specification.
 * @return true if problem was added successfully, false otherwise.
 */
bool add_problem(const unsigned int problem_number, const char *cc) {
    common::global_pn = problem_number;
    int positive = 0;
    int positive_start = 1;
    int move = 0;
    if (*cc == '|') {
        move = s2i(cc + 1, positive) + 1;
        if (*(cc + move) == ',') {
            positive_start = positive;
            positive = 0;
            move += (s2i(cc + move + 1, positive) + 2);
        } else if (*(cc + move) == '|') {
            move++;
        } else {
            cerr << "Incorrect problem syntax" << endl;
            abort();
        }
    }
    const char *c;
    string tmp;
    if (*cc == '/') {
        c = cc;
    } else {
        tmp = folder + string(cc + move);
        c = tmp.c_str();
    }

    FILE *file;
    char load_string[LOAD_STR_SIZE] = "none";
    string str;
    string str_p;
    size_t found;

    int n, m;
    file = fopen(c, "r");
    if (file == nullptr) {
        cerr << string(c) << endl << "Invalid problem filename, exiting" << endl;
        return false;
    }

    str = "";

    vector<vector<vector<t_index> > > local_symmetries;

    unsigned int dimension = 0;
    while (fgets(load_string, sizeof(load_string), file)) {
        if (load_string[1] == '\n' || load_string[0] == '\n' || load_string[0] == '\r' || load_string[1] == '\r') {
            for (auto & symb : str) {
                if ((symb == '\n') || (symb == '\r')) {
                    symb = ' ';
                }
            }
            string substr = ("ExampleDimension");
            if (str.substr(0,substr.size()) == substr) {
                found = str.find('=');
                str = str.substr(found + 2);
                s2i(str.c_str(), n);
                dimension = n;
                if (positive == 0) { positive = dimension; }

                vector<vector<t_index> > all = all_sectors(n, positive, positive_start);

                if (dimension > MAX_IND) {
#ifdef SMALL_POINT
                    cout << "Maximal dimension with --small_point is " << MAX_IND << endl;
#else
                    cout << "Maximal dimension is " << MAX_IND << endl;
#endif
                    abort();
                }

                common::dimension = n;
                common::orderings_fast = unique_ptr<unique_ptr<t_index[]>[]>{new unique_ptr<t_index[]>[1 << (n)]{}};

                size_t max_sectors = (1u << n);
                common::sector_numbers_fast = unique_ptr<unsigned short[]>{new unsigned short[max_sectors]()};

                point::preferred.reserve(max_sectors+1);
                point::preferred_fast.reserve(max_sectors+1);
                for (unsigned int i = 0; i != max_sectors+1; ++i) {
                    point::preferred.emplace_back();
                    point::preferred_fast.emplace_back();
                }

                // for fast getting sector numbers
                // we fill it with -1 now and set 0 for non-existing sectors
                for (uint32_t i = 0; i != max_sectors; ++i) {
                    if (i >> (n - positive_start +1))  {
                        // highest bits correspond to first indices
                        common::sector_numbers_fast[i] = 0;
                    }
                    else if (i & ((1 << (n - positive)) - 1))  {
                        common::sector_numbers_fast[i] = 0;
                    }
                    else common::sector_numbers_fast[i] = static_cast<unsigned short>(-1); // number will be given later
                }
                str = "";
                continue;
            }

            substr = ("SBasis0L");
            if (str.substr(0,substr.size()) == substr) {
                if (true) {
                    found = str.find('=');
                    str = str.substr(found + 2);
                    s2i(str.c_str(), n);
                    for (int i = 0; i < n; ++i) {
                        point::ibps.emplace_back();
                    }
                }
                str = "";
                continue;
            }

            substr = ("SBasisR");
            if (str.substr(0,substr.size()) == substr) {
                found = str.find("True");
                if (found != string::npos) {
                    found = str.find(',');
                    str_p = str.substr(found + 2);
                    SECTOR sf;
                    s2sf(str_p.c_str(), sf);

                    common::sector_numbers_fast[sf] = 0; // needed, unnumbered are -1
                }
                str = "";
                continue;
            }

            substr = ("SBasis0C");
            if (str.substr(0,substr.size()) == substr) {
                if (true) {
                    found = str.find(',');
                    str = str.substr(found + 2);
                    s2i(str.c_str(), m);

                    found = str.find(',');
                    str = str.substr(found + 3);
                    found = str.find(']');
                    str_p = str.substr(0, found - 1);
                    str = str.substr(found + 4);
                    vector<t_index> v;
                    while (found != string::npos) {
                        s2i(str_p.c_str(), n);
                        v.push_back(n);
                        found = str_p.find(',');
                        str_p = str_p.substr(found + 2);
                        while (str_p[0] == ' ') str_p = str_p.substr(1, str_p.size() - 1);
                    }

                    str = ReplaceAllVariables(str);

                    // now it should be just the whole coefficient, that should be better parsed here

                    // due to substitutions some of the coefficients might become 0. we have to get them away
                    found = str.find('{');
                    string str_res = str.substr(0, found + 1);
                    str = str.substr(found + 1); // to the global opening bracket
                    bool good_term = false;
                    while (true) {
                        found = str.find('{');
                        if (found == string::npos) break;
                        str_res += str.substr(0, found + 1);
                        str = str.substr(found + 1);
                        found = str.find(',');
                        string string_temp = str.substr(0, found);
                        string_temp = calc(string_temp);
                        if (string_temp != "0") {
                            str_res += string_temp;
                        }
                        str = str.substr(found);
                        found = str.find('}');
                        if (string_temp != "0") {
                            str_res += str.substr(0, found + 1);
                            good_term = true;
                        } else {
                            str_res = str_res.substr(str_res.length() - 1);
                            if (str[found + 1] == ',') { ++found; } // this term completely ignored
                        }
                        str = str.substr(found + 1);
                    }
                    str_res += str;
                    str = str_res;
                    if (good_term) {
                        point::ibps[m - 1].emplace_back(split_coeff(str),point_fast(v));
                    }
                }
                str = "";
                continue;
            }

            substr = ("SBasisS");
            if (str.substr(0,substr.size()) == substr) {
                found = str.find('{'); //"}"
                str = str.substr(found);
                s2vvv(str.c_str(), local_symmetries);
                str = "";
                continue;
            }

            substr = ("SBasisN");
            if (str.substr(0,substr.size()) == substr) {
                found = str.find(',');
                str = str.substr(found + 1);
                const char *pos = str.c_str();
                SECTOR sf;
                s2sf(pos, sf);
                if (common::sector_numbers_fast[sf]) {
                    common::lsectors.emplace(sf); // only if this is a non-zero sector (actually -1 at this point);
                }
                str = "";
                continue;
            }

            substr = ("SBasisO");
            if (str.substr(0,substr.size()) == substr) {
                found = str.find(',');
                str = str.substr(found + 2);
                vector<t_index> sec;

                int new_move = s2v(str.c_str(), sec);

                found = str.find('{', new_move); //"}"
                str = str.substr(found);

                vector<vector<t_index> > ord;
                s2vv(str.c_str(), ord);

                common::orderings_fast[sector_fast(sec)] = unique_ptr<t_index[]>(new t_index[common::dimension * common::dimension]);
                auto mat = common::orderings_fast[sector_fast(sec)].get();
                int row = 0;
                for (auto itr_row = ord.begin(); itr_row != ord.end(); ++itr_row, ++row) {
                    // now we stopped using vectors, moving to arrays, it's much faster
                    // and array of vectors leads to crashes
                    int column = 0;
                    for (auto itr_column = itr_row->begin(); itr_column != itr_row->end(); ++itr_column, ++column) {
                        mat[row * common::dimension + column] = *itr_column;
                    }
                }
                //}
                str = "";
                continue;
            }
            str = "";
        } else {
            str += load_string;
            if (str == string("Null\n")) { str = string(""); }
        }
    }

    vector<t_index> uuuu;
    uuuu.push_back(positive);
    for (auto &local_symmetry : local_symmetries) {
        local_symmetry.push_back(uuuu);
    }

    common::symmetries = local_symmetries;
    if (dimension == 0) {
        cout << "Something weird in add_problem - dimension wasn't set!" << endl;
        abort();
    }
    vector<vector<t_index> > all0 = all_sectors(dimension, positive, positive_start);
    vector<vector<t_index> > all;
    all.reserve(all0.size());

    for (const auto & sector_in_all : all0) {
        if (common::sector_numbers_fast[sector_fast(sector_in_all)] == static_cast<unsigned short>(-1)) {
            all.push_back(sector_in_all);
        }
        // non-zero sectors
        // we will enumerate later, so not changing sector_numbers_fast
    }

    std::sort(all.begin(), all.end(), sector_sort_function);

    common::ssectors.emplace_back();
    common::ssectors.emplace_back();
    unsigned int current_sector = 2;

    for (const auto &l_sector : all) {
        // we do not enumerate sectors with level 16+
        if (positive_index(l_sector) <= 15) {
            vector<vector<t_index> > orbit;
            Orbit(l_sector, orbit, local_symmetries);
            vector<t_index> *lowest = &(*orbit.begin());
            auto itr = orbit.begin();
            itr++;
            while (itr != orbit.end()) {
                if (sector_sort_function(*itr, *lowest)) lowest = &(*itr);
                itr++;
            }

            if (common::sector_numbers_fast[sector_fast(*lowest)] == static_cast<unsigned short>(-1)) {
                common::sector_numbers_fast[sector_fast(*lowest)] = current_sector;
                common::ssectors.push_back(*lowest);
                if (positive_index(*lowest) < common::abs_min_level) { common::abs_min_level = positive_index(*lowest); }

                ++current_sector;
                if (current_sector == MAX_SECTORS - 2) {
                    cerr << "Too many non-zero sectors!";
                    fclose(file);
                    return false;
                }
            }
        }
    };

    common::virtual_sector = current_sector;
    for (const auto &l_sector : all) {
        // we do not enumerate sectors with level 16+
        if (positive_index(l_sector) <= 15) {
            if (common::sector_numbers_fast[sector_fast(l_sector)] == static_cast<unsigned short>(-1)) {
                common::sector_numbers_fast[sector_fast(l_sector)] = common::virtual_sector;
                // all sectors in higher symmetry orbits are marked as common::virtual_sector for faster comparing
            }
        }
    }

    //virtual
    vector<t_index> *v;
    if (all.empty()) {
        v = &(*all0.rbegin());
    } else {
        v = &(*all.rbegin());
    }
    common::ssectors.push_back(*v);

    common::abs_max_sector = current_sector - 1;
    if (common::abs_min_level == 0) {
        cout << "Non-zero zero sector (no restrictions at all)" << endl;
        abort();
    }


    // invert orderings
    int sn = 0;
    common::iorderings.resize(common::abs_max_sector + 2);

    for (const auto &ssector : common::ssectors) {
        if (sn > 1) {
            if (common::orderings_fast[sector_fast(ssector)] == nullptr) {
                common::orderings_fast[sector_fast(ssector)] = unique_ptr<t_index[]>(new t_index[common::dimension * common::dimension]);

                make_ordering(common::orderings_fast[sector_fast(ssector)].get(), ssector);
            }
            if (dimension == 0) {
                cout << "Something weird in add_problem - dimension wasn't set!" << endl;
                abort();
            }
            double matrix[MAX_IND][MAX_IND];
            for (unsigned int i = 0; i != dimension; ++i) {
                for (unsigned int j = 0; j != dimension; ++j) {
                    matrix[i][j] = (common::orderings_fast[sector_fast(ssector)][i * common::dimension + j]);
                }
            }

            invers(matrix, dimension);
            vector<vector<t_index> >& iord = common::iorderings[sn];
            for (unsigned int i = 0; i != dimension; ++i) {
                vector<t_index> temp;
                for (unsigned int j = 0; j != dimension; ++j) {
                    if (double(t_index(matrix[i][j])) != matrix[i][j]) {
                        cout << "Bad inverse!" << endl;
                        abort();
                    }
                    temp.push_back(t_index(matrix[i][j]));
                }
                iord.push_back(temp);
            }
        }
        ++sn;
    }
    fclose(file);
    return true;
}

/**
 * @brief Initialies dependency relations by sectors
 * Knowing the sector numbers creates proper matrices
 * For normal sectors writes out the IBP go down relations
 */
void initialize_dependencies() {
    dependencies.reserve(common::abs_max_sector + 1);
    needed_sectors.reserve(common::abs_max_sector + 1);
    for (int i = 0; i!=common::abs_max_sector + 1; ++i) {
        vector<bool> dependency;
        if (i>1) {
            dependency.reserve(common::abs_max_sector + 1);
            for (int j = 0; j!=common::abs_max_sector + 1; ++j) {
                dependency.push_back(false);
            }
            if (in_lsectors(i)) {
                // normal sector, lets set down dependency right now
                vector<t_index> v = common::ssectors[i];
                for (int k=0; k!=common::dimension;++k) {
                    if (v[k]==1) {
                        auto v2(v);
                        v2[k] = -1;
                        auto num = common::sector_numbers_fast[sector_fast(v2)];
                        if (num == static_cast<unsigned short>(-1)) {
                            // that's a non lowest sector in an orbit
                            vector<vector<t_index> > orbit;
                            Orbit(v2, orbit, common::symmetries);
                            vector<t_index> *lowest = &(*orbit.begin());
                            auto itr = orbit.begin();
                            itr++;
                            while (itr != orbit.end()) {
                                if (sector_sort_function(*itr, *lowest)) lowest = &(*itr);
                                itr++;
                            }
                            num = common::sector_numbers_fast[sector_fast(*lowest)];
                        }
                        if (num) dependency[num] = true;
                    }
                }
            }
        }
        dependencies.emplace_back(dependency);
        needed_sectors.push_back(false);
    }
}

/**
 * Advanced recursive mkdir.
 * @param s path to create.
 * @param mode file permission.
 * @return error code, 0 on success.
 */
int mkpath(std::string s, mode_t mode) {
    size_t pre = 0, pos;
    std::string dir;
    int mdret;
    if (s[s.size() - 1] != '/') {
        // force trailing / so we can handle everything in loop
        s += '/';
    }
    while ((pos = s.find_first_of('/', pre)) != std::string::npos) {
        dir = s.substr(0, pos++);
        pre = pos;
        if (dir == "") continue; // if leading / first time is 0 length
        if ((mdret = mkdir(dir.c_str(), mode)) && errno != EEXIST) {
            return mdret;
        }
    }
    return 0;
}

int parse_config(const string &filename, set<point, indirect_more> &points, string &output) {
    
    FILE *config_file = fopen(filename.c_str(), "r");
    if (config_file == nullptr) {
        cerr << "No " << filename << ", exiting!" << endl;
        return 1;
    }

    string pid_folder;
    char hostname[64];
    gethostname(hostname, 64);
    pid_folder = "/" + string(hostname) + "-" + to_string(getpid()); // master job in mpi_mode appends pid to database folder
    
    string cdatabase;
    bool started = false;
    char finput[LOAD_STR_SIZE] = "none";
    int bucket_value = 20;

    bool loaded_rules = false; // we should not load rules before preferred

    auto & fermat = common::fermat;
    auto & variables = common::variables;
    while (fgets(finput, sizeof(finput), config_file)) {
        string str = finput;
        if (str.substr(0, 2) == "##") {
            continue; // comment
        } else if (str.substr(0, 7) == "#fermat") {
            size_t pos = 7;
            while (str[pos] == ' ') pos++;
            fermat = str.substr(pos);
            fermat[fermat.find('\n')] = '\0';
        } else if (str.substr(0, 8) == "#threads") {
            continue; // skip
        } else if (str.substr(0, 9) == "#sthreads") {
            continue; // skip
        } else if (str.substr(0, 9) == "#lthreads") {
            continue; // skip
        } else if (str.substr(0, 9) == "#fthreads") {
            continue; // skip
        } else if (str.substr(0, 11) == "#compressor") {
            continue; // skip
        } else if (str.substr(0, 8) == "#storage") {
            continue; // skip
        } else if (str.substr(0, 7) == "#bucket") {
            continue; // skip
        } else if (str.substr(0, 5) == "#wrap") {
            continue; // skip
        } else if (str.substr(0, 7) == "#memory") {
            continue; // skip
        } else if (str.substr(0, 6) == "#clean") {
            continue; // skip
        } else if (str.substr(0, 9) == "#pos_pref") {
            size_t pos = 10;
            while (str[pos] == ' ') pos++;
            str = str.substr(pos);
            int pos_pref;
            s2i(str.c_str(), pos_pref);
            common::pos_pref = pos_pref;
        } else if (str.substr(0, 10) == "#variables") {
            size_t pos = 10;
            while (str[pos] == ' ') pos++;
            char variables_temp[COEFF_BUF_SIZE];
            strcpy(variables_temp, (str.substr(pos)).c_str());
            char *begin = variables_temp;
            char *now = variables_temp;
            bool mode_right = false;
            string left;
            unsigned short count_replaced_vars = 0;
            while (*now != '\0') {
                if (*now == '\n') { *now = ','; }
                if (*now == ',') {
                    *now = '\0';
                    if (mode_right) { // now making a variable replacement rule
                        string right = begin;
                        mode_right = false;
                        auto itr = common::var_values_from_arv.find(left);
                        if (itr != common::var_values_from_arv.end()) right = itr->second;
                        common::prt_replace[left] = right;
                        common::variable_replacements.emplace(left, right);
                        left = "";
                    } else { // just adding a variable
                        left = begin;
                        auto itr = common::var_values_from_arv.find(left);
                        if (itr == common::var_values_from_arv.end()) {
                            #if defined(FMPQ) || defined(FloatR) || defined(PRIME)
                            cout << "variable: " << left << " does not get a numeric value." << endl;
                            exit(0);
                            #endif
                            variables += begin;
                            variables += '\n';
                        } else {
                            string right = itr->second;
                            common::prt_replace[left] = right;
                            common::variable_replacements.emplace(left, right);
                        }
                        left = "";
                    }
                    now++;
                    begin = now;
                } else if (*now == ' ') {
                    now++;
                    begin++;
                } else if ((*now == '-') && (now[1] == '>')) { // left side of a rule
                    *now = '\0';
                    left = begin;
                    now++;
                    now++;
                    begin = now;
                    mode_right = true;
                } else now++;
            }
         } else if (str.substr(0, 9) == "#database") {
            if (common::path == "") {
                size_t pos = 9;
                while (str[pos] == ' ') pos++;
                common::path = str.substr(pos);
                common::path.erase(common::path.find('\n'));
                if (common::suffix != "") {
                    common::path = common::path + "-" + ReplaceAll(common::suffix, "/", "_");
                }
                if (common::run_mode) {
                    if(access(common::path.c_str(), 0)) mkpath(common::path, 0777);
                    // clean on database
                    if(!common::run_sector) {
                        for (const auto & file : std::filesystem::directory_iterator(common::path)) {
                            if(common::run_mode==1) file_remove(file.path());
                            else {
                                int status = 0;
                                fstream db(file.path(), fstream::in);
                                db >> status;
                                db.close();
                                if(status==0) file_remove(file.path());
                            }
                        }
                    }
                }
                
            } else {
                cout<<"#database setting ignored"<<endl;
            }
        } else if (str.substr(0, 6) == "#prime") {
            #if defined(PRIME) || defined(FloatR)
            if (!COEFF::prime) {
                int pos = 7;
                while (str[pos] == ' ') pos++;
                str = str.substr(pos);
                sscanf(str.c_str(), "%hu", &COEFF::prime_number);
                if (COEFF::prime_number > 127) {
                    cout << "Option #prime ignored: index for a prime should be in range from 0 to 127, refer to primes.cpp" << endl;
                } else {
                    COEFF::init_prime(COEFF::prime_number);
                }
            }
            #endif
        } else if (str.substr(0, 8) == "#keepall") {
            common::keep_all = true;
        } else if (str.substr(0, 6) == "#small") {
            continue;
        } else if (str.substr(0, 7) == "#allIBP") {
            common::all_ibps = true;
        } else if (str.substr(0, 12) == "#no_presolve") {
            common::disable_presolve = true;
        } else if (str.substr(0, 6) == "#start") {
            if (started) {
                cerr << "extra #start directive, exiting" << endl;
                fclose(config_file);
                return 1;
            }
            
            #if defined(FMPQ)
            COEFF::init_ctx_q();
            #elif defined(FloatR)
            COEFF::init_ctx_f();
            #elif defined(FlintM)
            auto vars = variables;
            int n = std::count(vars.begin(), vars.end(), '\n');
            if(n==0) { n = 1; vars = "_"; }
            COEFF::xs = (const char **)malloc(sizeof(const char *)*n);
            char cvs[vars.length()+1];
            strcpy(cvs,vars.c_str());
            char *pos = cvs;
            char *end = pos;
            for(int i=0; i<n; i++) {
                while(*end != '\n') ++end;
                *end = '\0';
                int len = std::distance(pos,end);
                char* nstr = (char*)malloc(sizeof(char)*(len+1));
                strcpy(nstr,pos);
                COEFF::xs[i] = nstr;
                COEFF::vs.push_back(nstr);
                pos = end;
                pos++;
            }
            COEFF::init_ctx_m(n);
            #endif
            
            if (fermat == "") {
                cerr << "No proper #fermat setting, exiting" << endl;
                fclose(config_file);
                return 1;
            }
            if (common::path == "") {
                common::path = "db"+common::config_file;
            }
            if (common::t1 <= 0) {
                cerr << "No proper #threads setting, exiting" << endl;
                fclose(config_file);
                return 1;
            }

            {
                if (iniCalc(fermat.c_str(), variables.c_str()) != 0) {
                    cerr << "Error starting fermat, exiting" << endl;
                    fclose(config_file);
                    return 1;
                }
            }

            started = true;
        } else {
            if (!started) {
                cerr << "#start directive missing, exiting" << endl;
                fclose(config_file);
                return 1;
            }

            if (str.substr(0, 7) == "#folder") {
                size_t pos = 7;
                while (str[pos] == ' ') pos++;
                str = str.substr(pos);
                str.erase(str.find('\n'));
                folder = str;
            } else if (str.substr(0, 8) == "#problem") {
                // problem should be parsed even for substitutions, however coefficients should not be loaded
                size_t pos = 8;
                unsigned int pn;
                while (str[pos] == ' ') pos++;
                str = str.substr(pos);

                str.erase(str.find('\n'));
                const char *r = str.c_str();

                int move = s2u(r, pn);
                while (*(r + move) == ' ') move++;

                if (!add_problem(pn, (string(r + move)).c_str())) {
                    fclose(config_file);
                    return 1;
                }

            } else if (str.substr(0, 5) == "#hint") {
                size_t pos = 6;
                while (str[pos] == ' ') pos++;
                common::hint = true;
                common::hint_path = str.substr(pos);
                common::hint_path.erase(common::hint_path.find('\n'));
                if (common::hint_path[0] != '/') {
                    common::hint_path = folder + common::hint_path;
                }
            } else if (str.substr(0, 6) == "#rules") {
                {
                    // we need to load preferred from rules anyway
                    size_t pos = 6;
                    while (str[pos] == ' ') pos++;
                    str = str.substr(pos);
                    str.erase(str.find('\n'));
                    const char *r;
                    string tmp;
                    if (str[0] == '/') {
                        r = str.c_str();
                    } else {
                        tmp = folder + str;
                        r = tmp.c_str();
                    }
                    if (!add_rules(r)) {
                        fclose(config_file);
                        return 1;
                    }
                }
                loaded_rules = true;
            } else if (str.substr(0, 7) == "#lbases") {
                {  // only for particular sector work, but still for the main thread to produce warnings
                    size_t pos = 7;
                    while (str[pos] == ' ') pos++;
                    str = str.substr(pos);
                    str.erase(str.find('\n'));
                    const char *r;
                    string tmp;
                    if (str[0] == '/') {
                        r = str.c_str();
                    } else {
                        tmp = folder + str;
                        r = tmp.c_str();
                    }
                    if (!add_lbases(r)) {
                        return 1;
                    }
                }
            } else if (str.substr(0, 7) == "#output") {
                {
                    size_t pos = 7;
                    while (str[pos] == ' ') pos++;
                    bool norepeat = false;
                    if (str[pos] == '!') {
                        ++pos;
                        norepeat = true;
                    }
                    str = str.substr(pos);
                    if (str.find('\n') != string::npos) { str.erase(str.find('\n')); }
                    if (str[0] == '/') {
                        output = str;
                    } else {
                        output = folder + str;
                    }
                    common::only_masters = false;
                    if (common::suffix != "") {
                        output = ReplaceAll(output, ".tables", "-" + ReplaceAll(common::suffix, "/", "_") + ".tables");
                    }
                    if (common::split_masters) {
                        char mnb[32];
                        sprintf(mnb,"(%u-%u)",common::master_number_min,common::master_number_max);
                        output = ReplaceAll(output, ".tables", "-" + string(mnb) + ".tables");
                    }
                    if (norepeat) {
                        FILE* file = fopen(output.c_str(), "r");
                        if (file != nullptr) {
                            fgetc(file);
                            if (feof(file)) {
                                printf("Tables are reserved! %s\n", output.c_str());
                            } else {
                                printf("Tables are already created! %s\n", output.c_str());
                            }
                            fclose(file);
                            return -1;
                        }
                    }
                }
//            } else if (str.substr(0, 9) == "#dep_file") {
//                if (sector == 0) {
//                    size_t pos = 9;
//                    while (str[pos] == ' ') pos++;
//                    common::dep_file = str.substr(pos);
//                    common::dep_file.erase(common::dep_file.find('\n'));
//                    if (common::dep_file[0] != '/') {
//                        common::dep_file = folder + common::dep_file;
//                    }
//                    if (!common::silent) { cout << "Dependency file will be saved to " << common::dep_file << endl; }
//                }
            } else if (str.substr(0, 8) == "#masters") {
                size_t pos = 8;
                while (str[pos] == ' ') pos++;
                str = str.substr(pos);
                if (str.find('\n') != string::npos) { str.erase(str.find('\n')); }
                if ((str[0] == '|') || (common::master_number_min)) {
                    // that's the split masters mode
                    if (lbases_loaded) {
                        cout<<"Lbases should be loaded after masters in split master mode"<<endl;
                        abort();
                    }

                    common::split_masters = true;
                    if (loaded_rules) {
                        cout << "Masters file for split master mode should be specified before rules!"<<endl;
                        abort();
                    }
                    char buf[256];
                    if (common::master_number_min) {
                        const char* poss = str.c_str();
                        if (*poss != '|') {
                            // no vertical line in config, just taking the path
                        } else {
                            ++poss;
                            while ((*poss) && (*poss!='|')) ++poss;
                            if (!(*poss)) {
                                cout << "Incorrect syntax for master file with master numbers passed (no second |)"<<endl;
                                abort();
                            }
                            str = str.substr(poss+1-str.c_str());
                        }
                        sscanf(str.c_str(),"%255s",buf);
                    } else {
                        const char* poss = str.c_str();
                        ++poss;
                        while ((*poss) && (*poss!='|') && (*poss!='-')) ++poss;
                        if (!(*poss)) {
                            cout << "Incorrect syntax for master file (no second |)"<<endl;
                            abort();
                        }
                        if (*poss == '-') {
                            sscanf(str.c_str(),"|%u-%u|%255s",&common::master_number_min,&common::master_number_max,buf);
                        } else {
                            sscanf(str.c_str(),"|%u|%255s",&common::master_number_min,buf);
                            common::master_number_max = common::master_number_min;
                        }
                        if (!common::master_number_min) {
                            cout << "Incorrect syntax for master file"<<endl;
                            abort();
                        }
                        if (common::master_number_max < common::master_number_min) {
                            cout << "Incorrect range of master integrals"<<endl;
                            abort();
                        }
                    }
                    string masters_list_file(buf);
                    if (masters_list_file[0] != '/') {
                        masters_list_file = folder + masters_list_file;
                    }
                    const char* r = masters_list_file.c_str();

                    initialize_dependencies();

                    set_masters_to_zero(r,true);

                } else {
                    {
                        if (str[0] == '/') {
                            output = str;
                        } else {
                            output = folder + str;
                        }
                    }
                    common::only_masters = true;
                }
            } else if ((str.substr(0, 9) == "#prefered") || (str.substr(0, 10) == "#preferred")) { // parse preferred right here
                { // on main for rules and in Laporta. it is needed for points creation
                    if (str.substr(0, 9) == "#prefered") {
                        cout << "WARNING: Please fix config file, this option should be named '#preferred'" << endl
                             << "For now '#prefered' is supported, but will be removed later" << endl;
                    }
                    if (loaded_rules) {
                        cout << "List of preferred master-integrals should be loaded before rules!";
                        return 1;
                    }
                    size_t pos = 10;
                    while (str[pos] == ' ') pos++;
                    str = str.substr(pos);
                    if (str.find('\n') != string::npos) str.erase(str.find('\n'));
                    if (str[0]!='/') str = folder + str;
                    const char* r = str.c_str();
                    int pref_count = load_preferred(r);
                    if (pref_count == -1) {
                        abort();
                    }
                }
            } else if (str.substr(0, 10) == "#integrals") { // parse integrals right here
                for (int i = 2; i < common::abs_max_sector; ++i) {
                    if (point::preferred[i].empty()) {
                        // no preferred, have to add them
                        vector<t_index> v = common::ssectors[i];
                        vector<t_index> c = corner(v);
                        point Corner = point_reference(c);
                        for (int j = 0; j != abs(common::pos_pref) + 1; ++j) {
                            set<point_fast> new_pref = level_points_fast(Corner, (common::pos_pref > 0) ? j : 0, (common::pos_pref > 0) ? 0 : j);
                            for (const auto &fast_pnt : new_pref) {
                                point p = point(fast_pnt);
                                vector<t_index> vec = p.get_vector();
                                bool to_add = true;
                                if (common::pos_pref > 0) {
                                    for (const t_index &index : vec) {
                                        if (index > static_cast<t_index>(common::pos_pref)) {
                                            to_add = false;
                                            break;
                                        }
                                    }
                                }
                                if (to_add) {
                                    add_single_preferred(vec);
                                }
                            }
                        }
                    }
                }

                {
                    common::abs_max_level = common::abs_min_level;

                    size_t pos = 10;
                    while (str[pos] == ' ') pos++;
                    str = str.substr(pos);
                    if (str.find('\n') != string::npos) { str.erase(str.find('\n')); }
                    FILE *integral_file;
                    if (str[0] == '/') {
                        integral_file = fopen((str).c_str(), "r");
                    } else {
                        integral_file = fopen((folder + str).c_str(), "r");
                    }
                    char load_string[LOAD_STR_SIZE] = "none";
                    string contents;
                    if (integral_file == nullptr) {
                        cerr << "File with integral list could not be opened, exiting" << endl;
                        return -1;
                    }
                    while (fgets(load_string, sizeof(load_string), integral_file)) {
                        contents += load_string;
                    }
                    contents = contents.substr(contents.find('{') + 1); //"}"
                    for (unsigned int i = 0; i != contents.size(); ++i) {
                        if (contents[i] == '\n') contents[i] = ' ';
                    }
                    for (unsigned int i = 0; i != contents.size(); ++i) {
                        if (contents[i] == '\r') contents[i] = ' ';
                    }
                    const char *all = contents.c_str();
                    int move = 0;
                    while (all[move] != '}') {
                        while (all[move] == ' ') move++;
                        if (all[move] != '{') {
                            cout << "error in integrals: " << all[move] << endl;
                            abort();
                        }
                        move++;
                        while (all[move] == ' ') move++;
                        unsigned int i;
                        move += s2u(all + move, i);
                        if (i != common::global_pn) {
                            cerr << "Unknown problem number " << i << endl;
                            fclose(integral_file);
                            return -1;
                        }

                        while (all[move] == ' ') move++;
                        if (all[move] != ',') {
                            cout << "error in integrals: ";
                            abort();
                        }
                        move++;
                        while (all[move] == ' ') move++;
                        vector<t_index> v;
                        move += s2v(all + move, v);

                        point p = point_reference(v);
                        if (!p.is_zero()) {
                            points.insert(p);
                            int level = positive_index(common::ssectors[p.s_number()]);
                            int s_number = p.s_number();
                            if (level > common::abs_max_level) { common::abs_max_level = level; }
                            if (s_number > common::abs_max_sector) {
                                cout << "Integral in a sector that does not exist!"<<endl;
                                abort();
                            }
                        }
                        equation::initial.emplace(v, p); //mapping from requested vectors to our points

                        while (all[move] == ' ') move++;
                        if (all[move] != '}') {
                            cout << "error in integrals: ";
                            abort();
                        }
                        move++;
                        while (all[move] == ' ') move++;
                        if (all[move] == '}') break;
                        if (all[move] != ',') {
                            cout << "error in integrals: ";
                            abort();
                        }
                        move++;
                    }
                    fclose(integral_file);
                }
            } else if (str.length() > 1 && str[0] != '\n' && str[1] != '\n'){
                cout << "Met strange option, ignoring it. Line:" << endl << str << endl;
            }
        }
    }
    if (!common::split_masters && common::master_number_min) {
        cout << "Master numbers set from command line, but #master option missing"<<endl;
        abort();
    }

    if (common::split_masters) {
        //recursvely building dependencies
        while (true) {
            bool changed = false;
            for (int i = 2; i!=common::abs_max_sector + 1; ++i) {
                for (int j = 2; j!=common::abs_max_sector + 1; ++j) {
                    if (dependencies[i][j]) {
                        // i depends on j, let's look further
                        for (int k = 2; k!=common::abs_max_sector + 1; ++k) {
                            if (dependencies[j][k] && !dependencies[i][k]) {
                                dependencies[i][k] = true;
                                changed = true;
                            }
                        }
                    }
                }
            }
            if (!changed) break;
        }

        if (common::dep_file != "") {
            std::ofstream out;
            out.open(common::dep_file);
            out<<"{{";
            for (int i = 2; i!=common::abs_max_sector + 1; ++i) {
                out<<"{"<<i<<",";
                vector<t_index>& v = common::ssectors[i];
                auto it = v.begin();
                out << "{" << int(*it);
                for (it++; it != v.end(); it++) {
                    out << "," << int(*it);
                }
                out << "}}";
                if (i!=common::abs_max_sector) out<<", "<<endl;
            }
            out<<"},"<<endl<<"{";
            for (int i = 2; i!=common::abs_max_sector + 1; ++i) {
                out<<"{";
                for (int j = 2; j!=common::abs_max_sector + 1; ++j) {
                    if (i==j) out<<1; else out<<dependencies[i][j];
                    if (j!=common::abs_max_sector) out<<",";
                }
                out<<"}";
                if (i!=common::abs_max_sector) out<<","<<endl;
            }
            out<<"}}"<<endl;
            out.close();
            cout<<"Dependency file saved!"<<endl;
        }

        // create the list of sectors where reduction is needed, set others to zero
        for (int i = 2; i!=common::abs_max_sector + 1; ++i) {
            bool needed = false;
            if (needed_sectors[i]) {
                needed = true;
            } else {
                for (int j = 2; j!=common::abs_max_sector + 1; ++j) {
                    if (dependencies[i][j] && needed_sectors[j]) {
                        needed = true;
                        break;
                    }
                }
            }
            if (needed) {
                //if (!sector) cout<<"Needed sector: "<<i<<endl;
            } else {
                if (!common::silent) cout<<"Setting sector to zero: "<<i<<endl;
                vector<vector<t_index> > orbit;
                vector<t_index> v = common::ssectors[i];
                Orbit(v, orbit, common::symmetries);
                for (auto v2: orbit) {
                    common::sector_numbers_fast[sector_fast(v2)] = 0;
                }
            }
        }

        // resetting some of requested points to zero
        for (auto &vp : equation::initial) {
            if (!common::sector_numbers_fast[sector_fast(common::ssectors[vp.second.s_number()])]) {
                points.erase(vp.second);
                vp.second = point();
            }
        }
    }

    fclose(config_file);
    if (common::t2 == 0) common::t2 = common::t1;
    #if defined(PRIME) || defined(FloatR)
    if (!COEFF::prime) {
        cerr << "no proper #prime or -prime found, exit!" << endl;
        fclose(config_file);
        return 1;
    }
    #endif
    return 0;
}

vector<COEFF> split_coeff(const string &s) {
    // and now here we get to splitting a string into coeffs at individual a[i]
    vector<COEFF> cc(MAX_IND + 1, CO_0);
    const char *pos = s.c_str();
    while ((*pos != '\0') && (*pos != '{')) ++pos; //}
    if (*pos != '{') {
        cout << "No opening bracket at coefficient start, perhaps old form of start file" << endl;
        abort();
    }
    ++pos;
    while (true) { // a cycle to find all pairs in coeff - coefficient and number of a
        while ((*pos != '\0') && (*pos != '{') && (*pos != '}')) ++pos; //{
        if (*pos == '}') break; // closing bracket of all the list
        if (*pos != '{') {
            cout << "No opening bracket inside coefficient";
            abort();
        } // for internal pair
        ++pos;
        string coeff;
        while (*pos != ',') {
            coeff += *pos;
            ++pos;
        }
        ++pos; //,
        int n;
        sscanf(pos, "%d", &n); //{
        while (*pos != '}') ++pos;
        if (coeff.empty()) {
            cc[n] = CO_0; // n is 0 means free coeff
        } else {
            coeff += "(1)";
            coeff = ReplaceAll(coeff, " ", "");
            if (*coeff.begin() == '+') {
                coeff = "0" + coeff;
            }
            coeff = calc(coeff);
            cc[n] = COEFF(coeff);
        }
        ++pos; // passing closing bracket; //{
        while ((*pos != ',') && (*pos != '}')) ++pos;
        if (*pos == ',') ++pos;
    }
    return cc;
}

void parseArgcArgv(int argc, char *argv[]) {
    bool t1 = false, t2 = false, lt1 = false, lt2 = false;
    for (int i = 1; i < argc; ++i) {
        if ((i + 1 != argc) && (!strcmp(argv[i],"-c"))) {
            common::config_file = string(argv[i + 1]);
            i++;
        } else if((i + 1 != argc) && (!strcmp(argv[i],"-sector"))) {
            string str(argv[i + 1]);
            s2i(str.c_str(), common::run_sector);
            if(common::run_sector) common::run_mode = 2;
            else { common::run_mode = 3; common::silent = true; }
            i++;
        } else if((i + 1 != argc) && (!strcmp(argv[i],"-t"))) {
            string str(argv[i + 1]);
            unsigned int t;
            s2u(str.c_str(), t);
            if(!t1) common::t1 = t;
            if(!t2) common::t2 = t;
            i++;
        } else if((i + 1 != argc) && (!strcmp(argv[i],"-t1"))) {
            string str(argv[i + 1]);
            s2u(str.c_str(), common::t1);
            t1 = true;
            i++;
        } else if((i + 1 != argc) && (!strcmp(argv[i],"-t2"))) {
            string str(argv[i + 1]);
            s2u(str.c_str(), common::t2);
            t2 = true;
            i++;
        } else if((i + 1 != argc) && (!strcmp(argv[i],"-lt"))) {
            string str(argv[i + 1]);
            unsigned int lt;
            s2u(str.c_str(), lt);
            if(!lt1) common::lt1 = lt;
            if(!lt2) common::lt2 = lt;
            i++;
        } else if((i + 1 != argc) && (!strcmp(argv[i],"-lt1"))) {
            string str(argv[i + 1]);
            s2u(str.c_str(), common::lt1);
            lt1 = true;
            i++;
        } else if((i + 1 != argc) && (!strcmp(argv[i],"-lt2"))) {
            string str(argv[i + 1]);
            s2u(str.c_str(), common::lt2);
            lt2 = true;
            i++;
        } else if((i + 1 != argc) && (!strcmp(argv[i],"-lmt1"))) {
            string str(argv[i + 1]);
            s2u(str.c_str(), common::lmt1);
            i++;
        } else if((i + 1 != argc) && (!strcmp(argv[i],"-lmt2"))) {
            string str(argv[i + 1]);
            s2u(str.c_str(), common::lmt2);
            i++;
        }  else if((i + 1 != argc) && (!strcmp(argv[i],"-ltm"))) {
            string str(argv[i + 1]);
            s2i(str.c_str(), common::ltm);
            i++;
        } else if((i + 1 != argc) && (!strcmp(argv[i],"-ft"))) {
            string str(argv[i + 1]);
            int t = 0;
            s2i(str.c_str(), t);
            if(t>1) flint_set_num_threads(t);
            i++;
        } else if (!strcmp(argv[i],"-oo")) {
            if(common::run_mode<1) common::o_output = 1;
        } else if (!strcmp(argv[i],"-silent")) {
            common::silent = true;
        } else if (!strcmp(argv[i],"-skip")) {
            common::skip_if_exist = true;
        } else if (!strcmp(argv[i],"-dbo")) {
            if(common::run_mode<1) common::run_mode = 1;
        } else if (!strcmp(argv[i],"-re")) {
            if(common::run_mode<2) common::run_mode = 2;
        } else if((i + 1 != argc) && (!strcmp(argv[i],"-suffix"))) {
            common::suffix = string(argv[i + 1]);
            i++;
        } else if ((i + 1 != argc) && (!strcmp(argv[i],"-var") || !strcmp(argv[i],"-variables"))) {
            char argv_tmp[strlen(argv[i+1])+1];
            strcpy(argv_tmp, argv[i+1]);
            char *curpos = argv_tmp;
            while (true) {
                string key, val;
                char *pos = curpos;
                while ((*pos != '\0') && (*pos != '=')) ++pos;
                if (*pos == '=') {
                    *pos = '\0';
                    key = curpos;
                    curpos = ++pos;
                }
                while ((*pos != '\0') && (*pos != '|') && (*pos != ',')) ++pos;
                bool done = (*pos == '\0');
                *pos = '\0';
                val = curpos;
                curpos = ++pos;
                common::var_values_from_arv.emplace(key, val);
                if(done) break;
            }
            i++;
        } else if((i + 1 != argc) && (!strcmp(argv[i],"-prime"))) {
            #if defined(PRIME) || defined(FloatR)
            sscanf(argv[i+1], "%hu", &COEFF::prime_number);
            COEFF::init_prime(COEFF::prime_number);
            #endif
            i++;
        } else if((i + 1 != argc) && (!strcmp(argv[i],"-fp"))) {
            string str(argv[i + 1]);
            s2u(str.c_str(), COEFF::fp);
            i++;
        }
    }
    if (common::config_file == "") {
        cout << "Config file not specified"<<endl;
        abort();
    }
}
