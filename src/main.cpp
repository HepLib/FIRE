#include "common.h"
#include "parser.h"
#include "functions.h"
#include <cstdio>
#include <ifaddrs.h>
#include <cstring>
#include <arpa/inet.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <sys/socket.h>
#include <omp.h>
#include <filesystem>

static const pc_pair_ptr_vec empty_terms;
static map<point, pc_pair_ptr_vec> temp_db;
void p_set_temp(const point &p, const pc_pair_ptr_vec &terms) {
    temp_db.emplace(p, terms);
}
void p_get_temp(const point &p, pc_pair_ptr_vec &terms) {
    auto itr = temp_db.find(p);
    if (itr == temp_db.end()) {
        terms.clear();
    } else {
        terms = itr->second;
    }
}

int mkpath(std::string s, mode_t mode);

inline bool file_exists(string fn) {
    return (access(fn.c_str(),F_OK)!=-1);
}

void run_child(int sector);
int main(int argc, char *argv[]) {

    // parser silent first
    for(int i=1; i<argc; i++) {
        if(!strcmp(argv[i], "-silent")) {
            common::silent = true;
            break;
        }
    }
    
    // parser -h/--help then
    for(int i=1; i<argc; i++) {
        if(!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
            cout << "---------------------------------------------------------------" << endl;
            cout << "Version: " << common::version << endl;
            cout << "------------------" << endl;
            cout << "Supported options: " << endl;
            cout << "------------------" << endl;
            cout << "-c <fn>: to provide config file name." << endl;
            cout << "-prime <num>: set prime index, higher priority than config file" << endl;
            cout << "-t1 <num>: set threads for forward-stage." << endl;
            cout << "-t2 <num>: set threads for backward-stage." << endl;
            cout << "-t <num>: to set both -t1/-t2, lower priority than -t1/-t2." << endl;
            cout << "-lt1 <num>: set 2nd threads in each t1 run." << endl;
            cout << "-lt2 <num>: set 2nd threads in each t2 run." << endl;
            cout << "-lt <num>: to set both -lt1/-lt2, lower priority than -lt1/-lt2." << endl;
            cout << "-lmt1 <num>: set limit to activate -lt1." << endl;
            cout << "-lmt2 <num>: set limit to activate -lt2." << endl;
            cout << "-ltm <numb>: to set lt mode in forward-stage." << endl;
            cout << "-dbo: DB exported to and imported from file to save memory." << endl;
            cout << "-re: will continue to run from the saving point (imply -dbo)." << endl;
            cout << "-skip: just by pass if the .tables exists." << endl;
            cout << "-var <x1=v1,x22=v2>: to set x1->v1 and x2->v2, seperated by ',' or '|'," << endl;
            cout << "    e.g., -var \"d=11/5,m=1\", higher priority than config file." << endl;
            cout << "-variables: just the same as -var." << endl;
            cout << "-suffix <abc>: append '-abc' to filename of db/tables." << endl;
            cout << "-sector <num>: only run the sector where num can be 0, >0 or <0." << endl;
            cout << "    using num=0 to list the needed sector in next run." << endl;
            cout << "-ft <num>: call flint_set_num_threads(num) if num>1." << endl;
            cout << "-oo: the same output as official FIRE, i.e., 4/2 may appear." << endl;
            #if defined(FloatR)
            cout << "-fp <num>: set 2^fp as float precision." << endl;
            #endif
            cout << "---------------------------------------------------------------" << endl;
            return 0;
        }
    }
    
    char current[PATH_MAX];
    if (!getcwd(current, PATH_MAX)) {
        cout << "Can't get current dir name" << endl;
        return 1;
    }
    string scurrent = string(current);
    string srun = string(argv[0]);
    auto rfound = srun.rfind("/");

    if (rfound != std::string::npos) {  // running with full path
        srun = srun.substr(0, rfound+1);
    } else { //relative path, using current dir
        srun = scurrent + "/" + srun;
    }
    #if defined(__APPLE__) && defined(__aarch64__)
    common::fermat = srun + "/../usr/Ferm7a/fer64";
    #elif defined(__APPLE__)
     common::fermat = srun + "/../usr/Ferm7i/fer64";
    #else
    common::fermat = srun + "/../usr/Ferl7/fer64";
    #endif

    static_assert(sizeof(point) == POINT_SIZE, "Strange size of point class");

    parseArgcArgv(argc, argv);

    set<point, indirect_more> points;
    string output;
    int res = parse_config(common::config_file + ".config", points, output);
    if (res > 0) return -1;
    else if (res == -1) return 0;
    
    if(common::skip_if_exist && file_exists(output)) {
        if(!common::silent) cout << "Skip: " << output << " found and exit." << endl;
        if(common::run_mode) system(("rm -rf "+common::path).c_str());
        return 0;
    }
    
    if (!common::silent) {
        cout << "----------------------------------------" << endl;
        cout << "Version: " << common::version << endl;
        if(COEFF::prime) {
            cout << "Prime Index: " << COEFF::prime_number << endl;
            cout << "Prime Number: " << COEFF::prime << endl;
        }
        #if defined(FloatR)
        cout << "Float Precision: 2^" << COEFF::fp << endl;
        #endif
        cout << "Forward : T1/LT1/LMT1 = " << common::t1 << "/" << common::lt1 << "/" << common::lmt1 << endl;
        cout << "Backward: T2/LT2/LMT2 = " << common::t2 << "/" << common::lt2 << "/" << common::lmt2 << endl;
        if(common::prt_rule_counter) cout << "Rules: " << common::prt_rule_counter << endl;
        cout << "Mode: DB-" << common::run_mode << ", LT1-" << common::ltm << endl;
        if(!common::prt_replace.empty()) {
            cout << "Parameters: " << endl;
            for(auto kv : common::prt_replace) cout << "  " << kv.first << " -> " << kv.second << endl;
        }
        if(COEFF::vs.size()>0 && COEFF::vs[0]!="_") {
            cout << "Variables: " << COEFF::vs[0];
            for(int i=1; i<COEFF::vs.size(); i++) cout << ", " << COEFF::vs[i];
            cout << endl;
        }
        // TODO: to print variables HERE
        if(common::run_mode) cout << "Database: " << common::path << endl;
        if(common::only_masters) cout << "Masters: " << output << endl;
        else cout << "Output: " << output << endl;
    }
    
    if(common::run_sector) {
        run_child(common::run_sector);
        return 0;
    }
    
    if(!common::silent) cout << "----------------------------------------" <<  endl;

    map<sector_count_t, set<point> > needed;
    for (const auto &pnt : points) add_needed(needed, pnt);

    for (const auto &item : needed) {
        if(common::run_mode>1 && read_status(item.first)!=0) continue;
        open_database(item.first);
        for (const auto &pnt : item.second) {
            pc_pair_ptr_vec t;
            p_set(pnt, std::move(t), 126);
        }
        close_database(item.first);
    }

    if (!points.empty()) Evaluate();
    
    if (!common::silent) cout << "----------------------------------------" << endl;
    
    set<point, indirect_more> masters;

    if (!common::only_masters) {
        // clean on databases
        if (!common::silent) cout << "Clean DB ..." << endl;
        for(auto & kv : DBM) {
            auto itr = needed.find(kv.first);
            if(itr == needed.end()) kv.second.clear();
            else {
                open_database(itr->first);
                DBM[itr->first].keep = true;
                vector<point> vp;
                kv.second.lower.swap(vp);
                kv.second.umap.clear();
                set<point> pts;
                auto & pmap = kv.second.pmap;
                for (auto & kv2 : pmap) {
                    if(itr->second.find(kv2.first)==itr->second.end()) pts.insert(kv2.first);
                }
                for(auto & pnt : pts) pmap.erase(pnt);
                close_database(itr->first);
            }
        }
    
        //reading expressions from all databases, filling temporary tables, creating list of masters;
        for (const auto &item : needed) {
            for (const auto &pnt : item.second) {
                pc_pair_ptr_vec terms;
                p_get(pnt, terms);
                if (terms.empty()) {
                    masters.insert(pnt);
                } else {
                    p_set_temp(pnt, terms);
                    for (const auto &term : terms) {
                        if (!((*term).first == pnt)) {
                            masters.insert((*term).first);
                            p_set_temp((*term).first, empty_terms);
                        }
                    }
                }
            }
        }
    } else { // now the only masters way
        // we need to locate masters and add them to the final list
        if(!common::silent) cout << "Identifying Master Integrals ..." << endl;
        for (sector_count_t test_sector = 2; test_sector <= common::abs_max_sector; ++test_sector) {
            if (!database_exists(test_sector)) continue;
            auto const & pm = DBM[test_sector].pmap;
            for(auto const & kv : pm) {
                if(kv.second.first == 127) masters.insert(kv.first);
            }
        }
    }
    
    if(!common::silent) cout << "Master Integrals: " << masters.size() << endl;
    points.clear();
    points = masters;

    //replacing vectors in points with their original vectors in orbits
    if (!common::only_masters) {
        for (auto &integral : equation::initial) {
            point p = integral.second;
            if (!p.is_zero()) {  // the original vector was not mapped into zero
                if (integral.first != p.get_vector()) { // the vectors differ. we do not have such a point
                    point new_p(integral.first);
                    pc_pair_ptr_vec terms;
                    p_get_temp(p, terms);
                    if (!terms.empty()) {   // there is some table entry for our point, we are making a copy
                        const COEFF & c = terms.back()->second;
                        terms.pop_back();
                        terms.push_back(make_pc_ptr(new_p, c));
                    } else {    // there is no entry for our point, it is a master, sending the new one to the old one
                        terms.push_back(make_pc_ptr(p, CO_1m));
                        terms.push_back(make_pc_ptr(new_p, CO_1));
                    }
                    p_set_temp(new_p, terms);
                    points.insert(new_p);
                } else { //such a point already exists
                    points.insert(p);
                }
            } else {  // it is a zero point, mapping it to zero
                vector<t_index> vv = integral.first;

                point new_p(vv,0,-2); // get it to sector 1 without changes
                pc_pair_ptr_vec terms;
                terms.push_back(make_pc_ptr(new_p, CO_1));
                
                p_set_temp(new_p, terms);
                points.insert(new_p);
            }
        }
    }

    // everything done, saving tables
    if (!common::silent && !common::only_masters) cout << "Saving Output: " << output << endl;
    
    ostringstream oss;

    if(common::o_output) { // original FIRE output
        oss << "{" << endl << "    {" << endl;
        for (auto itr = points.begin(); itr != points.end(); ++itr) {
            oss << "        {" << itr->number() << "," << endl;
            oss << "            {" << endl;
            pc_pair_ptr_vec terms;
            p_get_temp(*itr, terms);
            if (terms.empty()) {
                oss << "{" << itr->number() << ",\"1\"}}}";
            } else {
                for (unsigned int i = 0; i != terms.size() - 1; ++i) {
                    oss << "                {" << terms[i]->first.number() << "," << "\"";
                    #ifdef PRIME
                    COEFF res;
                    _div_neg_(res, terms[i]->second, terms.back()->second);
                    oss << res;
                    #else
                    oss << "-(" << terms[i]->second << ")/(" << terms.back()->second << ")";
                    #endif
                    oss << "\"" << "}";
                    if (i + 2 != terms.size()) {
                        oss << ",";
                    }
                    oss << endl;
                }
                oss << "            }" << endl;
                oss << "        }";
            }
            itr++;
            if (itr != points.end()) {
                oss << "," << endl;
            }
            itr--;
        }
        oss << endl << "    }," << endl;

        oss << "    {" << endl;
        for (auto itr = points.begin(); itr != points.end(); ++itr) {
            oss << "        {" << itr->number() << "," << *itr << "}";
            itr++;
            if (itr != points.end()) {
                oss << ",";
            }
            itr--;
            oss << endl;
        }
        oss << "    }" << endl << "}" << endl;
    } else {
        oss << "{" << endl << "    {" << endl;
        for (auto itr = points.begin(); itr != points.end(); ++itr) {
            oss << "        {" << itr->number() << "," << endl;
            oss << "            {" << endl;
            pc_pair_ptr_vec terms;
            p_get_temp(*itr, terms);
            if (terms.empty()) {
                oss << "{" << itr->number() << ",\"1\"}}}";
            } else {
                for (unsigned int i = 0; i != terms.size() - 1; ++i) {
                    oss << "                {" << terms[i]->first.number() << "," << "\"";
                    COEFF res;
                    _div_neg_(res, terms[i]->second, terms.back()->second);
                    oss << res;
                    oss << "\"" << "}";
                    if (i + 2 != terms.size()) {
                        oss << ",";
                    }
                    oss << endl;
                }
                oss << "            }" << endl;
                oss << "        }";
            }
            itr++;
            if (itr != points.end()) {
                oss << "," << endl;
            }
            itr--;
        }
        oss << endl << "    }," << endl;

        oss << "    {" << endl;
        for (auto itr = points.begin(); itr != points.end(); ++itr) {
            oss << "        {" << itr->number() << "," << *itr << "}";
            itr++;
            if (itr != points.end()) {
                oss << ",";
            }
            itr--;
            oss << endl;
        }
        oss << "    }" << endl << "}" << endl;
    }
    
    {
        fstream out;
        out.open(output, fstream::out);
        out << oss.str() << flush;
        out.flush();
        out.close();
        
        bool ok = true;
        if(file_exists(output)) {
            string str_1st;
            fstream fs_1st(output, fstream::in);
            fs_1st >> str_1st;
            fs_1st.close();
            if(str_1st != "{") ok = false;
        }
        
        if((out.rdstate() & std::ofstream::failbit) || !file_exists(output) || !ok) {
            string tdir = std::filesystem::temp_directory_path().string();
            string tfn = tdir + "/" + output;
            fstream out;
            out.open(tfn, fstream::out);
            out << oss.str() << flush;
            out.flush();
            out.close();
            std::filesystem::rename(tfn, output);
        }
    }
    
    flint_cleanup();
    
    closeCalc();
    if(common::run_mode) system(("rm -rf "+common::path).c_str());
    if(!common::silent) cout << "----------------------------------------" <<  endl;
    return 0;
}
