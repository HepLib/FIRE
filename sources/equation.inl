/** @file equation.inl
 *  @author Alexander Smirnov
 *
 *  This file is a part of the FIRE package.
 */

#ifdef MPQBin
size_t mpq_size(const mpq_class & q);
void mpq_write(const mpq_class & q, char* &pos);
#endif

#ifdef FlintC3
slong str_len3(const fmpz_poly_q_t p);
void write_str3(const fmpz_poly_q_struct* p, char* &pos);
#endif

template<class I>
void p_set(const point &p, unsigned int n, I termsB, I termsE, unsigned char level, unsigned short fixed_database_sector) {
    unsigned short dsector = (fixed_database_sector == 0) ? p.s_number() : fixed_database_sector;
    unsigned int string_size = 0;
    size_t nterms = 0;
    for (auto itr = termsB; itr != termsE; ++itr) ++nterms;
#ifdef PRIME
    string_size = n * sizeof(unsigned long long); // we will just put our number into the buffer
#elif defined(MPQStr)
    vector<string> svec(nterms);
    int svec_i = 0;
    for (auto itr = termsB; itr != termsE; ++itr) {
        svec[svec_i] = itr->second.s.get_str(62);
        string_size += (svec[svec_i].size() + 1); // note the 62 base
        svec_i++;
    }
    svec_i = 0;
#elif defined(MPQBin)
    for (auto itr = termsB; itr != termsE; ++itr) {
        string_size += mpq_size(itr->second.s);
    }
#elif defined(FlintX)
    vector<string> svec(nterms);
    int svec_i = 0;
    for (auto itr = termsB; itr != termsE; ++itr) {
        svec[svec_i] = itr->second.s.to_string();
        string_size += (svec[svec_i].size() + 1);
        svec_i++;
    }
    svec_i = 0;
#elif defined(FlintC3)
    for (auto itr = termsB; itr != termsE; ++itr) {
        string_size += str_len3(itr->second.s[0])+1;
    }
#elif defined(FMPQ)
    vector<string> svec(nterms);
    int svec_i = 0;
    for (auto itr = termsB; itr != termsE; ++itr) {
        char * cs = fmpq_get_str(NULL,62,itr->second.s[0]); // note the 62 base
        svec[svec_i] = cs;
        flint_free(cs); 
        string_size += (svec[svec_i].size() + 1);
        svec_i++;
    }
    svec_i = 0;
#elif defined(FlintC)
    vector<string> svec(nterms);
    int svec_i = 0;
    for (auto itr = termsB; itr != termsE; ++itr) {        
        svec[svec_i] = fmpz_poly_q_get_str2(itr->second.s[0]);
        string_size += (svec[svec_i].size() + 1);
        svec_i++;
    }
    svec_i = 0;
#elif defined(FlintM)
    vector<string> svec(nterms);
    int svec_i = 0;
    for (auto itr = termsB; itr != termsE; ++itr) {  
        svec[svec_i] = fmpz_mpoly_q_get_str(itr->second.s[0]);
        string_size += (svec[svec_i].size() + 1);
        svec_i++;
    }
    svec_i = 0;
#else
    for (auto itr = termsB; itr != termsE; ++itr) {
        string_size += (itr->second.s.size() + 1);
    }
#endif

    char *buf = static_cast<char*>(malloc(3 + n * sizeof(point) + string_size + 1));
    // the structure of the buffer
    // 1) n, 2 bytes
    // 2) level, 1 byte
    // 3) types of lower points, n/4 bytes in case of small, nothing otherwise
    // 4) points, differs in case of small, 24*n otherwise
    // 5) coeffs, 8*n in case of prime, string_size otherwise
    // 6) the 0 symbol, just in case

    if (!buf) {
        cout<<"Cannot malloc in p_set"<<endl;
        abort();
    }

    reinterpret_cast<unsigned short*>(buf)[0] = n; // bytes 1 and 2... can n be longer?
    buf[2] = level; // byte 3

    unsigned int points_size;
    unsigned short j = 0;
    for (auto itr = termsB; itr != termsE; ++itr) {
        reinterpret_cast<point *>(buf + 3)[j] = itr->first;
        ++j;
    }
    points_size = n * sizeof(point);

    uint32_t buf_size = 3 + points_size + string_size + 1;
    if (n) { // just some checks
        auto terms_itr = termsE;
        --terms_itr;
        if (p != terms_itr->first) {
            cout << "p_set error" << endl;
            cout << p << endl;
            cout << n << endl;
            for (auto itr = termsB; itr != termsE; ++itr) { cout << itr->first << endl; }
            abort();
        }
    }
    char *pos = buf + 3 + points_size; // preparing a string of coeffs
    for (auto itr = termsB; itr != termsE; ++itr) {
#ifdef PRIME
        *reinterpret_cast<unsigned long long *>(pos) = itr->second.n;
        pos += sizeof(unsigned long long);
#elif defined(MPQStr)
        auto const & s = svec[svec_i];
        strncpy(pos, s.c_str(), s.size());
        pos += s.size();
        svec[svec_i].clear();
        *pos = '|';
        ++pos;
        svec_i++;
#elif defined(FlintC3)
        write_str3(itr->second.s[0], pos);
        *pos = '|';
        ++pos;
#elif defined(FlintX) || defined(FMPQ) || defined(FlintC) || defined(FlintM)
        auto const & s = svec[svec_i];
        strncpy(pos, s.c_str(), s.size());
        pos += s.size();
        svec[svec_i].clear();
        *pos = '|';
        ++pos;
        svec_i++;
#elif defined(MPQBin)
        mpq_write(itr->second.s, pos);
#else
        strncpy(pos, itr->second.s.c_str(), itr->second.s.size());
        pos += itr->second.s.size();
        *pos = '|';
        ++pos;
#endif
    }
    *pos = '\0'; // do we really need it? we won't be interpreting the contents as a string any longer


    class VisitorImpl : public kyotocabinet::DB::Visitor {
        const char *visit_full(const char *kbuf, size_t ksiz, const char *vbuf, size_t vsiz, size_t *sp) {
            const point p = *reinterpret_cast<const point *>(kbuf);

            *sp = buf_size;

            unsigned short n = reinterpret_cast<unsigned short *>(buf)[0];
            char level = buf[2];

            unsigned short old_n = reinterpret_cast<const unsigned short *>(vbuf)[0];
            char old_level = vbuf[2];

            if ((old_n != 0) && (n == 0)) {
                //this means that we have a rule in a lower sector and now it is being marked as needed
                //we should not change the table in this case... but only check the level it is needed up to
                return kyotocabinet::DB::Visitor::NOP;
            }

            if ((old_level == 127) && (n) && (level != 127)) {
                // the point was marked as a master, but now we have a table
                cout << "SOMEHOW " << p << " IS NO LONGER A MASTER!!!" << endl;
                old_level = 126;
            }
            char new_level = old_level;
            if (level > new_level) new_level = level;
            buf[2] = new_level;
            return buf;
        }

        const char *visit_empty(const char *kbuf, size_t ksiz, size_t *sp) {
            *sp = buf_size;
            return buf;
        }
    public:
        char * buf;
        uint32_t buf_size;
    } visitor;

    visitor.buf = buf;
    visitor.buf_size = buf_size;
    common::points[dsector]->accept(p.ww, sizeof(point), &visitor, true);
    free(buf);
}
