/** @file common.h
 *  @author Alexander Smirnov
 *
 *  This file is a part of the FIRE package.
 *
 *  It contains multiple definitions of static variables, gathered in class common.
 */


#ifndef COMMON_H_INCLUDED
#define COMMON_H_INCLUDED

#if !defined(PRIME) && !defined(MPQ) && !defined(FlintX) && !defined(FMPQ) && !defined(FlintC) && !defined(FlintM)
#define PolyMode
#endif

#ifdef MPQ
#include <gmpxx.h>
#elif defined(FlintX)
#include "flint/fmpz_poly_qxx.h"
#elif defined(FMPQ)
#include <flint/fmpq.h>
inline void fmpq_set_si(fmpq_t res, slong p) { fmpq_set_si(res,p,1); }
#elif defined(FlintC)
#include "flint/fmpz_poly_q.h"
#elif defined(FlintM)
#include "fmpz_mpoly_q.h"
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <list>
#include <set>
#include <vector>
#include <unistd.h>
#include <algorithm>
#include <kcdb.h>
#include <kchashdb.h>
#include <kccachedb.h>
#include <kcstashdb.h>
#include <cstdint>
#include <lz4.h>
#include <lz4hc.h>
#include <fcntl.h>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>

/** Used for storing indices of integrals */
typedef char t_index;

#ifndef SMALL_POINT
constexpr size_t MAX_IND = {22};
#else
constexpr size_t MAX_IND = {18};
constexpr size_t BITS_PER_INDEX = {6};
#endif
/**< @brief
 * 22 indices should be enough for most problems;
 * 5-loop propagator = 20;
 * 6-loop bubble = 21;
 */

constexpr size_t MAX_THREADS = {64};
///<Maximal number of threads, 64 is just a number
constexpr size_t MAX_SOCKET_THREADS = {64};
///<Maximal number of socket threads (child communication), 64 is just a number
constexpr unsigned short MAX_SECTORS  = {128*256};
///<Number of sectors should fit into 2^15

/**
 * 128-but signed int with a standard name
 */
typedef __int128 int128_t;

/**
 * 128-but unsigned int with a standard name
 */
typedef unsigned __int128 uint128_t;


/** SECTOR type uses a bit for each coordinate, 1 being positing, 0 - negative.
 *   Virtual sectors also have a preceding 1 bit, corresponding to sector 1 and used for right-hand sides of rules
 */
typedef uint32_t SECTOR;


using namespace std;

/**
 * Type of compressor used in database.
 */
enum class t_compressor {
    C_NONE = 0,
    C_SNAPPY = 1,
    C_ZLIB = 2,
    C_LZ4FAST = 3,
    C_LZ4 = 4,
    C_LZ4HC = 5,
    C_ZSTD = 6
};

/**
 * Sectors that use the given sector in tables.
 */
static set<unsigned short> needed_for[MAX_SECTORS + 1];


#ifdef WITH_SNAPPY

#include <snappy.h>

/**
 * @brief Snappy Compressor extension for the kyotocabinet.
 */
class SnappyCompressor : public kyotocabinet::Compressor {
private:
    char *compress(const void *buf, size_t size, size_t *sp) {
        _assert_(buf && size <= MEMMAXSIZ && sp);
        *sp = snappy::MaxCompressedLength(size);
        char *out = new char[*sp + 1];
        snappy::RawCompress(static_cast<const char *>(buf), size, out, sp);
        out[*sp] = '\0';
        return out;
    }

    char *decompress(const void *buf, size_t size, size_t *sp) {
        _assert_(buf && size <= MEMMAXSIZ && sp);
        if (!snappy::GetUncompressedLength(static_cast<const char *>(buf), size, sp)) {
            return nullptr;
        }
        char *out = new char[*sp + 1];
        if (!snappy::RawUncompress(static_cast<const char *>(buf), size, out)) {
            delete[] out;
            return nullptr;
        }
        out[*sp] = '\0';
        return out;
    }
};
#endif

/**
* @brief LZ4 Fast Compressor extension for the kyotocabinet.
*/
class LZ4FastCompressor : public kyotocabinet::Compressor {
private:
    char *compress(const void *buf, size_t size, size_t *sp) override {
        _assert_(buf && size <= MEMMAXSIZ && sp);
        *sp = LZ4_COMPRESSBOUND(size);
        char *out = new char[*sp + 2];
        int written_bytes = LZ4_compress_fast(static_cast<const char *>(buf), out + 1, size, *sp, 5);
        if (written_bytes == 0) {
            delete[] out;
            fprintf(stderr, "Can't compress\n");
            return nullptr;
        }
        *sp = static_cast<size_t>(written_bytes);
        auto ratio = static_cast<unsigned char>(static_cast<unsigned int>(size) / static_cast<unsigned int>(*sp));
        ++ratio;
        ++(*sp);
        out[0] = ratio;
        out[*sp + 1] = '\0';
        return out;
    }

    char *decompress(const void *buf, size_t size, size_t *sp) override {
        _assert_(buf && size <= MEMMAXSIZ && sp);
        auto lim = static_cast<unsigned int>(size * static_cast<const unsigned char *>(buf)[0]);
        char *out = new char[lim];
        int decompressed_bytes = LZ4_decompress_safe(static_cast<const char *>(buf) + 1, out, size - 1, lim);
        if (decompressed_bytes < 0) {
            delete[] out;
            fprintf(stderr, "Can't decompress\n");
            return nullptr;
        }
        *sp = static_cast<size_t>(decompressed_bytes);
        out[*sp] = '\0';
        return out;
    }
};

/**
 * @brief LZ4 Compressor extension for the kyotocabinet.
 */
class LZ4Compressor : public kyotocabinet::Compressor {
private:
    char *compress(const void *buf, size_t size, size_t *sp) override {
        _assert_(buf && size <= MEMMAXSIZ && sp);
        *sp = LZ4_COMPRESSBOUND(size);
        char *out = new char[*sp + 2];
        int written_bytes = LZ4_compress_default(static_cast<const char *>(buf), out + 1, size, *sp);
        if (written_bytes == 0) {
            delete[] out;
            fprintf(stderr, "Can't compress\n");
            return nullptr;
        }
        *sp = static_cast<size_t>(written_bytes);
        auto ratio = static_cast<unsigned char>(static_cast<unsigned int>(size)/static_cast<unsigned int>(*sp));
        ++ratio;
        ++(*sp);
        out[0] = ratio;
        out[*sp + 1] = '\0';
        return out;
    }

    char *decompress(const void *buf, size_t size, size_t *sp) override {
        _assert_(buf && size <= MEMMAXSIZ && sp);
        auto lim = static_cast<unsigned int>(size * static_cast<const unsigned char *>(buf)[0]);
        char *out = new char[lim];
        int decompressed_bytes = LZ4_decompress_safe(static_cast<const char *>(buf) + 1, out, size - 1, lim);
        if (decompressed_bytes < 0) {
            delete[] out;
            fprintf(stderr, "Can't decompress\n");
            return nullptr;
        }
        *sp = static_cast<size_t>(decompressed_bytes);
        out[*sp] = '\0';
        return out;
    }
};

/**
 * @brief LZ4HC Compressor extension for the kyotocabinet.
 */
class LZ4HCCompressor : public kyotocabinet::Compressor {
private:
    char *compress(const void *buf, size_t size, size_t *sp) override {
        _assert_(buf && size <= MEMMAXSIZ && sp);
        *sp = LZ4_COMPRESSBOUND(size);
        char *out = new char[*sp + 2];
        int written_bytes = LZ4_compress_HC(static_cast<const char *>(buf), out + 1, size, *sp, 9);
        if (written_bytes == 0) {
            delete[] out;
            fprintf(stderr, "Can't compress\n");
            return nullptr;
        }
        *sp = static_cast<size_t>(written_bytes);
        auto ratio = static_cast<unsigned char>(static_cast<unsigned int>(size) / static_cast<unsigned int>(*sp));
        ++ratio;
        ++(*sp);
        out[0] = ratio;
        out[*sp + 1] = '\0';
        return out;
    }

    char *decompress(const void *buf, size_t size, size_t *sp) override {
        _assert_(buf && size <= MEMMAXSIZ && sp);
        auto lim = static_cast<unsigned int>(size * static_cast<const unsigned char *>(buf)[0]);
        char *out = new char[lim];
        int decompressed_bytes = LZ4_decompress_safe(static_cast<const char *>(buf) + 1, out, size - 1, lim);
        if (decompressed_bytes < 0) {
            delete[] out;
            fprintf(stderr, "Can't decompress\n");
            return nullptr;
        }
        *sp = static_cast<size_t>(decompressed_bytes);
        out[*sp] = '\0';
        return out;
    }
};

// declarations of functions
/**
 * Calculate number of positive indices in a sector.
 * @param v vector of indices
 * @return the number of positive indices
 */
int positive_index(vector<t_index> v);

/**
 * Calculate the sector corresponding the current integral.
 * @param v vector of indices
 * @return corresponding vector of 1 and -1
 */
vector<t_index> sector(const vector<t_index> &v);

/**
 * Calculate corner - the corner integral corresponding the current integral.
 * @param v vector of indices
 * @return corresponding vector of 1 and 0
 */
vector<t_index> corner(const vector<t_index> &v);

/**
 * Calculate degree. Degree is the shift from the corner integral.
 * @param v vector of indices
 * @return corresponding vector of degrees
 */
vector<t_index> degree(const vector<t_index> &v);

/**
 * Calculate level.
 * @param v vector of indices
 * @return complexity of a point - the number of dots and the number of irreducible denominators.
 */
pair<unsigned int, unsigned int> level(const vector<t_index> &v);

/**
 * Print vector of indices.
 * @param v vector of indices
 */
void print_vector(const vector<t_index> &v);

/**
 * Print indices stored in SECTOR variable.
 * @param sf compressed vector of indices.
 */
void print_sector_fast(const SECTOR &sf);

/**
 * Generate symmetry orbit of a point.
 * @param v point for which we generate orbit
 * @param orbit orbit to be updated
 * @param sym symmetries which we use for generation
 */
void Orbit(const vector<t_index> &v, vector<vector<t_index> > &orbit, vector<vector<vector<t_index> > > &sym);

/**
 * Create an ordering for the sector. Called in the initialization stage.
 */
void make_ordering(t_index *, const vector<t_index> &);

/**
 * List all sectors to be considered during solution.
 * @param d dimension
 * @param positive maximal index that can be positive
 * @param positive_start minimal index that can be positive
 * @return vector of all sectors
 */
vector<vector<t_index> > all_sectors(unsigned int d, int positive, int positive_start);

/**
 * Self-made conversion from int to string with leading zeroes.
 * @param i number
 * @return resulting string
 */
string int2string(int i);

// fast versions
/**
 * Calculate the sector corresponding the current integral and store it in compressed manner.
 * @param v vector of indices
 * @return compressed sector
 */
SECTOR sector_fast(const vector<t_index> &v);

/**
 * @brief Contains only static members used globally in FIRE.
 */
class common {
public:
    /**
     * Path to the file to save sector depenencies
     */
    static string dep_file;
    /**
     * Setting that allows to tweak the choice of master integrals. See #pos_pref in configuration file.
     */
    static int pos_pref;
    /**
     * Array of sector numbers, takes SECTOR as index.
     */
    static unique_ptr<unsigned short[]> sector_numbers_fast;
    /**
     * Array of pointers to ordering matrices, NULL initially.
     */
    static unique_ptr<unique_ptr<t_index[]>[]> orderings_fast;
    /**
     * Number of indices.
     */
    static unsigned short dimension;

    /**
     * True if this is Ftool executable.
     */
    static bool ftool;

    /**
     * True to disable IBP presolving without index substitution
     */
    static bool disable_presolve;

    /**
     * Self-describing. To pass them to FLAME.
     */
    static bool variables_set_from_command_line;
    /**
     * Prefix added to tables in PRIME mode.
     */
    static string tables_prefix;
    /**
     * Values of variables from config or command line.
     */
    static vector<unsigned short> var_values_from_arv;

    /**
     * Compressor choice for database entries.
     */
    static t_compressor compressor;
    /**
     * Compressor level for database entries.
     */
    static int compressor_level;
    /**
     * Compressor for database entries.
     */
    static unique_ptr<kyotocabinet::Compressor> compressor_class;
    /**
     * The number of a non-sector - the highest one that is used for all global symmetry mappings.
     */
    static int virtual_sector;

    /** @name Limits for levels and sectors.
    *  Maximal and minimal levels and sectors encountered.
    */
    /**@{*/
    /**
     * Positive index of lowest level.
     */
    static int abs_min_level;
    /**
     * Highest level.
     * abs_max_level cannot be higher than 15 - this is the maximum for the 6-loop bubble
     */
    static int abs_max_level;
    /**
     * Maximal sector number. Minimal sector number is always 2 in our enumeration
     */
    static int abs_max_sector;
    /**@}*/
    /**
     * Indication that only a part of masters will be used during reduction and substitutions.
     */
    static bool split_masters;
    /**
     * The minimal number of the master-integral that is not set to zero in split_masters mode.
     */
    static unsigned int master_number_min;
    /**
     * The maximal number of the master-integral that is not set to zero in split_masters mode.
     */
    static unsigned int master_number_max;
    /** @name Database wrapper members.*/
    /**@{*/
    /**
     * Mutex that controls access to wrapper tar file
     */
    static mutex wrapper_mutex;
    /**
     * Database for storing other databases
     */
    static kyotocabinet::HashDB wrapper_database;
    /**
     * Flag that corresponds to \#wrap in config.
     * True if \#wrap is used, false otherwise.
     */
    static bool wrap_databases;
    /**@}*/

    /**
     * Flag that corresponds to selection of \#masters option in config.
     * True if we used \#masters, False if \#output.
     */
    static bool only_masters;

    /**
     * Flag controlling if we lock the database.
     */
    static bool nolock;

    /**
     * If set to false, FIRE will print much more verbose information about work being done.
     * False by default.
     */
    static bool silent;

    /**
     * Database handlers.
     */
#ifdef DISK_DB
    static kyotocabinet::HashDB *points[MAX_SECTORS + 1];
#elif defined(Stash_DB)
    static kyotocabinet::StashDB *points[MAX_SECTORS + 1];
#else
    static kyotocabinet::CacheDB *points[MAX_SECTORS + 1];
#endif

    /** @name Database bucket settings and sizes.*/
    /**@{*/
    /**
     * Property of database, see kyotocabinet documentation for details.
     */
    static int buckets[MAX_SECTORS + 1];
    /**
     * Value that stores extra information about database. It's related to
     * buckets in the following way: buckets_full[i] == 2^buckets[i].
     */
    static int64_t buckets_full[MAX_SECTORS + 1];
    /**@}*/

    /** True if we use RAM database. */
    static bool memory_db;

    /** @name Pipes for expression communication.*/
    /**@{*/
    /**
     * File stream that child is writing to.
     */
    static FILE *child_stream_from_child;
    /**
     * File stream that child is reading from.
     */
    static FILE *child_stream_to_child;
    /**
     * Flag that tells binary that it should receive answers from child. In other words, that it's a parent.
    */
    static bool receive_from_child;
    /**
     * Flag that tells binary that it should send answers to parent. In other words, that it's a child working in no-separate fermat mode.
     */
    static bool send_to_parent;
    /**@}*/

    /**
     * True if we keep all entries, false otherwise.
     */
    static bool keep_all;

    /** True if we are remote worker, that is FLAME binary called with sector 0. */
    static bool remote_worker;

    /**
     * True if FIRE was run in parallel mode.
     * That is a call from the MPI binary or simply by providing the -parallel option.
     * The result is separation of database paths and semapthore names.
     */
    static bool parallel_mode;

    /** @name Stored paths to folders and files.*/
    /**@{*/
    /**
     * Path FIRE folder with input and output.
     */
    static string FIRE_folder;

    /**
     * Path to configuration file.
     */
    static string config_file;

    /**
     * Path to databases.
     */
    static string path;

    /**
     * Path to the so-called storage, copies of databases (if we use them).
     */
    static string cpath;

    /**
     * Path to folder with hints.
     */
    static string hint_path;
    /**@}*/


    /**
     * Database tuning, see kyotocabinet documentation for details.
     */
    static int msiz;

    /** @name Variables for statistics.*/
    /**@{*/
    /**
     * Time spend for expression simplification
     */
    static atomic<long long> fermat_time;
    /**
     * Total time spent in thread.
     */
    static atomic<long long> thread_time;
    /**@}*/

    /**
     * This flag corresponds to usage of \#storage option in config.
     */
    static bool cpath_on_substitutions;

    /**
     * Whether to clean temporary databases after work
     */
    static bool clean_databases;

    /**
     * True if we use all IBPs. See \#allIBP in configuration file.
     */
    static bool all_ibps;

    /**
     * Inverse orderings in sectors (matrices)
     */
    static vector<vector<vector<t_index> > > iorderings;

    /**
     * Maps numbers to sectors (as vectors).
     */
    static vector<vector<t_index> > ssectors;

    /**
     * Global symmetries.
     */
    static vector<vector<vector<t_index> > >  symmetries;

    /**
     * Set of sectors lower than others in their level.
     */
    static set<SECTOR> lsectors;

    /**
     * Lee bases.
     */
    static map<unsigned short,
            vector<
                    pair<
                            vector<
                                    pair<
                                            vector<t_index>,
                                            pair<short, bool>
                                    >
                            >,
                            vector<
                                    pair<
                                            string,
                                            vector<
                                                    pair< vector<t_index>, short>
                                            >
                                    >
                            >
                    >
            >
    > lbases;

    /**
     * The diagram number in use.
     */
    static unsigned int global_pn;

    /** @name Various thread counts.*/
    /** By default threads_number == fthreads_number == sthreads_number*/
    /**@{*/
    /**
     * Number of threads.
     */
    static unsigned int threads_number;

    /**
     * Number of level workers inside a sector. Equals to 1 by default.
     */
    static unsigned int lthreads_number;

    /**
     * Number of threads working during substitution stage.
     * should be decreased in case of
     * terminate called after throwing an instance of 'std::runtime_error'
     * what():  pthread_key_create
     */
    static unsigned int sthreads_number;

    /**
     * Number of fermat processes.
     */
    static unsigned int fthreads_number;

    /**
     * Number of fermat separate queues. Equals to 1 by default, which means all sectors use same fermat queue.
     */
    static unsigned int f_queues;
    /**@}*/

    /**
     * Number of iterations between printing information during reduce and substitution.
     * If equal to 0 we don't print anything.
     */
    static int print_step;
    /**
     * Prime number we use for modular arithmetic.
     */
    static uint64_t prime;
    /**
     * Index of prime number in primes array in primes.cpp.
     */
    static unsigned short prime_number;

    /**
     * Map of variable substitutions.
     */
    static map<string, string> variable_replacements;
    /**
     * True if \#small is in configuration file.
     */
    static bool small;
    /**
     * True if we use hints. See \#hint in configuration file.
     */
    static bool hint;

    /**
     * Thread which we use to listen for incoming connections.
     */
    static std::thread socket_listen;
};


// some more function declarations
/**
 * Check existence of database by its number.
 * @param number sector number
 * @return true if exists on disk
 */
bool database_exists(int number);

/**
 * Remove the database by its number.
 * @param number sector number
 */
void clear_database(int number);

/**
 * Open a database (either on disk, or in RAM).
  * @param number sector number
 */
void open_database(int number);

/**
 * Reopen database by number.
 * @param number sector number
 */
void reopen_database(int number);

/**
 * Close database by number.
 * @param number sector number
 */
void close_database(int number);

/**
 * Copy database by number to or from storage.
 * @param number
 * sector number
 * @param from_storage copy direction
 */
void copy_database(unsigned short number, bool from_storage);

/**
 * @brief Call database wrapper on specific numbers and write from database to file or backwards.
 *
 * number = 0 and back = false is used to start the wrapper thread
 * number = 0 and back = true is used to stop the wrapper thread
 *
 * Use number of sector directly to write sector,
 * number|65536 to check existence of sector.
 * @param number sector number
 * @param back direction; true for adding file to storage
 * @param remove_file whether to remove the file from the storage after getting, important only in case of back=false
 * @return successfullness of the operation
 */
bool database_to_file_or_back(int number, bool back, bool remove_file=true);

/**
 * Check whether a sector is in the list of sectors without Lee external symmetries (lower).
 * @param test_sector sector number
 * @return check result
 */
bool in_lsectors(unsigned short test_sector);

/**
 * Read from stream when communicating via pipe
 * @param buf buffer to write to, can be reallocated if size is not enough
 * @param buf_size location of the size to be written to, can be allocated
 * @param stream_from_child stream to read from
 */
void read_from_stream(char **buf, int *buf_size, FILE *stream_from_child);

/**
 * Instruct operating system to stop caching the file at path.
 * @param path full path to the file.
 */
void file_not_needed(string path);

/**
 * Used in inverting the denominator of fraction when working in PRIME mode,
 * It is essential for a and b to be relative primes.
 * @param a first number
 * @param b second number
 * @return a/b modular common::prime
 */
int128_t mul_inv(int128_t a, int128_t b);

/**
 * @brief Bring number stored in s to the mod range by common::prime.
 * @param s input string
 * @param positive_number_in_range indication whether we know that the number in the string is positive and does not exceed common::prime
 * If positive_number_in_range is set to true, that means we are reading a POSITIVE number, and
 * hence if it is read as negative, we should add 2^64 to it (autmatically).
 * If positive_number_in_range is false, then the number is small, negative is negative, so we are adding common::prime.
 * @return resulting number in proper range
 */
unsigned long long mod(string &s, bool positive_number_in_range = false);


/**
 * Replace all occurrences of from string in str to to string.
 * @param str string that will be updated
 * @param from etalon of substring to be replaced
 * @param to replacement string
 * @return copy of updated string
 */
string ReplaceAll(string str, const string &from, const string &to);

/**
 * Substitute all variables from common::variable_replacements in str, using ReplaceAll().
 * @param str input string
 * @return copy of updated string
 */
string ReplaceAllVariables(string str);

/**
 * @brief Wrapper for coefficient in monoms.
 *
 * Is a number, if used in prime version of FIRE, string in normal version.
 */
class COEFF {
    public:
#if defined(PRIME)
    unsigned long long n; ///< Exact value of coefficient modulo selected prime, used in PRIME mode.
#elif defined(MPQ)
    mpq_class s; // MPQ mode
#elif defined(FlintX)
    flint::fmpz_poly_qxx s;
    static string x;
#elif defined(FMPQ)
    vector<fmpq_t> s;
#elif defined(FlintC)
    vector<fmpz_poly_q_t> s;
    static string x;
#elif defined(FlintM)
    vector<fmpz_mpoly_q_t> s;
    static fmpz_mpoly_ctx_t ctx;
    static const char** xs;
#else
    std::string s; ///< String representation of coefficient, used in normal mode.
#endif

    COEFF() {
#ifdef PRIME
        n = 0;
#elif defined(MPQ)
        s = 0;
#elif defined(FlintX)
        s.set_zero();
#elif defined(FMPQ)
        s = vector<fmpq_t>(1);  
        fmpq_init(s[0]); 
        fmpq_zero(s[0]);
#elif defined(FlintC)
        s = vector<fmpz_poly_q_t>(1);  
        fmpz_poly_q_init(s[0]); 
        fmpz_poly_q_zero(s[0]);
#elif defined(FlintM)
        s = vector<fmpz_mpoly_q_t>(1);  
        fmpz_mpoly_q_init(s[0],ctx); 
        fmpz_mpoly_q_zero(s[0],ctx);
#else
        s = "";
#endif
    }

#if defined(FMPQ)
    COEFF(const COEFF & c) {
        s = vector<fmpq_t>(1);  
        fmpq_init(s[0]);
        fmpq_set(s[0],c.s[0]);
    }
    COEFF & operator=(const COEFF & c) {
        s = vector<fmpq_t>(1);  
        fmpq_init(s[0]);
        fmpq_set(s[0],c.s[0]);
        return *this;
    }
    COEFF(COEFF && c) : s(std::move(c.s)) { }
    COEFF & operator=(COEFF && c) {
        s = std::move(c.s);
        return *this;
    }
    ~COEFF() {
        if(s.size()>0) fmpq_clear(s[0]);
    }
#elif defined(FlintC)
    COEFF(const COEFF & c) {
        s = vector<fmpz_poly_q_t>(1);  
        fmpz_poly_q_init(s[0]);
        fmpz_poly_q_set(s[0],c.s[0]);
    }
    COEFF & operator=(const COEFF & c) {
        s = vector<fmpz_poly_q_t>(1);  
        fmpz_poly_q_init(s[0]);
        fmpz_poly_q_set(s[0],c.s[0]);
        return *this;
    }
    COEFF(COEFF && c) : s(std::move(c.s)) { }
    COEFF & operator=(COEFF && c) {
        s = std::move(c.s);
        return *this;
    }
    ~COEFF() {
        if(s.size()>0) fmpz_poly_q_clear(s[0]);
    }
#elif defined(FlintM)
    COEFF(const COEFF & c) {
        s = vector<fmpz_mpoly_q_t>(1);  
        fmpz_mpoly_q_init(s[0],ctx);
        fmpz_mpoly_q_set(s[0],c.s[0],ctx);
    }
    COEFF & operator=(const COEFF & c) {
        s = vector<fmpz_mpoly_q_t>(1);  
        fmpz_mpoly_q_init(s[0],ctx);
        fmpz_mpoly_q_set(s[0],c.s[0],ctx);
        return *this;
    }
    COEFF(COEFF && c) : s(std::move(c.s)) { }
    COEFF & operator=(COEFF && c) {
        s = std::move(c.s);
        return *this;
    }
    ~COEFF() {
        if(s.size()>0) fmpz_mpoly_q_clear(s[0],ctx);
    }
#endif

    /**
     * Constructor from a number, now supports 1 or -1
     * @param number a number
     * @return mew coefficient
     */
    COEFF(int64_t number) {
#ifdef PRIME
        if (number == 1) n = 1;
        else if (number == -1) n = common::prime - 1;
        else abort();
#elif defined(MPQ)
        if (number == 1) s = 1;
        else if (number == -1) s = -1;
        else abort();
#elif defined(FlintX)
    if(number==1) s=1;
    else if(number==-1) s=-1;
    else abort();
#elif defined(FMPQ)
    s = vector<fmpq_t>(1);  
    fmpq_init(s[0]);
    if(number==1 || number==-1) fmpq_set_si(s[0],number);
    else abort();
#elif defined(FlintC)
    s = vector<fmpz_poly_q_t>(1);  
    fmpz_poly_q_init(s[0]);
    if(number==1 || number==-1) fmpz_poly_q_set_si(s[0],number);
    else abort();
#elif defined(FlintM)
    s = vector<fmpz_mpoly_q_t>(1);  
    fmpz_mpoly_q_init(s[0],ctx);
    if(number==1 || number==-1) fmpz_mpoly_q_set_si(s[0],number,ctx);
    else abort();
#else
        if (number == 1)
            s = "1";
        else if (number == -1)
            s = "-1";
        else
            abort();
#endif
    }

    /**
    * Check whether the coefficient is zero or non-set, universal for both versions of FIRE
    * @return result of the check
    */
    bool empty() const {
        #ifdef PRIME
            return (!n);
        #elif defined(MPQ)
            return s==0;
        #elif defined(FlintX)
            return s.is_zero();
        #elif defined(FMPQ)
            return fmpq_is_zero(s[0]);
        #elif defined(FlintC)
            return fmpz_poly_q_is_zero(s[0]);
        #elif defined(FlintM)
            return fmpz_mpoly_q_is_zero(s[0],ctx);
        #else
            return ((s == "") || (s == "0"));
        #endif
    }

    /** @fn operator==(const COEFF &, const COEFF &)
     * @brief Checks whether two coefficients are equal
     * @param c1 first coefficient
     * @param c2 second coefficient
     * @return Check result
     */
    friend bool operator==(const COEFF &c1, const COEFF &c2) {
        //COEFF result;
    #ifdef PRIME
        return (c1.n == c2.n);
    #elif defined(FMPQ)
        return fmpq_equal(c1.s[0],c2.s[0]);
    #elif defined(FlintC)
        return fmpz_poly_q_equal(c1.s[0],c2.s[0]);
    #elif defined(FlintM)
        return fmpz_mpoly_q_equal(c1.s[0],c2.s[0],ctx);
    #else
        return (c1.s == c2.s);
    #endif
    }

    /** @fn operator+(const COEFF &, const COEFF &)
     * @brief Add two coefficients
     * @param c1 first coefficient
     * @param c2 second coefficient
     * @return Their sum
     */
    friend COEFF operator+(const COEFF &c1, const COEFF &c2) {
        COEFF result;
#ifdef PRIME
        uint128_t n1 = c1.n;
        uint128_t n2 = c2.n;
        uint128_t n_new = n1 + n2;
        n_new = n_new % common::prime;
        result.n = n_new;
#elif defined(MPQ)
	if (c1.s == 0) result.s = c2.s;
        else if (c2.s == 0) result.s = c1.s;
        else result.s = c1.s + c2.s;
#elif defined(FlintX)
    result.s = (c1.s+c2.s).evaluate();
#elif defined(FMPQ)
    fmpq_add(result.s[0],c1.s[0],c2.s[0]);
#elif defined(FlintC)
    fmpz_poly_q_add(result.s[0],c1.s[0],c2.s[0]);
#elif defined(FlintM)
    fmpz_mpoly_q_add(result.s[0],c1.s[0],c2.s[0],ctx);
#else
        if (c1.s == "")
            result.s = c2.s;
        else if (c2.s == "")
            result.s = c1.s;
        else
            result.s = "(" + c1.s + ") + (" + c2.s + ")";
#endif
        return result;
    }

    /** @fn operator-(const COEFF &, const COEFF &)
     * @brief Substract two coefficients
     * @param c1 first coefficient
     * @param c2 second coefficient
     * @return Their difference
     */
    friend COEFF operator-(const COEFF &c1, const COEFF &c2) {
        COEFF result;
#ifdef PRIME
        uint128_t n1 = c1.n;
        uint128_t n2 = c2.n;
        uint128_t p = common::prime;
        n2 = p - n2;
        uint128_t n_new = n1 + n2;
        n_new = n_new % common::prime;
        result.n = n_new;
#elif defined(MPQ)
        if ((c1.s == 0) && (c2.s == 0)) result.s = 0;
        else if (c1.s == 0) result.s = -c2.s;
        else if (c2.s == 0) result.s = c1.s;
        else result.s = c1.s - c2.s;
#elif defined(FlintX)
        result.s = (c1.s-c2.s).evaluate();
#elif defined(FMPQ)
        fmpq_sub(result.s[0],c1.s[0],c2.s[0]);
#elif defined(FlintC)
        fmpz_poly_q_sub(result.s[0],c1.s[0],c2.s[0]);
#elif defined(FlintM)
        fmpz_mpoly_q_sub(result.s[0],c1.s[0],c2.s[0],ctx);
#else
        if ((c1.s == "") && (c2.s == ""))
            result.s = "";
        else if (c1.s == "")
            result.s = "- (" + c2.s + ")";
        else if (c2.s == "")
            result.s = c1.s;
        else
            result.s = "(" + c1.s + ") - (" + c2.s + ")";
#endif
        return result;
    }

    /** @fn operator*(const COEFF &, const COEFF &)
     * @brief Multiply two coefficients
     * @param c1 first coefficient
     * @param c2 second coefficient
     * @return Their product
     */
    friend COEFF operator*(const COEFF &c1, const COEFF &c2) {
        COEFF result;
#ifdef PRIME
        uint128_t n1 = c1.n;
        uint128_t n2 = c2.n;
        uint128_t n_new = n1 * n2;
        n_new = n_new % common::prime;
        result.n = n_new;
#elif defined(MPQ)
        if (c1.s == 0) result.s = 0;
        else if (c2.s == 0) result.s = 0;
        else result.s = c1.s * c2.s;
#elif defined(FlintX)
        result.s = (c1.s * c2.s).evaluate();
#elif defined(FMPQ)
        fmpq_mul(result.s[0], c1.s[0], c2.s[0]);
#elif defined(FlintC)
        fmpz_poly_q_mul(result.s[0], c1.s[0], c2.s[0]);
#elif defined(FlintM)
        fmpz_mpoly_q_mul(result.s[0], c1.s[0], c2.s[0], ctx);
#else
        if (c1.s == "")
            result.s = "";
        else if (c2.s == "")
            result.s = "";
        else
            result.s = "(" + c1.s + ") * (" + c2.s + ")";
#endif
        return result;
    }
};

#if defined(FlintC)

inline string fmpz_poly_q_get_str2(const fmpz_poly_q_struct* p) {
    auto cs = fmpz_poly_get_str(fmpz_poly_q_numref(p));
    string s(cs);
    flint_free(cs);
    if(!fmpz_poly_is_one(fmpz_poly_q_denref(p))) {
        cs = fmpz_poly_get_str(fmpz_poly_q_denref(p));
        string ds(cs);
        flint_free(cs);
        s = s+"/"+ds;
    }
    return s;
}

inline void fmpz_poly_q_set_str2(fmpz_poly_q_struct* p, char* pos) {
    char *end = pos;
    while (*end != '/' && *end != '\0') ++end;
    if('/'==*end) { /* numerator / denominator */
        *end = '\0';
        fmpz_poly_set_str(fmpz_poly_q_numref(p),pos);
        pos = end+1;
        fmpz_poly_set_str(fmpz_poly_q_denref(p),pos);
    } else {
        fmpz_poly_set_str(fmpz_poly_q_numref(p),pos);
        fmpz_poly_set_si(fmpz_poly_q_denref(p),1);
    }
}

inline void fmpz_poly_q_set_str2(fmpz_poly_q_struct* p, const char* pos) {
    int n = strlen(pos);
    char buff[n+1];
    strcpy(buff,pos);
    fmpz_poly_q_set_str2(p, buff);
}

#elif defined(FlintM)

inline string fmpz_mpoly_q_get_str(const fmpz_mpoly_q_struct* mp) {
    auto cs = fmpz_mpoly_get_str_pretty(fmpz_mpoly_q_numref(mp),COEFF::xs,COEFF::ctx);
    string s(cs);
    flint_free(cs);
    if(!fmpz_mpoly_is_one(fmpz_mpoly_q_denref(mp),COEFF::ctx)) {
        cs = fmpz_mpoly_get_str_pretty(fmpz_mpoly_q_denref(mp),COEFF::xs,COEFF::ctx);
        string ds(cs);
        flint_free(cs);
        bool is_nc = fmpz_mpoly_is_fmpz(fmpz_mpoly_q_numref(mp),COEFF::ctx);
        bool is_dc = fmpz_mpoly_is_fmpz(fmpz_mpoly_q_denref(mp),COEFF::ctx);
        if(!is_nc) s = "("+s+")";
        if(!is_dc) ds = "("+ds+")";
        s = s+"/"+ds;
    }
    return s;
}

inline void fmpz_mpoly_q_set_str(fmpz_mpoly_q_struct* mp, char* pos) {
    char *end = pos;
    while (*end != '/' && *end != '\0') ++end;
    if('/'==*end) { /* numerator / denominator */
        *end = '\0';
        fmpz_mpoly_set_str_pretty(fmpz_mpoly_q_numref(mp),pos,COEFF::xs,COEFF::ctx);
        pos = end+1;
        fmpz_mpoly_set_str_pretty(fmpz_mpoly_q_denref(mp),pos,COEFF::xs,COEFF::ctx);
    } else {
        fmpz_mpoly_set_str_pretty(fmpz_mpoly_q_numref(mp),pos,COEFF::xs,COEFF::ctx);
        fmpz_mpoly_set_si(fmpz_mpoly_q_denref(mp),1,COEFF::ctx);
    }
}

inline void fmpz_mpoly_q_set_str(fmpz_mpoly_q_struct* mp, const char* pos) {
    int n = strlen(pos);
    char buff[n+1];
    strcpy(buff,pos);
    fmpz_mpoly_q_set_str(mp, buff);
}

inline void fmpz_mpoly_q_addmul(fmpz_mpoly_q_struct* rop, const fmpz_mpoly_q_struct* op1, const fmpz_mpoly_q_struct* op2, const fmpz_mpoly_ctx_t ctx) {
    fmpz_mpoly_q_t mp;
    fmpz_mpoly_q_init(mp,ctx);
    fmpz_mpoly_q_mul(mp,op1,op2,ctx);
    fmpz_mpoly_q_add(rop,rop,mp,ctx);
    fmpz_mpoly_q_clear(mp,ctx);
}
    
#endif


#ifdef WITH_ZSTD

#include <zstd.h>
/**
 * @brief ZStandard Compressor extension for the kyotocabinet.
 */
class ZstdCompressor : public kyotocabinet::Compressor {
private:
    char *compress(const void *buf, size_t size, size_t *sp) {
        _assert_(buf && size <= MEMMAXSIZ && sp);
        *sp = ZSTD_compressBound(size);
        char *out = new char[*sp + 1];
        if (common::lthreads_number == 1) {
            *sp = ZSTD_compressCCtx(ccontext, out, *sp, buf, size, common::compressor_level);
        } else {
            *sp = ZSTD_compress(out, *sp, buf, size, common::compressor_level);
        }
        out[*sp] = '\0';
        return out;
    }

    ZstdCompressor(const ZstdCompressor&) = delete;
    void operator=(ZstdCompressor const &x) = delete;

    char *decompress(const void *buf, size_t size, size_t *sp) {
        _assert_(buf && size <= MEMMAXSIZ && sp);
        if (!sp) {
            cout<<"Nullptr passed to decompress"<<endl;
            abort();
        }
        *sp = ZSTD_getDecompressedSize(buf, size);
        if (!*sp) {
            return nullptr;
        }
        char *out = new char[*sp + 1];
        if (common::lthreads_number == 1) {
            if (!ZSTD_decompressDCtx(dcontext, out, *sp, buf, size)) {
                delete[] out;
                return nullptr;
            }
        } else {
            if (!ZSTD_decompress(out, *sp, buf, size)) {
                delete[] out;
                return nullptr;
            }
        }
        out[*sp] = '\0';
        return out;
    }
    public:
    ZstdCompressor() {
        ccontext = ZSTD_createCCtx();
        dcontext = ZSTD_createDCtx();
    }
    ~ZstdCompressor() {
        ZSTD_freeCCtx(ccontext);
        ZSTD_freeDCtx(dcontext);
    }

    private:
        ZSTD_CCtx* ccontext;
        ZSTD_DCtx* dcontext;
};

#endif


#endif // COMMON_H_INCLUDED
