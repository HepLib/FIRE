/** @file point.h
 *  @author Alexander Smirnov
 *
 *  This file is a part of the FIRE package
 *  It contains the basic point class corresponding to a Feynman integral
 */


#ifndef _points_h_
#define _points_h_

#include "common.h"

/** Used for enumerating virtual points (masked expressions) */
typedef uint32_t virt_t;

#include <iomanip>

#ifndef DOXYGEN_DOCUMENTATION
class point_fast;
#endif

/** Used for storing ibps without substituted indices */
typedef vector<pair<vector<COEFF>, point_fast > > ibp_type;

/**
 * @brief Ordered point class.
 *
 * Contains point's degrees already multiplied by ordering matrix and calculated sector number.
 * There also are virtual points - there are not real integrals, but some masking expressions of tails in lower sectors.
 * Point is designed to fit into the 24 = 3*8 byte size
 */
class point {
public:
    /** Globally known preferred points. */
    static vector<set<vector<t_index> > > preferred;
    /** Globally known preferred fast-points. */
    static vector<set<point_fast> > preferred_fast;
    /** Flag for Mathematical point printing in FTool module. */
    static bool print_g;

    /**
     * Default constructor, fills with zeros.
     */
    point();

    /**
     * Constructor from another point.
     * @param p another point
     * @return new point
     */
    point(const point &p);

    /**
     * @brief Constructor from point_fast.
     * Constructs only real points, no sector 1 or masking
     * @param pf the fast version of the point
     * @param ssector sector (not number) if we know it, -1 otherwise.
     */
    point(const point_fast &pf, SECTOR ssector = -1);

    /**
     * Constructor from point and added shift point_fast in case sector is not changed.
     * We must know and specify ssector.
     * @param p initial point
     * @param v shift to be added
     * @param ssector sector number
     * @param not_preferred whether the new point is not preferred
     */
    point(const point &p, const point_fast &v, SECTOR ssector, bool not_preferred);

    /** @fn operator=(const point &) const
     * @brief Set one point equal to enother
     * @return The new point
     */
    point &operator=(const point &);

    /**
     * Constructor from vector of indices in case sector is not changed.
     * @param v the vector of indices
     * @param virt the masking number, normal integral if zero
     * @param ssector sector number if we know it, -1 otherwise.
     */
    point(const vector<t_index> &v, virt_t virt = 0, SECTOR ssector = -1);

    /** Constructor from safe buffer.
     *  @param buf the buffer filled by safe_string().
     * */
    point(char *buf);

    /** For normal points it stores the point degrees already multiplied by the ordering matrix.
    *  Bytes are stored in inverted manner, to speedup comparison.
    *
    *  If point is virtual, then the ww[MAX_IND - 1] is equal to zero and 4 bytes before that are used to store
    *  virtual number point.
    *  Last two bytes are the so called h1.
     */
    #ifndef SMALL_POINT
    char ww[24];
    #else
    char ww[16];
    #endif

    /**
    * 15 top bits are for the sector number.
    * 1 bit for the "shift" value at the moment of point initializations, hence some points are lower (virtual among them).
    * Virtual points are not real integrals, but some masking expressions of tails in lower sectors.
     * @return two last bytes as an unsigned short.
     */

    unsigned short h1() const {
        return *reinterpret_cast<const unsigned short*>(ww + sizeof(point) - 2);
    };

    /**
     * @return pointer to the last bytes, referred to as h1 (see h1() description)
     */
    unsigned short* h1p() {
        return reinterpret_cast<unsigned short*>(ww + sizeof(point) - 2);
    };

    /**
    * Create a string that is written in tables
    * @return a string to be written in tables
     */
    string number() const;

    /**
    * Retrieve the sector number from h1.
    * @return sector number of the point.
     */
    unsigned short s_number() const { return (h1() >> 1); }

    /**
    * Check if the point is virtual.
    * @return true if point is virtual, false otherwise.
     */
    bool virt() const { return (ww[sizeof(point) - 3] == 0); }

    /**
    * Get the level of point - number of non-negative indices
    * @return level of point.
     */
    int level() const;


    /**
    * Fill buffer with information about without any special symbols that can't be piped.
    * @param buf pointer allocated with 48 bytes.
     */
    void safe_string(char *buf) const;

    /**
     * Get the real point indices, using inverse ordering.
     * @return vector of real point indices.
     */
    vector<t_index> get_vector() const;

    /**
     * Check if point data is filled with zeroes.
     * Point data can be filled with zeroes only as an exception, for example
     * if we where searching for a point and found none.
     * @return true if point data is filled with zeroes, false otherwise.
     */
    bool is_zero() const;

    /** @fn operator==(const point &, const point &) const
     * @brief Compare two points
     * @param p1 first point
     * @param p2 second point
     * @return True if they are equal
     */
    friend bool operator==(const point &p1, const point &p2) {
        #ifndef SMALL_POINT
            if (reinterpret_cast<const uint64_t *>(p1.ww)[2] ^ reinterpret_cast<const uint64_t *>(p2.ww)[2]) return false;
            if (reinterpret_cast<const uint64_t *>(p1.ww)[1] ^ reinterpret_cast<const uint64_t *>(p2.ww)[1]) return false;
            if (reinterpret_cast<const uint64_t *>(p1.ww)[0] ^ reinterpret_cast<const uint64_t *>(p2.ww)[0]) return false;
            return true;
        #else
            return !(reinterpret_cast<const uint128_t *>(p1.ww)[0] ^ reinterpret_cast<const uint128_t *>(p2.ww)[0]);
        #endif
    }

    /** @fn operator!=(const point &, const point &) const
     * @brief Compare two points
     * @param p1 first point
     * @param p2 second point
     * @return True if they are not equal
     */
    friend bool operator!=(const point &p1, const point &p2) {
        #ifndef SMALL_POINT
            return (!(p1 == p2));
        #else
            return reinterpret_cast<const uint128_t *>(p1.ww)[0] ^ reinterpret_cast<const uint128_t *>(p2.ww)[0];
        #endif
    }

    /** @fn operator<=(const point &, const point &) const
     * @brief Compare two points
     * @param p1 first point
     * @param p2 second point
     * @return True if they are equal or first is less
     */
    friend bool operator<=(const point &p1, const point &p2) {
        #ifndef SMALL_POINT
            return (p1 == p2 || p1 < p2);
        #else
            return reinterpret_cast<const uint128_t *>(p1.ww)[0] <= reinterpret_cast<const uint128_t *>(p2.ww)[0];
        #endif
    }

    /**
     * @fn operator<(const point &, const point &) const
     * Compare points first by h1 field (including virtuality), then by the ww vector.
     * @param p1 first point.
     * @param p2 second point.
     * It's a tricky compare:
     * vectors are interpreted by ints, so direction is backwards to agree with memory orderings.
     * h1 that used to be in compare stands in memory after buf[22], so plays the role of top 2 bytes.
     * MAX_IND = 22 is principally important.
     * @param p1 first point
     * @param p2 second point
     * @return true if p1 < p2, false otherwise
     */
    friend bool operator<(const point &p1, const point &p2) {
        #ifndef SMALL_POINT
        if (reinterpret_cast<const uint64_t *>(p1.ww)[2] ^ reinterpret_cast<const uint64_t *>(p2.ww)[2]) {
            return reinterpret_cast<const uint64_t *>(p1.ww)[2] < reinterpret_cast<const uint64_t *>(p2.ww)[2];
        }
        if (reinterpret_cast<const uint64_t *>(p1.ww)[1] ^ reinterpret_cast<const uint64_t *>(p2.ww)[1]) {
            return reinterpret_cast<const uint64_t *>(p1.ww)[1] < reinterpret_cast<const uint64_t *>(p2.ww)[1];
        }
        return reinterpret_cast<const uint64_t *>(p1.ww)[0] < reinterpret_cast<const uint64_t *>(p2.ww)[0];
        #else
            return reinterpret_cast<const uint128_t *>(p1.ww)[0] < reinterpret_cast<const uint128_t *>(p2.ww)[0];
        #endif
    }


    friend ostream &operator<<(ostream &out, const point &p);

    /** Globally known Lee internal symmetries for points. */
    static map<unsigned short,
            vector< // because there can be multiple
                    pair<vector<t_index>,
                            vector<
                                    pair<
                                            vector<
                                                    pair<COEFF, point_fast>
                                            >,
                                            t_index>
                            >
                    >
            >
    > ibases;

    /** Globally known Lee external symmetries for points. */
    static map<unsigned short,
            pair<vector<t_index>,
                    vector<
                            pair<
                                    vector<
                                            pair<COEFF, point_fast>
                                    >,
                                    t_index>
                    >
            >
    > dbases;

    /**
     * Integration by part relations (without substituting indices). Each ibp is a vector of pairs (split coefficient and point shift).
     */
    static vector<ibp_type> ibps;
}
#ifndef SMALL_POINT
__attribute__((aligned (8))); // size is 24, so align by 8
#else
__attribute__((aligned (16))); // size is 16, let's alighn by 16 too
#endif

/** @brief Primitive point class.
 *
 *  Contains point's indices as is.
 */
class point_fast {
public:
    ~point_fast();

    point_fast();

    /**
     * Copy constructor
     * @param v initial point
     */
    point_fast(const point_fast &v);

    /**
     * @fn point_fast(const point &)
     * Constructor from a full point class
     * @param p point
     * @brief Uses the inverse ordering matrix to get a set if indices from the point class
     */
    point_fast(const point &p);

    /**
     * Constructor from vector of indices
     * @param v vector of indices
     */
    point_fast(const vector<t_index> &v);

    /** @fn operator=(const point_fast &) const
     * @brief Make two points equal
     * @return The new point
     */
    point_fast &operator=(const point_fast &);

    /** @fn operator==(const point_fast &, const point_fast &) const
     * @brief Compare two points
     * @param p1 first point
     * @param p2 second point
     * @return True if they are equal
     */
    friend bool operator==(const point_fast &p1, const point_fast &p2) {
        return !memcmp(p1.buf, p2.buf, MAX_IND);
    }

    /** @fn operator!=(const point_fast &, const point_fast &) const
     * @brief Compare two points
     * @param p1 first point
     * @param p2 second point
     * @return True if they are not equal
     */
    friend bool operator!=(const point_fast &p1, const point_fast &p2) {
        return memcmp(p1.buf, p2.buf, MAX_IND);
    }

    /** @fn operator<=(const point_fast &, const point_fast &) const
     * @brief Compare two points
     * @param p1 first point
     * @param p2 second point
     * @return True if they are equal or first is less
     */
    friend bool operator<=(const point_fast &p1, const point_fast &p2) {
        return (memcmp(p1.buf, p2.buf, MAX_IND) <= 0);
    }

    /** @fn operator<(const point_fast &, const point_fast &) const
     * @brief Compare two points
     * @param p1 first point
     * @param p2 second point
     * @return True if first is less. This is just an alphabetic comparisson having no relation with the ordering
     */
    friend bool operator<(const point_fast &p1, const point_fast &p2) {
        return (memcmp(p1.buf, p2.buf, MAX_IND) < 0);
    }

    /** @fn operator+(const point_fast &, const point_fast &) const
     * @brief Add points as vectors. Used in IBP generation.
     * @param p1 first point
     * @param p2 second point
     * @return Resulting point
     */
    friend point_fast operator+(const point_fast &p1, const point_fast &p2) {
        point_fast result;
        const t_index *pos1 = p1.buf;
        const t_index *pos2 = p2.buf;
        for (unsigned short i = 0; i != MAX_IND; ++i, ++pos1, ++pos2) result.buf[i] = *pos1 + *pos2;
        return result;
    }

    /**
     * The degree of the point
     * @return - the absolute values of the shift from the corresponding sector corner
     */
    point_fast degree() const;

    /**
     * Sector corresponding to a point
     * @return the corresponding sector (1 and -1 indices)
     */
    SECTOR sector_fast() const;

    /** Buffer for point's indices. */
    t_index buf[MAX_IND];
};

/**
 * Compare fast points l and r degree wise
 * This is a partial sorting function, for example, points in different sectors of same level cannot be compared
 * @param l first point
 * @param r second point
 * @return true if l is above r or equal, meaning all indices are greater or equal pairwise
 */
bool over_fast(const point_fast &l, const point_fast &r);

/**
 * Get set of points on a level
 * @param s starting point
 * @param pos number of dots
 * @param neg number of numerators
 * @return set of primitive points on a level
 */
set<point_fast> level_points_fast(point_fast s, const unsigned int pos, const unsigned int neg);

#endif
