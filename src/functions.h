#pragma once

#include "equation.h"

void apply_table(const pc_pair_ptr_vec &terms,
                 sector_count_t fixed_database_sector, unsigned short sector_level);

set<point, indirect_more>::reverse_iterator expressed_by(set<point, indirect_more> &to_test, sector_count_t sector_number);

void pass_back(const set<point, indirect_more> &cur_set, set<point, indirect_more>::const_reverse_iterator ritr,
          sector_count_t fixed_database_sector);

pc_pair_ptr_lst::iterator split(pc_pair_ptr_lst &terms, sector_count_t sector_number, uint64_t virts_number);
void mul_add_to(pc_pair_ptr_lst &terms1, const pc_pair_ptr_lst &terms2, const COEFF &coeff, bool skip_last);
void add_to(pc_pair_ptr_lst &terms1, const pc_pair_ptr_lst &terms2, bool skip_last);

equation apply(const vector<pair<vector<COEFF>, point_fast > >& ibp, point_fast v, const SECTOR ssector_fast);

point_fast lowest_in_sector_orbit_fast(const point_fast &p, SECTOR s, const vector<vector<vector<t_index> > > &sym);

void forward_stage(sector_count_t ssector_number);

void perform_substitution(sector_count_t ssector_number);

void Evaluate();

void work_with_master();

unsigned int sort_ibps(const point &p, const set<pair<unsigned int, unsigned int> > &current_levels,
                       vector<pair<point, pair<point_fast, unsigned short> > > &ibps_vector,
                       const vector<point_fast> &IBPdegree,
                       const vector<point_fast> &IBPdegreeFull);

void make_master(const point &);

void add_ibps(const COEFF& first_mul, const COEFF& second_mul, const ibp_type& first, const ibp_type& second, const SECTOR sector_fast, ibp_type& result);

int matching_index(const vector<COEFF>& first, const vector<COEFF>& second);

void mark_master_integrals(const point &Corner, const unsigned int pos, const unsigned int neg);

vector<pair<unsigned int, unsigned int> > under_levels(const unsigned int p0, const unsigned int m0);

void add_needed(map<sector_count_t, set<point> > &needed_lower, const point &p);

void clear_sector(sector_count_t sn, set<point, indirect_more> *ivpl);

bool sort_pair_point_coeff_by_point(const pair<point, COEFF> &lhs, const pair<point, COEFF> &rhs);

pc_pair_vec group_equal_in_sorted(const pc_pair_vec &mon);
