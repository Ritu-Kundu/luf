/**
   luf: Tool to compute Longest Unbordered Factor (LUF) Array
  Copyright (C) 2017 Ritu Kundu
   This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

/** Implements class Twinset
 */
#include "../include/Twinset.hpp"

namespace luf {
Twinset::Twinset(const UINT len, const UINT ref, const UINT level,
                 const UINT luf)
    : _length(len), _ref(ref), _level(level), _ref_luf(luf) {}

// Twinset::~Twinset(void) {}

// returns total valid blocks (for short query)
// Assumes reference is in the Twinset
void Twinset::reset(const UINT ref, const SessionIter ref_session,
                    UINT &short_sensible) {
  /* delete all previous sessions */
  // These will never be used
  auto sit = _sessions.begin();
  while (sit != ref_session) {
    sit = _sessions.erase(sit);
  }
  /* Delete (virtually) all blocks before the ref block */
  // They will never be used
  sit = ref_session;
  UINT current_ind = sit->first_live_ind;
  Block &current_block = _blocks[current_ind];

  while (current_block.position != ref) {
    sit->num_blocks_alive -= 1;
    ++current_ind;
    current_block = _blocks[current_ind];
    sit->first_live_ind = current_ind;
  }

  // now sit points to the session it represents
  auto ref_ind = current_ind;
  auto &ref_block = current_block;

  // blocks in the same session will be skipped
  // If first in the session, calculate distance from the previous ref
  // If none, then its rank
  // If not the first one in the session, then only one query allowed
  if (ref_session->is_first_ref) { // first ref in this session

    ref_session->is_first_ref = false;
    _is_first_in_session = true;
    // +1 because it must check the partial prefix of the tick we consider
    // long query contains the whole word of the tick
    if (ref_block.recent_querying_ref ==
        0) { // this is the first one in this twin-set
      short_sensible = ref_ind + 1;
    } else {
      short_sensible = ref_ind - ref_block.recent_querying_ref + 1;
    }
    _success_queries_pvs.clear();
    _p_query_pos = _success_queries_pvs.begin();
    _queried_indices.clear();
    // insert ref as success-index
    _queried_indices.push_back(ref_ind);
  } else { // there has already been a reference in this session
    short_sensible = 1;
    _is_first_in_session = false;

    _distance = ref_ind - _querying_ref_ind;
    // Erase element from the previous to previos ref's successes
    _success_queries_pvs.erase(_p_query_pos, _success_queries_pvs.end());
    _p_query_pos = _success_queries_pvs.begin();
  }
  /* Reset the data-structures used for maintainance */
  _querying_ref_ind = ref_ind;
  _querying_ref_session_ind = ref_ind - ref_session->first_ind;
}

// returns the first position to query (if any) and number of blocks
// including & between hook of the last success and this position
bool Twinset::find_position_first_in_session(UINT &query_pos,
                                             UINT &num_blocks) {
  /* Do house-keeping */
  // delete success
  auto tick = _queried_indices.back();
  _queried_indices.pop_back();
  auto &tick_block = _blocks[tick];
  tick_block.is_delete = true;

  // only failures left in the list (if any)
  auto it = _queried_indices.begin();
  while (it != _queried_indices.end()) {
    _blocks[*it].nearest_tick_ind = tick;    
    it = _queried_indices.erase(it);
  }

  // Insert in the list of the success_query
  _success_queries_pvs.push_back(tick);

#ifdef DETAILS
  std::cout << "Finding position after success : deleted-pos : "
            << tick_block.position << std::endl;
#endif

  /* Find the session of this success */
  auto sit = tick_block.session;
  sit->num_blocks_alive -= 1;
#ifdef DETAILS
  std::cout << "Left in session: " << sit->num_blocks_alive << std::endl;
#endif

  if (sit->num_blocks_alive == 0) { // erase this session as well
    sit = _sessions.erase(sit);
#ifdef DETAILS
    std::cout << "This session deleted \n ";
#endif
  } else {
    // update if this was the live beginning of the session
    if (sit->first_live_ind == tick) {
      sit->first_live_ind = tick + 1;
    }
    // Go to the next session.
    ++sit;
  }
  /* Find position */
  _beta_query_session = sit;
  bool found = false;
  while (sit != _sessions.end()) { // session exists
    auto current_ind = sit->first_ind + _querying_ref_session_ind;
    if (sit->size > _querying_ref_session_ind &&
        !_blocks[current_ind].is_delete) { // block exist and is live
      found = true;
      query_pos = _blocks[current_ind].position;
      num_blocks = current_ind - _beta_query_session->first_ind + 1;

      _queried_indices.push_back(current_ind);
      _blocks[current_ind].recent_querying_ref = _querying_ref_ind;
#ifdef DETAILS
      std::cout << "Found: pos, num-blocks : " << query_pos << " " << num_blocks
                << std::endl;
#endif
      break;
    }
    ++sit;
  }
  return found;
}

bool Twinset::find_position_quick(UINT &query_pos) {
  /* Do house-keeping */

  auto &rit = _p_query_pos;
  auto &li = _success_queries_pvs;
  // pointed value is success (initially previous reference index)
  *rit += _distance;
  auto tick = *rit;
  auto &tick_block = _blocks[tick];
  tick_block.is_delete = true;

#ifdef DETAILS
  std::cout << "Finding quick position after success (not the first reference) "
               ": updated-pos : "
            << tick_block.position << std::endl;
#endif

  /* Find the session of this success */
  auto sit = tick_block.session;
  sit->num_blocks_alive -= 1;
#ifdef DETAILS
  std::cout << "Left in session: " << sit->num_blocks_alive << std::endl;
#endif

  if (sit->num_blocks_alive == 0) { // erase this session as well
    sit = _sessions.erase(sit);
#ifdef DETAILS
    std::cout << "This session deleted \n ";
#endif
  } else {
    // update if this was the live beginning of the session
    if (sit->first_live_ind == tick) {
      sit->first_live_ind = tick + 1;
    }
  }
  /* Find position */
  bool found = false;
  // Go to next position in the list
  ++rit;
  while (rit != li.end() && *rit + _distance < _blocks.size()) {
    auto &pvs_block = _blocks[*rit];
    auto current_ind = *rit + _distance;
    auto &current_block = _blocks[current_ind];
    if (!current_block.is_delete &&
        pvs_block.session ==
            current_block.session) { // it exists in the same session
      found = true;
      query_pos = current_block.position;
#ifdef DETAILS
      std::cout << "Found: pos : " << query_pos << std::endl;
#endif
      break;
    } else { // this position did not give anything, delte it, move beyond
      rit = li.erase(rit);
    }
  }
  return found;
}

// returns the first position to query (if any) and adds the number of blocks
// between last poistion queried and this position
bool Twinset::find_position_after_failure(UINT &query_pos, UINT &num_blocks) {
  auto last_query_ind = _queried_indices.back();
  auto block = _blocks[last_query_ind];
  auto sit = block.session;
  ++sit;
#ifdef DETAILS
  std::cout << "Finding position after failure : " << std::endl;
#endif
  bool found = false;
  while (sit != _sessions.end()) { // session exists
    auto current_ind = sit->first_ind + _querying_ref_session_ind;
    if (sit->size > _querying_ref_session_ind &&
        !_blocks[current_ind].is_delete) { // block exist and is live
      found = true;
      query_pos = _blocks[current_ind].position;
      num_blocks = current_ind - _beta_query_session->first_ind + 1;

      _queried_indices.push_back(current_ind);
      _blocks[current_ind].recent_querying_ref = _querying_ref_ind;
#ifdef DETAILS
      std::cout << "Found: pos, num-blocks : " << query_pos << " " << num_blocks
                << std::endl;
#endif
      break;
    }
    ++sit;
  }
  return found;
}

// returns the first position to make long-query (if any)
// At distance d from the position of the first trial
bool Twinset::find_position_long_query(UINT dist, UINT &query_pos) {
  auto ind = _beta_query_session->first_live_ind;
  auto &block = _blocks[ind];
  auto next_ind = dist;
  // if this is the first reference, then add LUF to the ditance
  // Else use the nearest tick

  next_ind += block.nearest_tick_ind;

#ifdef DETAILS
  std::cout << "Finding long query pos: dist : session-start(of pos)   : "
            << dist << " " << block.position << std::endl;
#endif
  bool found = false;
  if (block.nearest_tick_ind != ind &&
      next_ind < _blocks.size()) { // thesrs is a tick and d distance from it is
                                   // a valid position
    found = true;
    auto current_ind = next_ind;
    auto &current_block = _blocks[current_ind];
    query_pos = current_block.position;
#ifdef DETAILS
    std::cout << "Found : nearest-tick  query-pos :  : "
              << _blocks[block.nearest_tick_ind].position << " " << query_pos
              << std::endl;
#endif

    // Add this index to the list of queried-indices
    _queried_indices.push_back(current_ind);
    current_block.recent_querying_ref = _querying_ref_ind;
  }
  return found;
}

// inserts the index of the long_query poistion in the list of queried positions
// So that the next call to find_position sets everything correctly
// last poistion queried and this position
// Assumes that the given position is indeed in the session
void Twinset::handle_outside_long_query(const SessionIter sit,
                                        const UINT q_pos) {
  auto index = sit->first_ind + _querying_ref_session_ind;
  ;
  auto &block = _blocks[index];
  // A futile check
  if (block.position != q_pos) {
#ifdef DETAILS
    std::cout
        << "Handling outside: SOMEthing is wrong !!!! : query_pos, actual pos: "
        << q_pos << " " << block.position << std::endl;
#endif
  } else {
#ifdef DETAILS
    std::cout << "Handling outside: " << index << " " << block.position
              << std::endl;
#endif
  }
  _queried_indices.push_back(index);
}

// Insert a complete session and return its iterator
SessionIter Twinset::insert_session(std::list<UINT> &positions) {
  /* Create new session */
  Session s = {_blocks.size(), _blocks.size(), positions.size(),
               positions.size(), true};
  _sessions.push_back(s);
  auto sit = std::prev(_sessions.end());

  /* Add block for each position */
  for (auto p : positions) {
    Block b = {p, sit, false, 0, _blocks.size()};
    _blocks.push_back(b);
  }

  return sit;
}

bool Twinset::is_empty() const { return _sessions.empty(); }

bool Twinset::is_first_in_session() const { return _is_first_in_session; }

UINT Twinset::level() const { return _level; }

UINT Twinset::creator_ref() const { return _ref; }

bool Twinset::is_failed_ref() const {
  auto ref_block = _blocks[_querying_ref_ind];
  return !(ref_block.recent_querying_ref == 0);
}

// Prints only present elements
void Twinset::print() const {
  std::cout << "TWINSET : len: Creater-Ref : level: " << _length << " : "
            << _ref << " : " << _level << std::endl;
  std::cout << "{ ";
  auto sit = _sessions.begin();
  auto i = sit->first_live_ind;
  while (i < _blocks.size()) {
    auto &b = _blocks[i];
    if (sit != _sessions.end() &&
        i == sit->first_live_ind) { //  new session starts
      std::cout << " | ";
      ++sit;
    }
    if (!b.is_delete) {
      std::cout << b.position << " ";
    }
    i += 1;
  }
  std::cout << " } \n ";
}

// Prints the list of previous success position of the reference in the same
// session
void Twinset::print_session_succ_list() const {
  std::cout << "PREVIOUS REF LIST :  " << std::endl;
  std::cout << "[ ";
  for (auto i : _success_queries_pvs) {
    if (i == *_p_query_pos) {
      std::cout << "( ";
    }
    std::cout << _blocks[i].position << " ";
    if (i == *_p_query_pos) {
      std::cout << ") ";
    }
  }

  std::cout << " ] \n ";
}

} // end namespace
