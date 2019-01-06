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
void Twinset::reset(const UINT ref, UINT &short_num_blocks) {

  /* Find the session of this reference */
  auto sit = _sessions.begin();
  UINT current_ind = sit->first_ind;
  ;
  Block current_block = _blocks[current_ind];
 
  while (current_block.position != ref) {
    sit->first_ind =
        current_ind +
        current_block.next; // next of this becomes the first position
    sit->num_blocks_alive -= 1;
    if (sit->num_blocks_alive == 0) { // all deleted
      ++sit;                          // point to the next session
    }
    current_ind = sit->first_ind;
    current_block = _blocks[current_ind];    
  }

  // now sit points to the session it represents
  auto ref_session = sit;
  auto ref_ind = current_ind;
  auto ref_block = current_block;

  /* delete all previous sessions */
  // These will never be used
  sit = _sessions.begin();
  while (sit != ref_session) {
    sit = _sessions.erase(sit);
  }

  /* Reset the data-structures used for maintainance */
  _queried_indices.clear();
  _querying_ref_ind = ref_ind;

  // blocks in the same session will be skipped
  if (ref_block.recent_querying_ref ==
      0) { // this is the first one in this twin-set
    short_num_blocks = ref_ind + 1;
  } else {
    short_num_blocks = ref_ind - ref_block.recent_querying_ref;
  }

  // insert ref as success-index
  _queried_indices.push_back(ref_ind);
}

// returns the first position to query (if any) and number of blocks
// including&between
// hook of the last success and this position
bool Twinset::find_position(UINT &query_pos, UINT &num_blocks) {
  /* Do house-keeping */
  // update immediate neighbors of success
  auto tick = _queried_indices.back();
  _queried_indices.pop_back();
  auto tick_block = _blocks[tick];
  auto n = tick_block.next;
  auto p = tick_block.pvs;
  if (tick + n < _blocks.size()) { // there is a next
    _blocks[tick + n].pvs += p;
  }

  if (tick > p) { // there is a previous
    _blocks[tick - p].next += n;
  }

#ifdef DETAILS
  std::cout << "Finding position after success : deleted-pos : "
            << tick_block.position << std::endl;
#endif
  // only failures left in the list (if any)
  auto it = _queried_indices.begin();

  while (it != _queried_indices.end()) {
    _blocks[*it].nearest_tick_ind = tick;
    it = _queried_indices.erase(it);
  }

  /* Find the session of this success */
  auto sit = tick_block.session;
  sit->num_blocks_alive -= 1;
#ifdef DETAILS
  std::cout << "Left in session: " << sit->num_blocks_alive << std::endl;
#endif

  // blocks in the same session will be skipped
  auto skipped_blocks = sit->num_blocks_alive;
  if (sit->num_blocks_alive == 0) { // erase this session as well
    sit = _sessions.erase(sit);
#ifdef DETAILS
    std::cout << "This session deleted \n ";
#endif
  } else {
    // update if this was the beginning of the session
    if (sit->first_ind == tick) {      
      sit->first_ind = tick+n;
    }
    // Go to the next session, if it exists.
    ++sit;
  }
  /* Find position */
  bool found = false;
  if (sit != _sessions.end()) {
    found = true;
    // it's first block's position is the answer to the query
    auto current_ind = sit->first_ind;
    auto current_block = _blocks[current_ind];
    query_pos = current_block.position;
    num_blocks = current_ind - tick - skipped_blocks;
    // Add this index to the list of queried-indices
    _queried_indices.push_back(current_ind);
    current_block.recent_querying_ref = _querying_ref_ind;
    _beta_query_session = sit;
#ifdef DETAILS
    std::cout << "Found: pos, num-blocks : " << query_pos << " " << num_blocks
              << std::endl;
#endif
  }
  return found;
}

// returns the first position to query (if any) and adds the number of blocks
// between last poistion queried and this position
bool Twinset::find_position_after_failure(UINT &query_pos, UINT &num_blocks) {
  auto last_query_ind = _queried_indices.back();
  auto block = _blocks[last_query_ind];
  auto next_ind = last_query_ind + block.next;
#ifdef DETAILS
  std::cout << "Finding position after failure : " << std::endl;
#endif
  bool found = false;
  if (next_ind < _blocks.size()) { // valid position
    auto current_ind = next_ind;
    auto current_block = _blocks[current_ind];
    query_pos = current_block.position;
    num_blocks += block.next;

    // Add this index to the list of queried-indices
    _queried_indices.push_back(current_ind);
    current_block.recent_querying_ref = _querying_ref_ind;
    found = true;
#ifdef DETAILS
    std::cout << "Found: pos, num-blocks : " << query_pos << " " << num_blocks
              << std::endl;
#endif
  }
  return found;
}

// returns the first position to make long-query (if any)
// At distance d from the position of the first trial
bool Twinset::find_position_long_query(UINT dist, UINT &query_pos) {
  auto ind = _beta_query_session->first_ind;
  auto block = _blocks[ind];
  auto ref_block = _blocks[_querying_ref_ind];
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
    auto current_block = _blocks[current_ind];
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
  auto index = sit->first_ind;
  auto block = _blocks[index];
  while (block.position != q_pos) {
    index = index + block.next;
    block = _blocks[index];
  }
  std::cout << "Handling outside: " << index << " " << block.position
            << std::endl;
  _queried_indices.push_back(index);
}

// Insert a complete session and return its iterator
SessionIter Twinset::insert_session(std::list<UINT> &positions) {
  /* Create new session */
  Session s = {positions.size(), _blocks.size()};
  _sessions.push_back(s);
  auto sit = std::prev(_sessions.end());

  /* Add block for each position */
  for (auto p : positions) {
    Block b = {p, 1, 1, sit, _blocks.size(), 0};
    _blocks.push_back(b);
  }

  return sit;
}

bool Twinset::is_empty() const { return _sessions.empty(); }

UINT Twinset::level() const { return _level; }

UINT Twinset::creator_ref() const { return _ref; }

bool Twinset::is_failed_ref() const {
  auto ref_block = _blocks[_querying_ref_ind];
  return (ref_block.recent_querying_ref == 0);
}

// Prints only present elements
bool Twinset::print() const {
  std::cout << "TWINSET : len: Creater-Ref : level: " << _length << " : "
            << _ref << " : " << _level << std::endl;
  std::cout << "{ ";
  auto sit = _sessions.begin();
  auto i = sit->first_ind;
  while (i < _blocks.size()) {
    auto b = _blocks[i];
    if (sit != _sessions.end() && i == sit->first_ind) { //  new session starts
      std::cout << " | ";
      ++sit;
    }
    std::cout << b.position << " ";
    i += b.next;
  }
  std::cout << " } \n ";
}

} // end namespace
