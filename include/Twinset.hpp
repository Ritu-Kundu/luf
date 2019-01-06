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

/** Defines the class Twinset.
 *
 */

#ifndef TWINSET_HPP
#define TWINSET_HPP

//#define DETAILS

#include <cmath>
#include <cstdlib>

//#include <math.h>
//#include <numeric>

#include "globalDefs.hpp"
// TODO : Making deletion of twin-set at runtime to make better use of space
namespace luf {

// deletion from the twin-sets are virtual
class Twinset {

  /** Block in a twin-set.
  * They are added in the twin-set in decreasing order of positions in the text.
  * In other words, blocks of the twin-sets correspond to the decreasing order
  * of the text-positions.
  * Immediately following block will represent a position on the left.
  *  */
  using Block = struct blo {
    /** Position in the text at which this block starts */
    UINT position;

    /** Pointer to the session in which this block lies. */
    SessionIter session;

    /** Indicates if this block has been deleted from this twinset. */
    bool is_delete;

    /** Indicates the previous reference for which this position was queried for
      * beta.
      * Initialised to 0
      * Updated while computing by FindBeta
      */
    UINT recent_querying_ref;

    /** Indicates index of the nearest tick of the previous reference for
   * which this block was a failure
   * Initialised to itself
   * Updated while computing by FindBeta
   * */
    UINT nearest_tick_ind;
  };

public:
  Twinset(const UINT len, const UINT ref, const UINT level, const UINT luf);

  //~Twinset(void);

  void reset(const UINT ref, const SessionIter ref_session,
             UINT &short_sensible);
  bool find_position_first_in_session(UINT &query_pos, UINT &num_blocks);
  bool find_position_quick(UINT &query_pos);
  bool find_position_after_failure(UINT &query_pos, UINT &num_blocks);
  bool find_position_long_query(UINT dist, UINT &query_pos);
  void handle_outside_long_query(const SessionIter sit, const UINT q_pos);
  SessionIter insert_session(std::list<UINT> &positions);
  bool is_empty() const;
  bool is_first_in_session() const;
  bool is_failed_ref() const;
  UINT creator_ref() const;
  UINT level() const;
  void print() const;
  void print_session_succ_list() const;
  //////////////////////// private ////////////////////////
private:
  /** Length of the blocks in this twin-set */
  const UINT _length;

  /** Reference that created this twin-set */
  const UINT _ref;

  /** LUF of the Reference that created this twin-set */
  const UINT _ref_luf;

  /** Level in the stack-tree of this twin-set */
  const UINT _level;

  /** Sessions of this twin-set
   * Deletion of sessions is necessary in order to avoid keep testing deleted
   * sessions
  */
  std::list<Session> _sessions;

  /** Blocks of this twin-set */
  std::vector<Block> _blocks;

  //////////////////////////////////////////////////////////
  // Data-structures needed for maintainance of the twinset
  // afer operations by FindBeta
  //////////////////////////////////////////////////////////

  /** Index of the Reference from the twin-set which is trying to find beta */
  UINT _querying_ref_ind;

  /** Index of the Reference in its session
   * Set on the reset call
   */
  UINT _querying_ref_session_ind;

  /** Indicates if the reference is the first in its session
   * Set on the reset call
   */
  bool _is_first_in_session;

  

  /** List of indices of the blocks-vector which resulted in successful queries
   * for the previous querying reference.
   * The first ref in this session clears old and fills the list for the first time.
   * Following references update it, potentially shortening the list.
   * Will be used by the next reference within the same session.
   * p_query_pos differentiates between query-positions of the current ref (i.e.
   * updated ones) and the previous ones.
   * From beginning to just before the p_query_pos are the updated ones.
   * Remaining (from p_query_pos to end) are from the previous reference. These
   * will be cleared when the next reference will call reset.
   * current_ref keeps changing the indices up to position pointedd to by
   * p_query_pos.
   *
   * Reset call fto a reference
   *   - deletes everytihng after p_query_pos
   *   - reset p_query_pos to the beginning
  */
  std::list<UINT> _success_queries_pvs;

  /** Pointer in the list success_queries_pvs to show which position in the list
   * was recently used to calculate query-position
   *  If this fails, everything is deleted after (including this one) by the
   * next reference from the same session in its reset call.
   * If it passes, the index is changed accordingly and pointer shifts to the
   * next one.
   * Rset on the reset call
   */
  std::list<UINT>::iterator _p_query_pos;

  /** Distance of the current querying ref from the previous querying ref of the
   * session
   * Initially 0.
   * Keeps changing with different references making calls to FindBeta
   */
  UINT _distance;

  /** List of indices of consecutive failures. 
   * Indices added when the position to try is queried by FindLCP.
   * recent_querying_ref of the block is set before adding its index in the list
   * List is reset when call to reset is made
   * Indices are dealt in accordance if the last poition added was a success or failure.
   * These indices are used to set corresponding nearest_tick_offset for failures if a success is found.
   * The success is deleted (virtualy).
   * Next call to find_position will mean  that the last index added was a success.
   * The first index in this list is always that of the querying reference // Used to avoid testing boundary condition
  */ std::list<UINT>
      _queried_indices;

  /** Iterator to the session where we started finding beta
   * To be used to answer long queries.
   * Set by iniial call; later updated by every success
   */
  SessionIter _beta_query_session;
};

} // end namespace
#endif
