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

    /** offset until the next block present
     * Initially every block points to the immediate next (following).
     * As blocks will be deleted from the set (because of the successful lcp)
     * to be moved in another, it will keep changing.
     * It also tells how many blocks were initially present between the block
     * represented by the current block and that represented by its current
     * remaining neighbour .
     * Initialised with 1
    */
    UINT next;

    /** offset until the pvs block present
     * Initially every block points to the immediate previous (preceeding).
     * As blocks will be deleted from the set (because of the successful lcp)
     * to be moved in another, it will keep changing.
     * It also tells how many blocks were initially present between the block
     * represented by the current block and that represented by its current
     * remaining neighbour.
     * Initialised to 1.
     * First block of the twinset never needs previous. Thus, we do not need to
     * address boundary conditions
    */
    UINT pvs;

    /** Pointer to the session in which this block lies. */
    SessionIter session;

    /** Indicates index of the nearest tick of the previous reference for
   * which this block was a failure
   * Initialised to itself
   * Updated while computing by FindBeta
   * */
    UINT nearest_tick_ind;

    /** Indicates the previous reference for which this position was queried for
    * beta.
    * Initialised to 0
    * Updated while computing by FindBeta
    */
    UINT recent_querying_ref;
  };

public:
  Twinset(const UINT len, const UINT ref, const UINT level, const UINT luf);

  //~Twinset(void);

  void reset(const UINT ref, UINT &short_num_blocks);
  bool find_position(UINT &query_pos, UINT &num_blocks);
  bool find_position_after_failure(UINT &query_pos, UINT &num_blocks);
  bool find_position_long_query(UINT dist, UINT &query_pos);
  void handle_outside_long_query(const SessionIter sit, const UINT q_pos);
  SessionIter insert_session(std::list<UINT> &positions);
  bool is_empty() const;
  bool is_failed_ref() const;
  UINT creator_ref() const;
  UINT level() const;
  bool print() const;

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

  /** Sessions of this twin-set */
  std::list<Session> _sessions;

  /** Blocks of this twin-set */
  std::vector<Block> _blocks;

  //////////////////////////////////////////////////////////
  // Data-structures needed for maintainance of the twinset
  // afer operations by FindBeta
  //////////////////////////////////////////////////////////

  /** List of indices of consecutive failures. 
   * Indices added when the position to try is queried by FindLCP.
   * recent_querying_ref of the block is set before adding its index in the list
   * List is reset when call to reset is made
   * Indices are dealt in accordance if the last poition added was a success or failure.
   * These indices are used to set corresponding nearest_tick_offset for failures if a success is found.
   * The next and pvs neighbours of the success is also set.
   * Next call to find_position_after_success will mean  that the last index added was a success.
   * Or else, call to done will mean 
  */ std::list<UINT>
      _queried_indices;
  /** Index of the Reference from the twin-set which is trying to find beta */
  UINT _querying_ref_ind;

  /** Iterator to the session where we started finding beta
   * To be used to answer long queries.
   * Set by iniial call; later updated by every success
   */
  SessionIter _beta_query_session;
};

} // end namespace
#endif
