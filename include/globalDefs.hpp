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

/** Header file containing the definitions, macros and declarations used by each
 * of the other modules.
 */
#ifndef LUF_DEFS
#define LUF_DEFS

#include <cstdint>
#include <iostream>
#include <list>
#include <vector>

namespace luf {
using UINT_64 = uint64_t;
using UINT_32 = uint32_t;
using UINT_16 = uint16_t;
using UINT_8 = std::uint8_t;

using UINT = UINT_32;

/** Represent a session of some stack (identifying same blocks in the same
   * session) */
using Session = struct sess {

  /** Number of blocks in all the sessions preceeding and including this */
  // UINT cumulative_num_blocks;

  /** Index in the blocks-vector of the first block in this session */
  UINT first_ind;

  /** Index in the blocks-vector of the first (live i.e. not deleted) block in
   * this session */
  UINT first_live_ind;

  /** Number of blocks in this session (initially) */
  UINT size;

  /** Number of blocks in this session that have not been deleted */
  UINT num_blocks_alive;

  /** Indicates if the current querying reference is the first in this session
   * Initialised to false.
   */
  UINT is_first_ref; 
  
};

using SessionIter = std::list<Session>::iterator;

/** Enum for various possible states (success or errors) rsturned from a
 * function */
enum class ReturnStatus {
  SUCCESS,
  ERR_ARGS,
  ERR_FILE_OPEN,
  ERR_INVALID_INPUT,
  ERR_INVALID_INDEX,
  ERR_LIMIT_EXCEEDS,
  ERR_EXTERNAL,
  HELP
};
//////////////////////////////////////////////////////////
// Parameters collected for stats
//////////////////////////////////////////////////////////
using Stats = struct stat {
  UINT np;             //< total num of pushes
  UINT nlcp;           // total num of lcp
  UINT f_lcp_base;     //< num of failed lcp for base references
  UINT f_lcp_non_base; //< num of failed lcp for non-base references
  UINT s_lcp_base;     //< num of successful lcp for base references
  UINT s_lcp_non_base; //< num of successful lcp for non-base references
  UINT sq_uses;        //< num of times when sq was used to find beta
  UINT depth; //< maximum depth of the stack-tree amongst all base-refrences
};

} // end namespace
#endif
