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

/** Module defining the class that implements the main algorithm
 */
#ifndef BORDERLESS_HPP
#define BORDERLESS_HPP

#include "Twinset.hpp"
#include "globalDefs.hpp"

#define DETAILS

// NOTE: SDSL SA Construction doesn,t work with 0 as a character
#include <sdsl/lcp.hpp>
#include <sdsl/suffix_arrays.hpp>

#include <cstddef>
#include <fstream>
#include <iterator>
//#include <queue>
//#include <sdsl/rmq_support.hpp>
#include <stack>
#include <utility>

namespace luf {
/** Class DegeneratePatternMatch
 * This class provides for the following:
 *  - Returns the occurences of pattern(degenerate) in input-sequence.
 * It prepares for constant-time LCA/LCP queries while being constructed.
 */

using LIStack =
    std::stack<std::pair<UINT, UINT>>; // stack containing pairs (length, index)

class Borderless {
public:
  Borderless(const std::string &t, const std::string &a, Stats &s);

  /** Destructor
   */
  ~Borderless();

  /** Returns the indices of occurences in the given list of the pth pattern */
  void getLUF(std::vector<UINT> &luf);
  void naiveLUF(std::vector<UINT> &luf);

  //////////////////////////////////////////////////////////
  // PRIVATE
  //////////////////////////////////////////////////////////
private:
  const std::string &_cAlphabet; //< String representing alphabet

  const UINT _cAlphabet_size; //< Size of the alphabet

  const std::string &_cText; //< reference to the text (inout seuence)

  const UINT _cTLen; //< Length of the input sequence

  //////////////////////////////////////////////////////////
  // Data-structures for to help pre-processing & naive
  //////////////////////////////////////////////////////////
  sdsl::csa_bitcompressed<> _csa; //< Suffix Array and Rank array

  sdsl::lcp_bitcompressed<> _lcp; //< LCP array

  UINT *_rmq; //< For RMQs

  //////////////////////////////////////////////////////////
  // Data-structures filled by pre-processing stage
  //////////////////////////////////////////////////////////
  /** LSF (Longest successor Factor-length Array */
  std::vector<UINT> _lsf_l;
  /** LSF (Longest successor Factor-reference Array */
  std::vector<UINT> _lsf_r;
  /** Shortest Square-length Array */
  std::vector<UINT> _sq;

  // std::vector<std::list<INT>> rj;
  /** Array containng maximum referenced length for a reference
   * _max_len [j] = l => that for the reference j, l is the maximum length of
   * lsf of all i's that refere to j.
   * Used to decide if j is a complex-reference i.e. it;s hook info will be used
   * to find LUF for some i
   * */
  std::vector<UINT> _max_len;
  /** Hook Array
   * Initiakised to i while pre-processing
   * Keeps updating while computing by FindHook
  */
  std::vector<UINT> _hook;

  /** Array indicationg if a position is a reference */
  std::vector<bool> _isReference;

  /** Array indicationg if a position is a base-reference
   * Initialised to false
   * Updated while computing
  */
  std::vector<bool> _isBase;

  // std::vector<bool> isPvsUpdated;

  //////////////////////////////////////////////////////////
  // Data-structures for FindBeta (base-references)
  //////////////////////////////////////////////////////////
  /** 2-d Array indicationg indices/positions of the occurrences of each letter
   * in the text
   * _letter_occs [c] [i] = p => ith occurrence of character (mapped at position
   * c in alphabet) is at position p in the text.
   * #rows = alphabet-size
   * size of each row c = number of occurrences of character c (i.e. character
   * at position c in alphabet)
   * Arranged in ascending order of indices
   */
  std::vector<std::vector<UINT>> _letter_occs;

  /** Array indicationg index of the the letter at position i (say c in
   * alphabet) in its occurrences array _letter_occs[c] */
  std::vector<UINT> _inv_letter_occ;

  /** Array indicationg the previous base-reference for which this position was
   * a failure
   * Initialised to 0
   * Updated while computing by FindBeta
   */
  std::vector<UINT> _recent_failing_base;

  //////////////////////////////////////////////////////////
  // Data-structures for FindBeta (non-base-references)
  //////////////////////////////////////////////////////////
  /** Array contaning pointer to the current twinset for each position of the
   * text
   *
  */
  std::vector<Twinset*> _inv_twinsets;

  /** Array contaning iterator to the current session in the current twinset for
   * each position of the
   * text
   *
  */
  std::vector<SessionIter> _inv_session;

  /** Twinset corresponding to  each length in a particular stack
   * Used while creating twinsets for that stack
   * The maximum length of a twin-set can be luf of the corresponding reference.
   * It is kept global so that every reference doesn't do work equivalent to its
   * luf in initialising these pointers.
   * A reference uses it for creating twin-sets in FindHook call
   * It needs cleaning up at the end of FindHook call for the next reference.
   * */
  std::vector<Twinset*> _twinsets;

  /** Pointers to all the twinsets created
   * Used in deleting them later
   * */
  std::list<Twinset *> _twinsets_collection;

  /** Array of positions indicating the index of nearest luf-block to this
   * position
   * Used only if this position doesn't have luf-block cut
   * Used to answer long-query positions for references in twin-set not seen by
   * any other
   * Initially 0
   * Updated while popping stack
   * Positions inserted in the stack between two luf-blocks are collected in a
   * list
   * */
  std::vector<UINT> _nearest_luf_positions;

  /** Array contaning position for which a base refernce has made the testing
   * Initially set to itself
  */
  std::vector<UINT> _smallest_tested_pos;

  //////////////////////////////////////////////////////////
  // Parameters collected for stats
  //////////////////////////////////////////////////////////
  Stats &_stats;

  void find_hook(const UINT ref, const UINT luf_ref);
  UINT find_beta_base(const UINT pos, const UINT ref, const UINT luf_ref,
                      const UINT useSq_len, const UINT test_num_alpha);
  UINT find_beta(const UINT pos, const UINT ref, const UINT luf_ref,
                 const UINT short_sensible);
  void handle_popping(const UINT ref, const UINT luf, const UINT beta,
                      const UINT q, const UINT level, LIStack &st,
                      std::list<UINT> &len_list, std::list<UINT> &stack_pos);

  //////////////////////////////////////////////////////////
  // Function to preprocess and its helper functions
  //////////////////////////////////////////////////////////
  ReturnStatus preProcess();

  /** Preapres for constant time RMQ queries */
  ReturnStatus prepareForRMQ();

  /** Returns Longest Common Prefix of suffices starting at given l and r */
  UINT getLongestCommonPrefix(const UINT l, const UINT r);

  //////////////////////////////////////////////////////////
  // Helper functions for answering RMQ
  //////////////////////////////////////////////////////////
  void rmq_preprocess();
  UINT getRMQ(UINT i, UINT j);
  UINT flog2(UINT v);

  UINT getRightZeroBorder(UINT start);
};

} // end name space

#endif
