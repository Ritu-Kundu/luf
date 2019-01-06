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

/** Module implementing the class DegeneratePatternMatch
 */
#ifndef BORDERLESS_CPP
#define BORDERLESS_CPP

#include "Borderless.hpp"

namespace luf {
Borderless::Borderless(const std::string &t, const std::string &a, Stats &s)
    : _cText(t), _cTLen(t.size()), _cAlphabet(a), _cAlphabet_size(a.size()),
      _stats(s), _isReference(_cTLen, false), _isBase(_cTLen, true),
      _max_len(_cTLen, 0), _sq(_cTLen, 0), _inv_twinsets(_cTLen, nullptr),
      _inv_letter_occ(_cTLen, 0), _twinsets(_cTLen, nullptr),
      _recent_failing_base(_cTLen, 0), _nearest_luf_positions(_cTLen, 0),
      _letter_occs(_cAlphabet_size, std::vector<UINT>()), _inv_session(_cTLen),
      _twinsets_collection() {
  _lsf_l.reserve(_cTLen);
  _lsf_r.reserve(_cTLen);
  _hook.reserve(_cTLen);
  _smallest_tested_pos.reserve(_cTLen);

  /* Prepare for LCP queries */
  prepareForRMQ();
  preProcess();
}

Borderless::~Borderless(void) {

  for (auto p : _twinsets_collection) {
    // delete p;
  }
  delete[] _rmq;
}

void Borderless::getLUF(std::vector<UINT> &luf) {
#ifdef DETAILS
  std::cout << "LUF Construction started." << _cText << " " << _cTLen
            << std::endl;
#endif

  for (int i = _cTLen - 1; i >= 0; --i) {
    /* Set LUF value */
    if (_lsf_l[i] == 0) { // starting reference
      luf[i] = _cTLen - i;
    } else {
      UINT j = _lsf_r[i];

      auto curr_ind = _inv_letter_occ[i];
      if (_lsf_l[i] < luf[j]) { // Case a
        luf[i] = j + luf[j] - i;
      } else if (i >= _hook[j]) { // Case b
        luf[i] = luf[j];
      } else { // Case c
        luf[i] = _hook[j] - i;
      }
    }

    /* Update hook, if problematic reference */
    if (_isReference[i]) {
#ifdef DETAILS
      std::cout << "************************************Ref (i, luf[i]): " << i
                << " " << luf[i] << std::endl;
#endif
      if (_sq[i] == 0) { // it is not a square => Case Sq
                         // Do nothing as hook[j] = j already (initialisation)
#ifdef DETAILS
        std::cout << "Not Square: " << std::endl;
#endif
      }                                // case Sq
      else if (luf[i] > _max_len[i]) { // Case 2 but no pos refers to j with
                                       // len > LUF[j]
#ifdef DETAILS
        std::cout << "Special case: Not problematic. " << std::endl;
#endif
      }      // Special Case ends
      else { // Case when hook nees to be updated
#ifdef DETAILS
        std::cout << "Needs updating hook: " << std::endl;
#endif
        find_hook(i, luf[i]);
      }
    }
  }
  std::cout << "HOOK VALUES========\n";
  for (int i = 1; i < _cTLen; ++i) {
    std::cout << i <<  " ";
  }
  std::cout << std::endl;
  for (int i = 1; i < _cTLen; ++i) {
    std::cout << i << " "<< _hook[i] <<  "\n";
  }
  std::cout << std::endl;
}

void Borderless::find_hook(const UINT ref, const UINT luf_ref) {
  LIStack st; // Stack of (len, pos)
#ifdef DETAILS
  std::cout << "Set Hook(ref, luf_ref)::: " << ref << " " << luf_ref
            << std::endl;
#endif
  /* Setting up */
  std::list<UINT> len_list; // list to clean-up the _twinsets for the next call
  UINT q = _hook[ref];
  UINT beta = 0;
  UINT short_sensible = 0;

  // needed only by FindBeta (non-base reference)
  auto pTwin = _inv_twinsets[ref];
  UINT level = 0;
  UINT creator_ref;

  // Positions inserted in the stack between two luf-blocks are collected in
  // this list so that their _nearest_luf_positions can be set
  std::list<UINT> stack_pos;

  if (_isBase[ref]) { // establish blocks

    /* Count number of alphas to be used */
    auto ref_ind = _inv_letter_occ[ref];
    // If never failed, alphas in luf; otherwise alphas between this and its
    // recent-failing base
    // if (_recent_failing_base[ref] == 0) { // never failed

    // index of the letter just after luf
    auto next_to_luf = ref + luf_ref;
    if (next_to_luf < _cTLen) { // there is a character
      short_sensible = _inv_letter_occ[next_to_luf] - ref_ind;
    } else { // luf goes to the end
      auto letter = _cAlphabet.find(_cText[ref]);
      auto &letter_positions = _letter_occs[letter];
      short_sensible = letter_positions.size() - ref_ind;
    }
#ifdef DETAILS
    std::cout << "Never Failed base reference: num-alphas: " << short_sensible
              << std::endl;
#endif

    beta = find_beta_base(q, ref, luf_ref, luf_ref, short_sensible);
  } else {
#ifdef DETAILS
    std::cout << "Twinset at ref: " << ref << std::endl;
    pTwin->print();
#endif
    // Reset Twinset for this reference
    pTwin->reset(ref, _inv_session[ref], short_sensible);
    level = pTwin->level();
    creator_ref = pTwin->creator_ref();
#ifdef DETAILS
    std::cout << "num-blocs(to test) : " << short_sensible << std::endl;
#endif
    beta = find_beta(q, ref, luf_ref, short_sensible);
  }

#ifdef DETAILS
  std::cout << "Beta (q): " << beta << " " << q << std::endl;
#endif
  while (beta != 0) {
    auto si = q - beta; // chopping index pushed in stack
    handle_popping(ref, luf_ref, beta, q, level + 1, st, len_list, stack_pos);
    st.push({beta, si});
    ++(_stats.np);
    _isBase[si] = false;
    _nearest_luf_positions[si] = 0; // reset it again
#ifdef DETAILS
    std::cout << "::::::::::::PUSHED:::::Pushed in stack: " << beta << " " << si
              << std::endl;
#endif
    q = _hook[si];
    if (_isBase[ref]) { // establish blocks
      beta = find_beta_base(q, ref, luf_ref, beta, short_sensible);
    } else {
      beta = find_beta(q, ref, luf_ref, short_sensible);
    }
#ifdef DETAILS
    std::cout << "Beta (q): " << beta << " " << q << std::endl;
#endif
  }

  // set hook for the remaining elements of the stack
  if (!st.empty()) { // useful when beta is 0
    q = _hook[st.top().second];
#ifdef DETAILS
    std::cout << "Final setting hook: " << q << std::endl;
#endif
    // stack_pos in this case will not be of any use
    handle_popping(ref, luf_ref, luf_ref + 1, q, level + 1, st, len_list,
                   stack_pos); // len in stack <= luf_ref
  }
  _hook[ref] = q;

  // Clean-up the data-structures for the next call
  for (auto l : len_list) {
#ifdef DETAILS
    _twinsets[l]->print();
#endif
    _twinsets[l] = nullptr;
  }
}

/** Pops smaller elements than beta from the stack
 * Updates their hook
 * Pushes them in corresponding twinset
 **/
void Borderless::handle_popping(const UINT ref, const UINT luf, const UINT beta,
                                const UINT q, const UINT level, LIStack &st,
                                std::list<UINT> &len_list,
                                std::list<UINT> &stack_pos) {

  while (!st.empty() && st.top().first < beta) {
    auto ind = st.top().second;
    auto len = st.top().first;
    st.pop();
    // Needed to update nearest luf pos (required for lon-query)
    stack_pos.push_back(ind);

    // Create Twinset, if needed
    if (_twinsets[len] == nullptr) { // first element of this length
      //  create twin-setting
      // std::cout << "CREATING TWIN\n";
      auto pTwin = new Twinset{len, ref, level, luf};
      _twinsets[len] = pTwin;
      _twinsets_collection.push_back(pTwin);
      len_list.push_back(len);
      if (_stats.depth < level) {
        _stats.depth = level;
      }
    }
    auto pTwinset = _twinsets[len];
    // pTwinset->print();
    // Update hook and twin-set for every index of this len in this session
    _hook[ind] = q;
    _inv_twinsets[ind] = pTwinset;
    // collect all positions of this length => session of twinset
    std::list<UINT> same_len_pos;
    same_len_pos.push_back(ind);
    while (!st.empty() && st.top().first == len) {
      auto ind = st.top().second;
      st.pop();
      _hook[ind] = q;
      _inv_twinsets[ind] = pTwinset;
      same_len_pos.push_front(ind); // for deceasing order
    }
    // Add these positions in the twinset
    auto sit = pTwinset->insert_session(same_len_pos);

    // Update their session iterators
    for (auto p : same_len_pos) {
      _inv_session[p] = sit;
    }
  }

  // If this is luf-block, update nearest luf-block-position for all the indices
  // collected from the previous luf-block-position
  if (beta == luf) {
    auto luf_block_position = q - beta;
    for (auto p : stack_pos) {
      // blocks
      _nearest_luf_positions[p] = luf_block_position;
    }
    stack_pos.clear();
  }
}

/** Modifies the sessions of the twinset of ref;
 * Modifies the session iterator (The modified will be used in the next call)
 * Assumes session iter points to the session past hook (pos);
 * returns beta
 **/
UINT Borderless::find_beta(
    const UINT pos, const UINT ref, const UINT luf_ref,
    const UINT
        short_sensible) { // Finds the length of part to chop ending at pos
  bool done = false;
  UINT beta = 0;
  UINT stop = 0;
  if (pos > luf_ref) {
    stop = pos - luf_ref;
  }

#ifdef DETAILS
  std::cout << "Finding beta: pos, Stop, short_sensible: " << pos << " " << stop
            << " " << short_sensible << std::endl;
#endif
  auto pTwin = _inv_twinsets[ref];
  UINT query_pos = 0;
  UINT num_blocks = 0;
  bool is_first_in_session = pTwin->is_first_in_session();
  bool is_found = false;
  // Handle the case when not the first in the session
  // Only one query and return
  if (!is_first_in_session) {
#ifdef DETAILS
    pTwin->print_session_succ_list();
    ;
#endif
    if (!pTwin->find_position_quick(query_pos)) { // position not found
#ifdef DETAILS
      std::cout << "Nothing to check \n ";
#endif
      return beta;
    }
    auto l = getLongestCommonPrefix(query_pos, ref);
    ++(_stats.nlcp);
    if (query_pos + l >= pos) { // l reaches pos or go beyond: successful lcp
      beta = pos - query_pos;
#ifdef DETAILS
      std::cout << "%%% SUCCESS QUICK LCP: beta pos query_pos l " << beta << " "
                << pos << " " << query_pos << " " << l << std::endl;
#endif
      ++(_stats.s_lcp_non_base);
    } else { // failed query
      ++(_stats.f_lcp_non_base);
    }
    return beta;
  } // Case when not the first in the session done

  if (!pTwin->find_position_first_in_session(
          query_pos, num_blocks)) { // position not found
#ifdef DETAILS
    std::cout << "Nothing to check \n ";
#endif
    return beta;
  }

  /* Try to find short query */
  while (query_pos >= stop && num_blocks <= short_sensible) {
    auto l = getLongestCommonPrefix(query_pos, ref);
    ++(_stats.nlcp);
    if (query_pos + l >= pos) { // l reaches pos or go beyond: successful lcp
      beta = pos - query_pos;
#ifdef DETAILS
      std::cout << "%%% SUCCESS SHORT LCP: beta pos query_pos l " << beta << " "
                << pos << " " << query_pos << " " << l << std::endl;
      ;
#endif
      done = true;
      ++(_stats.s_lcp_non_base);
      break;
    } else { // failed query
#ifdef DETAILS
      std::cout << "%%% FAILED SHORT LCP: beta pos query_pos l " << beta << " "
                << pos << " " << query_pos << " " << l << std::endl;
      ;
#endif
      ++(_stats.f_lcp_non_base);
      if (!pTwin->find_position_after_failure(
              query_pos,
              num_blocks)) { // position not found
#ifdef DETAILS
        std::cout << "Nothing to check \n ";
#endif
        return beta;
      }
    }
  } // short-query attempts done

  /** Try long query */
  if (!done) {
    // Find query position
    if (pTwin->is_failed_ref()) { // it has failed for some previous reference
      if (!pTwin->find_position_long_query(short_sensible, query_pos)) {
#ifdef DETAILS
        std::cout << "Nothing to check \n ";
#endif
        return beta;
      }
    } else { // if this reference hasnt been seen
      auto creator_ref = pTwin->creator_ref();
      auto luf_block_pos = _nearest_luf_positions[pos];
      if (luf_block_pos == 0 ||
          luf_block_pos <
              (creator_ref - ref)) { // no luf-block or valid query-pos
#ifdef DETAILS
        std::cout << "Nothing to check \n ";
#endif
        return beta;
      }
      query_pos = luf_block_pos - (creator_ref - ref);
    }

    auto l = getLongestCommonPrefix(query_pos, ref);
    ++(_stats.nlcp);
    if (query_pos + l >= pos) { // l reaches pos or go beyond: successful lcp
      beta = pos - query_pos;
#ifdef DETAILS
      std::cout << "%%% SUCCESS LONG LCP: beta pos query_pos l " << beta << " "
                << pos << " " << query_pos << " " << l << std::endl;
      ;
#endif
      ++(_stats.s_lcp_non_base);
      if (!pTwin->is_failed_ref()) { // quer-pos calculated outside twinsets
        pTwin->handle_outside_long_query(_inv_session[query_pos], query_pos);
      }
    } else {
      ++(_stats.f_lcp_non_base);
    }
  } // long-query attempts done
  return beta;
}

UINT Borderless::find_beta_base(
    const UINT pos, const UINT ref, const UINT luf_ref, const UINT useSq_len,
    const UINT
        test_num_alpha) { // Finds the length of part to chop ending at pos
  UINT beta = 0;
  UINT stop = 0;
  if (pos > luf_ref) {
    stop = pos - luf_ref;
  }
#ifdef DETAILS
  std::cout << "Finding INITIAL  beta: pos, stop, test_num_alpha, sq[pos], "
               "useSq_len : "
            << pos << " " << stop << " " << test_num_alpha << " " << _sq[pos]
            << " " << useSq_len << std::endl;
#endif
  /* If never failed base reference */
  if (_recent_failing_base[ref] == 0 && _sq[pos] != 0 &&
      (useSq_len >= _sq[pos])) { // never failed. and a usable square
    beta = _sq[pos];
    ++(_stats.sq_uses);
#ifdef DETAILS
    std::cout << "Used Square: " << std::endl;
#endif
  }

  else if (_recent_failing_base[ref] != 0 &&
           pos >= _smallest_tested_pos[_recent_failing_base[ref]]) { // failed
    /* Otherwise failed base-reference (make use of sq until left-most tested
     * position) */
    beta = _sq[pos];
    ++(_stats.sq_uses);
#ifdef DETAILS
    std::cout << "Failed base ref: failed-for : " << _recent_failing_base[ref]
              << " : " << std::endl;
#endif
  }

  else if (_inv_letter_occ[pos] > 0) { // use LCP
    auto letter = _cAlphabet.find(_cText[pos]);
    auto &letter_positions = _letter_occs[letter];

    auto start_ind = _inv_letter_occ[pos] -
                     1; // letter occs are in ascending order in the array
    // possible starting position for beta
    auto current_ind = start_ind;
    auto current_pos = letter_positions[start_ind];
    UINT num_letters_checked = 0;
    while (current_pos >= stop && num_letters_checked < test_num_alpha) {
      auto l = getLongestCommonPrefix(current_pos, ref);
      ++(_stats.nlcp);

#ifdef DETAILS
      std::cout << "LCP: current, ref: " << current_pos << " " << ref
                << std::endl;
#endif
      if (current_pos + l >= pos) { // l reaches pos or go beyond
        beta = pos - current_pos;
#ifdef DETAILS
        std::cout << "%%% beta pos current l " << beta << " " << pos << " "
                  << current_pos << " " << l << std::endl;
        ;
#endif
        ++(_stats.s_lcp_base);
        break;
      } else { // fail-lcp
        _recent_failing_base[current_ind] = ref;
        _smallest_tested_pos[ref] = current_ind;
        ++(_stats.f_lcp_base);
      }
      if (current_ind == 0) { // last letter checked
        break;
      } else {
        --current_ind;
        current_pos = letter_positions[current_ind];
        ++num_letters_checked;
      }
    }
  }
  return beta;
}

ReturnStatus Borderless::preProcess() {
  const UINT len = _cTLen;

  for (UINT i = 0; i < len; ++i) {
    UINT max_l = 0; // length of the longest last successor factor index
    UINT lsf = 0;   // index of the longest last successor factor
    UINT lsquare = 0;
    for (UINT j = i + 1; j < len; ++j) {
      UINT l = getLongestCommonPrefix(i, j);
      if (l >= max_l) { // to find LSF
        max_l = l;
        lsf = j;
      }
      if (i + l >= j) { // Shortest square
        _sq[j] = j - i; // As we are moving left to right, it will eventually
                        // have correct values
      }
    } // checking at each j ends

    _lsf_l.push_back(max_l);
    _lsf_r.push_back(lsf);

    if (max_l != 0) {           // it refers to some position
      _isReference[lsf] = true; // mark the position i refers to as a reference
      if (max_l > _max_len[lsf]) { // this is longer length being referecned
        _max_len[lsf] = max_l;
      }
    } else { // starting reference
    }
    // Add pos to corresponding alphabet's occs
    auto letter = _cAlphabet.find(_cText[i]);
    if (letter == std::string::npos) { // invalid character
      std::cout << "Invalid character " << _cText[i] << "found at " << i
                << ". Valid alphabet: " << _cAlphabet << std::endl;
      return ReturnStatus::ERR_INVALID_INPUT;
    }
    _letter_occs[letter].push_back(i);
    _inv_letter_occ[i] = _letter_occs[letter].size() - 1;
    // initialise hook for this position
    _hook[i] = i;
    _smallest_tested_pos[i] = i;
  } // each i ends

  // Shrink the vectors
  for (UINT i = 0; i < _cAlphabet_size; ++i) {
    if (!_letter_occs[i].empty()) {
      _letter_occs[i].shrink_to_fit();
    }
  }

#ifdef DETAILS
  std::cout << "s\tchar\t\tlsf_l\t\tlsf_r\t\tsq\t\tisReference\tletter-"
               "occ\tHook\t\tMax_len "
            << std::endl;
  for (UINT s = 0; s < len; ++s) {
    std::cout << s << "\t" << _cText[s] << "\t\t" << _lsf_l[s] << "\t\t"
              << _lsf_r[s] << "\t\t" << _sq[s] << "\t\t" << _isReference[s]
              << "\t\t" << _inv_letter_occ[s] << "\t\t" << _hook[s] << "\t\t"
              << _max_len[s] << std::endl;
  }
#endif
  return ReturnStatus::SUCCESS;
}

UINT Borderless::getLongestCommonPrefix(const UINT l, const UINT r) {
  UINT l_rmq, r_rmq;
  UINT rx = _csa.isa[l];
  UINT ry = _csa.isa[r];
  if (rx < ry) {
    l_rmq = rx;
    r_rmq = ry;
  } else {
    l_rmq = ry;
    r_rmq = rx;
  }
  auto lcp_len = _lcp[getRMQ(l_rmq + 1, r_rmq)];
  // std::cout << "\nRMQ: " <<l <<" "<<r<<": "<<l_rmq+1<<" "<<r_rmq<<" : "<<" :
  // "<<getRMQ(l_rmq+1, r_rmq)<<" : "<<lcp_len << std::endl;
  return lcp_len;
}

ReturnStatus Borderless::prepareForRMQ() {
  const UINT len = _cTLen + 1;
  /* Compute Suffix Array */
  sdsl::construct_im(_csa, _cText, 1); // 1 for alphabet type
  sdsl::construct_im(_lcp, _cText, 1); // 1 for alphabet type

  /* Prepare LCP array for RMQ */
  UINT lgn = flog2(len + 1);
  _rmq = new UINT[(len)*lgn];
  rmq_preprocess();
  return ReturnStatus::SUCCESS;
}

void Borderless::rmq_preprocess() {
  const UINT len = _cTLen + 1;
  UINT lgn = flog2(len);
  // initialize rmq for the intervals with length $1$
  for (UINT i = 0; i < len; i++) {
    _rmq[i * lgn] = i;
  }
  // compute values from smaller to bigger intervals
  for (UINT_64 j = 1; 1 << j <= len; j++) {
    for (UINT_64 i = 0; i + (1 << j) - 1 < len; i++) {
      if (_lcp[_rmq[i * lgn + j - 1]] <
          _lcp[_rmq[(i + (1 << (j - 1))) * lgn + j - 1]]) {
        _rmq[i * lgn + j] = _rmq[i * lgn + j - 1];
      } else {
        _rmq[i * lgn + j] = _rmq[(i + (1 << (j - 1))) * lgn + j - 1];
      }
    }
  }
}

UINT Borderless::getRMQ(UINT i, UINT j) {
  const UINT_64 len = _cTLen + 1;
  UINT lgn = flog2(len);
  if (i > j) {
    UINT tmp = j;
    j = i;
    i = tmp;
  }
  if (i == j)
    return i;
  UINT k = flog2(j - i + 1);
  UINT a = _rmq[i * lgn + k];
  UINT b = _rmq[(j - (1 << k) + 1) * lgn + k];
  return _lcp[a] > _lcp[b] ? b : a;
}

UINT Borderless::flog2(UINT v) {
  union {
    UINT u[2];
    double d;
  } t;
  t.u[1] = 0x43300000;
  t.u[0] = v;
  t.d -= 4503599627370496.0;
  return (t.u[1] >> 20) - 0x3FF;
}

void Borderless::naiveLUF(std::vector<UINT> &luf) {
  for (UINT i = 0; i + 1 < _cTLen; ++i) {
    luf[i] = getRightZeroBorder(i);
    // std::cout << i << " "<<luf [i] << std::endl;;
  }
  // std::cout <<std::endl;
  luf[_cTLen - 1] = 1;
}
// #assumes start > _cTLen
UINT Borderless::getRightZeroBorder(UINT start) {
  UINT m = _cTLen - start;
  // Build border
  std::vector<int> border(m + 1);
  int i = 0;
  border[0] = -1;
  int j = border[0];
  while (i < m) {
    while (j > -1 && _cText[i + start] != _cText[j + start]) {
      j = border[j];
    }
    border[++i] = ++j;
    // std::cout << border[i] << " ";
  }
  // std::cout << std::endl;
  int len = 0;
  for (auto pos = m; pos > 0; --pos) {
    if (border[pos] == 0) {
      len = pos;
      break;
    }
  }
  return len;
}

} // end name space

#endif
