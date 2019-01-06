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

/** Module containing main() method which reads and preprocesses the input and
 * pattern files in the form required for identification of the ocuurences.
 */
#include "../include/Borderless.hpp"
#include "../include/globalDefs.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <random>

#define DETAILS
//#define RANDOM
//#define EXHAUSTIVE

using namespace luf;

int test(std::string &text);
const std::string cAlphabet = "abcd";
std::ofstream statfile("stats.txt");

int main(int argc, char **argv) {

  if (!statfile.is_open()) {
    std::cerr << "Cannot open  stats file \n";
    return 1;
  }

  statfile << "Len\t"
           << "depth \t"
           << "num_push\t"
           << "num_lcp \t"
           << "num_sq \t"
           << "num_Succ_lcp_base \t"
           << "num_Succ_lcp_non_base \t"
           << "num_fail_lcp_base \t"
           << "num_fail_lcp_non_base \t"
           << "time\n";

#ifdef DETAILS
// Individual strings that have caused problems in the past

  std::vector<std::string> text_array;
  /**
  // Example of base ref-dsitance concept => use previous failed for
  references
  // that are on right of the hook
  text_array.push_back("bcbbabbccabbcbbabcbbab");

  text_array.push_back("aabaabbaabaabbbaabaabbaabaabbbbaabaabbaabaabbbaabaabbaa"
                       "baabbbbaabaabbaabaabbbaabaabbaabaabbbbaabaabbaabaabbbaa"
                       "baabbaabaabbbbaabaabbaabaabbbaabaabbaabaabbbbaabaabbaab"
                       "aabbbaabaabbaabaabbbbaabaabbaabaabbbaabaabbaabaabbbbaab"
                       "aabbaabaabbbaabaabbaabaabbbbaabaabbaabaabbbaabaabbaabaa"
                       "bbbbaabaabbaabaabbbaabaabbaabaabbbb");

  text_array.push_back("aabaabbaabaabbbaaaaabaabbaabaabbbc");
  text_array.push_back("bbbbabbababaabaabbbabbbaababbabaaaabaababbbbaabbaabbbba"
                       "abaaabaabaabbbaaabbbabaabbaabbabbaaabbabababbabaabbbaaa"
                       "baaaaabaabbbbbbabbabbbaabbbbbabaabaaabaabaabbaaabababba"
                       "ababbaaabbbbbabbbbbbabaaaaaabaabbbbbaaaabaabbbbbaab");

  text_array.push_back(
      "ababbbbaaaabbbaabaaaababaabbaabbabaabbaaabbaabbabaabbaabbabaabbaaa");
  // // example do not delete refrence  in fail followed by success case
  text_array.push_back("bbabaaaabaaaabaaaaabbaaaaabaabababbbaababbbbaaabbbabbab"
                       "abbbbbbbabbabbbbbbbabba");
  text_array.push_back("bbababbaabbababbbaabbbbbbaabbabbababbaabbabbaba"); //
  // example for need of safest index
  text_array.push_back("bbaabbbbaabaabbbbababbbbaaaabaaaaaabbbbbababaababaaabba"
                       "babbaabbababbbbbbbbabbabbaabbabaaaaababbbabaabbabbbbbaa"
                       "aababbbbaabaabbbababaabaabbbbbaabbbabbabbaabaabaaabaaaa"
                       "baabbaabbbbbbabbaababbbbaabaababaaaaabaabaa");
  // // example all delete failed
  text_array.push_back("baababbbabaaaabbbbaababbaaaaaabbbaabbaabaababbbaaababab"
                       "bbbbbbabaabbabbaabaaabbabbbaaaabaaaaaabbaabababbaaaaabb"
                       "aaabbbabaabbabbaa");
  // // twin-set of a refernce becomes empty, so deleted before the ref is
  // considered
  text_array.push_back("aaaabaabbbbbbababbbaaaabbabaaaaabbabbbababbbbbbabb");
  //
  // Whole session containing ref deleted
  text_array.push_back("abbbabbbbbbbbbaababbabbbbabbbbbababaabababaabbaababba"
                       "ababaaababaabbaaababbbabbaaaaabaaabaabbbbaaabaababbab"
                       "baaaaaabbaaaabbbaaabaaaaaabababaababbbaabababbaaabbaa"
                       "aaaaaaabaabaababaabbababbabababbaabbbabbaaababbaaabab"
                       "baaa");

  text_array.push_back("abbbabbbbbbbbbaababbabbbbabbbbbababaabababaabbaababbaab"
                       "abaaababaabbaaababbbabbaaaaabaaabaabbbbaaabaababbabbaaa"
                       "aaabbaaaabbbaaabaaaaaabababaababbbaabababbaaabbaaaaaaaa"
                       "abaabaababaabbababbabababbaabbbabbaaababbaaababbaaa");
  text_array.push_back("abbbbababbabbbabbbbbabbbbababbabbbabbbbbababbabbbabbbba"
                       "babbabbbabbbbbababbabbbabbbbababbabbbabbbbbabbbbababbab"
                       "bbabbbbb");

  text_array.push_back("abbbbababbabbbabbbbbabbbbababbabbbabbbbbababbabbbabbb"
                       "bababbabbbabbbbbababbabbbabbbbababbabbbabbbbbabbbbaba"
                       "bbabbbabbbbbabbbbababbabbbabbbbbabbbbababbabbbabbbbba"
                       "babbabbbabbbbababbabbbabbbbbababbabbbabbbbababbabbbab"
                       "bbbbabbbbababbabbbabbbbb");
  text_array.push_back(
      "abbbbababbabbbabbbbbabbbbababbabbbabbbbbababbabbbabbbbababbabbbabbbbbaba"
      "bbabbbabbbbababbabbbabbbbbabbbbababbabbbabbbbbabbbbababbabbbabbbbbabbbba"
      "babbabbbabbbbbababbabbbabbbbababbabbbabbbbbababbabbbabbbbababbabbbabbbbb"
      "abbbbababbabbbabbbbbabbbbababbabbbabbbbbabbbbababbabbbabbbbbababbabbbabb"
      "bbababbabbbabbbbbababbabbbabbbbababbabbbabbbbbabbbbababbabbbabbbbb");

  text_array.push_back("bbbbababaabbbbaabbbbbbbbabbabbbbbaabbbbabba");
  text_array.push_back("bbbaaabbabbabaaabbbbbabbbbabbababbbbbbaaabbabaababaabaa"
                       "aaabbaaaaabaabbbbbaabaabbbbbaabaaa");
  */
  //text_array.push_back("daabbabaabbaababbabab");
  text_array.push_back("daabaabbaabaabbbaabaabbaabaabbb");

  for (auto i = 0; i < text_array.size(); ++i) {
    auto num_errors = test(text_array[i]);
    if (num_errors > 0) {
      std::cout << "FOUND ERROR :::::::::::::\n";
      break;
    }
  }

/**
  const int N = 20;
  for (auto n = 1; n <= N; ++n) {
    std::string text(2 * n + 7, 'a');

    text[n + 2] = 'b';
    text[n + 3] = 'b';
    text[2 * n + 6] = 'b';
    auto num_errors = test(text);
    if (num_errors > 0) {
      std::cout << "FOUND ERROR :::::::::::::\n";
      break;
    }
  }
  **/
/**
  const int N = 20;
  for (auto n = 1; n <= N; ++n) {
    std::string text(7 * n + 10, 'a');
    text[n] = 'b';
    text[2 * n + 2] = 'b';
    text[3 * n + 3] = 'b';
    text[4 * n + 6] = 'b';
    text[5 * n + 7] = 'b';
    text[6 * n + 9] = 'b';
    auto num_errors = test(text);
    if (num_errors > 0) {
      std::cout << "FOUND ERROR :::::::::::::\n";
      break;
    }
  }
  **/

// Stack-size Worst case examples
/**
std::vector<std::string> ui = {"a", "ab","abb", "abbb", "abbbb", "abbbbb",
"abbbbbb", "abbbbbbb", "abbbbbbbb", "abbbbbbbbb","abbbbbbbbbb"};
auto block_type = 3;

std::string str = "";
std::string block = "a";
for (auto i =0; i < block_type; ++i){
  str = str + block + str;
  block = block + "b";
}
str = str+std::string(block, 0, block.size()-1);

// Another way of generating the Stack-size Worst case examples
/**
std::string str = "";
std::string block = "a";
auto block_type = 3;
for (auto i =0; i < block_type-1; ++i){
  str = str + block + str;
  block = block + "b";
}
str = str+block;
std::string text(str+str);
**/

//////////////////////////////////////////////////////////////////////////
/** Reviwer's Counter Example */
// std::string text("abdabcabdaaaabdabcabdaaaabdabcabdabcabdaaaabd");
// std::string text("caabbaabaabbaabaabbabaabba");
// std::string text("dabaabaabbab");
// std::string text("daabbaabaabbabaabbbabba");
// std::string text("eabacadabacadabacad");
// std::string
// text("dabbbaabaabaabaabaabbaababbbaabaabaabaabaabbaababbbaabaabaabaabaabbaab");
/**
auto num_errors = test(text);
    if (num_errors > 0) {
      std::cout << "FOUND ERROR :::::::::::::\n";
    }
**/

/** Reviwer's Counter Example pattern */
/**
int K = 100;
for (auto k = 5; k <= K; ++k) {
  std::string b_k(k, 'b');
  std::string aa_bk = "";
  for (auto i = k; i > 0; --i) {
    aa_bk += "aa" + b_k.substr(0, i);
  }
  std::string aab_k2 = "";
  int k2 = k * k;
  for (int i = 0; i < k2; ++i) {
    aab_k2 += "aab";
  }
  std::string w = "a" + b_k + "b" + aab_k2 + aa_bk;
  std::string text(w + w + w);
  std::cout << "\n\n$$$$$$$$$$$$$$$$$$$$$$$$$$ " << k << " $$$$$$$\n\n";
  td::cout << "w: " << w << std::endl;
  statfile << k << "\t";
  auto num_errors = test(text);
    if (num_errors > 0) {
      std::cout << "FOUND ERROR :::::::::::::\n";
      break;
    }
}
**/

/** MANAL Example */
/**
int K = 100;
for (auto k = 5; k <= K; ++k) {
  std::string b_k(k, 'b');
  std::string aa_bk = "";
  for (auto i = 5; i > 0; --i) {
    aa_bk += "aa" + b_k.substr(0, k - i);
  }
  aa_bk += "a" + b_k;
  std::string a_k(k * k, 'a');
  std::string w = a_k + aa_bk;
  std::string text(w + w + w);
  std::cout << "\n\n$$$$$$$$$$$$$$$$$$$$$$$$$$ " << k << " $$$$$$$\n\n";
  std::cout << "w: " << w << std::endl;
  statfile << k << "\t";
  auto num_errors = test(text);
    if (num_errors > 0) {
      std::cout << "FOUND ERROR :::::::::::::\n";
      break;
    }
}
**/

/** My example : Only distance for non-base reference will not work */
/**
int K = 100;
for (auto k = 0; k <= K; ++k) {
  std::string ab_k = "";
  for (auto i = 0; i < k; ++i) {
    ab_k += "ab";
  }
  std::string w (ab_k + "a" + "abb" + ab_k + "abbb");

  std::string text(w + w);
  std::cout << "\n\n$$$$$$$$$$$$$$$$$$$$$$$$$$ " << k << " $$$$$$$\n\n";
  std::cout << "w: " << w << std::endl;
  statfile << k << "\t";
  auto num_errors = test(text);
    if (num_errors > 0) {
      std::cout << "FOUND ERROR :::::::::::::\n";
      break;
    }
}
**/


#endif

//########################################################################################################
#ifdef RANDOM

  // Generate Strings
  const int cNum_strings = 1000000;
  const int cMax_len = 200;

  std::random_device rd;
  std::mt19937 gen_len(rd());
  std::uniform_int_distribution<> dis_len(2, cMax_len);
  std::mt19937_64 gen_char{std::random_device()()};
  std::uniform_int_distribution<size_t> dist_char{0, cAlphabet.length() - 1};

  int failed_strings = 0;
  for (int i = 0; i < cNum_strings; ++i) {
    // Choose length
    auto len = dis_len(gen_len);
    // Generate string
    std::string text;
    std::generate_n(std::back_inserter(text), len,
                    [&] { return cAlphabet[dist_char(gen_char)]; });
    std::cout << i << std::endl;

    // Find LUF
    auto num_errors = test(text);
    if (num_errors > 0) {
      std::cout << "FOUND ERROR :::::::::::::\n";
      ++failed_strings;
      break;
    }
  }
  if (failed_strings == 0) {
    std::cout << "All successfull" << std::endl;
  } else {
    std::cout << "num Failed Strings: " << failed_strings << std::endl;
  }

#endif
//########################################################################################################
#ifdef EXHAUSTIVE
  // Generate Strings
  std::cout << "************Starting" << std::endl;
  const int max_len = 20;
  int failed_strings = 0;
  const int num_rows = 1 << (max_len - 1); // 2^(max_len-1)
  const int num_col = max_len;
  std::vector<std::vector<char>> all_strings(num_rows,
                                             std::vector<char>(num_col, 'a'));
  ;
  // Fill the strings array
  int block_size = 1;
  for (auto i = 1; i < num_col; ++i) {
    const int start_block_ind = block_size;
    const int next_block_ind = 2 * block_size;
    // std::cout << "col bs sind nind: " << i << " "<<block_size << " " <<
    // start_block_ind << " " <<next_block_ind<<std::endl;
    for (auto j = start_block_ind; j < num_rows; j = j + next_block_ind) {
      for (auto k = 0; k < block_size; ++k) {
        all_strings[j + k][i] = 'b';
      }
    }
    block_size = 2 * block_size;
  }
  for (auto i = 0; i < num_rows; ++i) {
    for (auto j = 0; j < num_col; ++j) {
      std::cout << all_strings[i][j] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "Starting building LUF\n";

  for (int str_len = 0; str_len < max_len; ++str_len) {
    auto len = str_len;
    auto num_str = 1 << str_len - 1; // 2^(str_len-1)
    for (auto row = 0; row < num_str; ++row) {
      // Generate string
      std::string text(all_strings[row].data(), str_len);
      std::cout << "row col " << row << " " << str_len << std::endl;

      // Find LUF
      auto num_errors = test(text);
      if (num_errors > 0) {
        std::cout << "FOUND ERROR :::::::::::::\n";
        ++failed_strings;
        break;
      }
    }
  }
  if (failed_strings == 0) {
    std::cout << "All successfull" << std::endl;
  } else {
    std::cout << "num Failed Strings: " << failed_strings << std::endl;
  }
#endif

  return 0;
}
// Returns number of errors found
int test(std::string &text) {
  std::cout << "================\n" << text << std::endl;

  std::vector<UINT> luf(text.size());
  std::vector<UINT> naiveLuf(text.size());
  Stats stats{};
  Borderless borderless(text, cAlphabet, stats);

  std::clock_t startTime = clock();
  borderless.getLUF(luf);
  std::clock_t endTime = clock();
  borderless.naiveLUF(naiveLuf);
  // testing
  int num_errors = 0;
  int n = text.size();
  int log_n = static_cast<int>(std::ceil(log2(text.size())));
  int n_log_n = n * log_n;
std::cout <<"LUF==== " <<std::endl;
  for (auto i = 0; i < text.size(); ++i) {
    std::cout << luf[i] << " " ;
    if (luf[i] != naiveLuf[i]) {
      std::cout << "(string, i, calc, naive): " << text << " " << i << " "
                << luf[i] << " " << naiveLuf[i] << "\n";
      ++num_errors;
    }
  }
  std::cout <<std::endl;
  std::cout << "\nNumber of errors: " << num_errors << "\n\n" << std::endl;
  if (stats.np > n_log_n) {
    std::cout << "ERROR ################## PUSH " << text << " " << text.size()
              << " " << stats.np << "\n";
  }
  if (stats.nlcp > n + n_log_n) {
    std::cout << "ERROR ################## LCP " << text << " " << text.size()
              << " " << stats.nlcp << "\n";
    ++num_errors;
  }
  if (stats.sq_uses > n) {
    std::cout << "ERROR ################## SQ BASE " << text << " "
              << text.size() << " " << stats.sq_uses << "\n";
    ++num_errors;
  }
  if (stats.s_lcp_base > n) {
    std::cout << "ERROR ################## SUCC LCP BASE " << text << " "
              << text.size() << " " << stats.s_lcp_base << "\n";
    ++num_errors;
  }
  if (stats.s_lcp_non_base > n_log_n) {
    std::cout << "ERROR ################## SUCC LCP NON-BASE " << text << " "
              << text.size() << " " << stats.s_lcp_non_base << "\n";
    ++num_errors;
  }
  if (stats.f_lcp_base > n) {
    std::cout << "ERROR ################## FAIL LCP BASE " << text << " "
              << text.size() << " " << stats.f_lcp_base << "\n";
    ++num_errors;
  }
  if (stats.f_lcp_non_base > n_log_n) {
    std::cout << "ERROR ################## FAIL LCP NON BASE " << text << " "
              << text.size() << " " << stats.f_lcp_non_base << "\n";
    ++num_errors;
  }
  if (stats.depth > log_n) {
    std::cout << "ERROR ################## FAIL DEPTH " << text << " "
              << text.size() << " " << stats.depth << "\n";
    ++num_errors;
  }
  statfile << text.size() << "\t" << stats.depth << "\t" << stats.np << "\t"
           << stats.nlcp << "\t" << stats.sq_uses << "\t" << stats.s_lcp_base
           << "\t" << stats.s_lcp_non_base << "\t" << stats.f_lcp_base << "\t"
           << stats.f_lcp_non_base << "\n";

  return num_errors;
}
