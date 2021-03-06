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

/** Defines some helper functions */

#include "../include/utilDefs.hpp"

namespace luf {

static struct option long_options[] = {
    {"mode", required_argument, NULL, 'm'},
    {"alphabet", required_argument, NULL, 'a'},
    {"selected-alphabet", optional_argument, NULL, 's'},
    {"input-file", required_argument, NULL, 'i'},
    {"antidictionary-file", required_argument, NULL, 'd'},
    {"help", no_argument, NULL, 'h'},
    {NULL, 0, NULL, 0}};

/** Decode the input flags
 */
ReturnStatus decodeFlags(int argc, char *argv[], struct InputFlags &flags) {
  int args = 0;
  int opt;
  std::string alph;
  std::string mode;

  /* initialisation */
  while ((opt = getopt_long(argc, argv, "m:a:s:i:d:h", long_options,
                            nullptr)) != -1) {
    switch (opt) {
    case 's':
      flags.selected_alphabet = std::string(optarg);
      args++;
      break;

    case 'i':
      flags.input_filename = std::string(optarg);
      args++;
      break;

    case 'd':
      flags.anti_dictionary_filename = std::string(optarg);
      args++;
      break;

    case 'h':
      return (ReturnStatus::HELP);
    }
  }
  if (args < 4) {
    std::cerr << "Invalid command: Too few arguments: " << std::endl;
    return (ReturnStatus::ERR_ARGS);
  }  else {
    return (ReturnStatus::SUCCESS);
  }
}

/*
 * Usage of the tool
 */
void usage(void) {
  std::cout << " Usage: mawcd <options>\n";
  std::cout << " Standard (Mandatory):\n";
  std::cout << "  -m, --mode \t\t\t<str> \t\t`AD' for anti-dictionary creation "
               "\n\t\t\t\t\t\t or `COM' for compression of single file "
               "\n\t\t\t\t\t\t or `DECOM' for decompression of single file "
               "\n\t\t\t\t\t\t or `BCOM' for compression of many files "
               "\n\t\t\t\t\t\t or `BDECOM' for decompression of many files. \n\n";
  std::cout << "  -a, --alphabet \t\t <str> \t \t `DNA' for nucleotide sequences"
               "\n\t\t\t\t\t\t or `PROT' for protein  sequences "
               "\n\t\t\t\t\t\t or `SEL' for user-defined "
               "\n\t\t\t\t\t\t or `GEN' for general (ASCII). \n\n";
  std::cout
      << "  -i, --input-file \t\t <str> \t \t Input file  name  "
         "\n\t\t\t\t\t\t(uncompressed file when mode is `COM'; compressed file when mode is `DECOM'; "
         "\n\t\t\t\t\t\t a file containing names of the files to compress or decompress in batch mode [one name on each line]).\n\n";

  std::cout << "  -d, --antidictionary-file \t <str> \t \t Anti-dictionary file  name "
               "\n\t\t\t\t\t\t(created when mode is `AD' and read when mode is any other).\n\n";

  std::cout << " Additional:\n";
  std::cout
      << "  -s, --selected-alphabet \t <str> \t \t case-sensitive alphabet  "
         "(required when alphabet is SEL). \n\n";
}

} // end namespace
