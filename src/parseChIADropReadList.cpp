/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include <unordered_map>
#include "SSP/common/util.hpp"

class posi {
 public:
  std::string chr;
  int32_t start;
 posi(): start(0) {}
  virtual ~posi(){}
  posi(const std::string &c, const std::string &s):
    chr(rmchr(c)), start(stoi(s))
  {}
  bool operator<(const posi &another) const
  {
    if (compare_chr(chr, another.chr) < 0) return 1;
    else if (compare_chr(chr, another.chr) == 0 && start < another.start) return 1;
    else return 0;
  };
};

class readpair {
 public:
  std::string barcode;
//  std::string chr1, chr2;
  int32_t p1, p2;
 readpair(): p1(0), p2(0) {}
  virtual ~readpair(){}
  readpair(const std::string &_barcode, const posi &posi1, const posi &posi2):
//    barcode(_barcode), chr1(posi1.chr), chr2(posi2.chr), p1(posi1.start), p2(posi2.start)
    barcode(_barcode), p1(posi1.start), p2(posi2.start)
  {}

  bool operator<(const readpair &another) const
  {
    if (p1 < another.p1) return 1;
    else if (p1 == another.p1 && p2 < another.p2) return 1;
    else return 0;
  };
/*  bool operator<(const readpair &another) const
  {
    if (compare_chr(chr1, another.chr1) < 0) return 1;
    else if (compare_chr(chr2, another.chr2) < 0) return 1;
    else if (compare_chr(chr1, another.chr1) > 0 || compare_chr(chr1, another.chr1) > 0) return 0;
    else if (p1 < another.p1) return 1;
    else return 0;
  };*/
};

using Variables = boost::program_options::variables_map;
using vBedMap = std::unordered_map<std::string, std::vector<posi>>;
using ChrPair = std::unordered_map<std::string, std::vector<readpair>>;
using MpPair = std::unordered_map<std::string, ChrPair>;

Variables argv_init(int argc, char* argv[])
{
  boost::program_options::options_description allopts("Options");
  allopts.add_options()
    ("bed,b",  boost::program_options::value<std::string>(), "Bed file")
    ("intra",  "output intrachromosomal pairs only")
    ("help,h", "print this message")
    ;

  Variables values;

  if (argc==1) {
    std::cout << "\n" << allopts << std::endl;
    exit(0);
  }
  try {
    boost::program_options::parsed_options parsed = parse_command_line(argc, argv, allopts);
    store(parsed, values);

    if (values.count("help")) {
      std::cout << "\n" << allopts << std::endl;
      exit(0);
    }
    for (auto x: {"bed"}) {
      if (!values.count(x)) {
	std::cerr << "specify --" << x << " option." << std::endl;
	exit(0);
      }
    }
    notify(values);

  } catch (std::exception &e) {
    std::cout << e.what() << std::endl;
    exit(0);
  }
  return values;
}

vBedMap parse_readlist(const std::string &fileName)
{
  vBedMap mp;
  std::ifstream in(fileName);
  if(!in) {
    std::cerr << "Error: BED file " << fileName << " does not exist." << std::endl;
    std::exit(1);
  }

  while (!in.eof()) {
    std::string lineStr;
    getline(in, lineStr);

    if(lineStr.empty() || lineStr[0] == '#') continue;
    std::vector<std::string> v;
    ParseLine(v, lineStr, ',');

    mp[v[0]].emplace_back(v[1], v[2]);
  }

  return mp;
}


int main(int argc, char* argv[])
{
  Variables values = argv_init(argc, argv);

  std::string filename(values["bed"].as<std::string>());
  vBedMap mp = parse_readlist(filename);

  MpPair mppair;

  for(auto &pair: mp) {
    int32_t nbed(pair.second.size());
    if(nbed == 1) continue;

    // sort posi vector
    std::sort(pair.second.begin(), pair.second.end());

    for (int32_t i=0; i<nbed; ++i) {
      for (int32_t j=i+1; j<nbed; ++j) {
	if (values.count("intra") && pair.second[i].chr != pair.second[j].chr) continue;
//	vpair.emplace_back(pair.first, pair.second[i], pair.second[j]);
	mppair[pair.second[i].chr][pair.second[j].chr].emplace_back(pair.first, pair.second[i], pair.second[j]);
      }
    }
  }

  std::vector<std::string> chrlist = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "M"};

  for(auto &x: chrlist) {
    if (mppair.count(x) == 0) continue;
    for(auto &y: chrlist) {
      if (mppair.at(x).count(y) == 0) continue;
      std::vector<readpair> vp(mppair.at(x).at(y));
      std::sort(vp.begin(), vp.end());
      for(auto &z: vp) {
	printf("%s\tchr%s\t%d\tchr%s\t%d\t.\t.\n",
	       z.barcode.c_str(), x.c_str(), z.p1, y.c_str(), z.p2);
      }
    }
  }

  return 0;
}
