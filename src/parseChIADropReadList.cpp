/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include <boost/program_options.hpp>
#include <unordered_map>
#include "SSP/common/BedFormat.hpp"

class posi {
 public:
  std::string chr;
  int32_t start;
 posi(): start(0) {}
  virtual ~posi(){}
  posi(const std::string &c, const std::string &s):
    chr(c), start(stoi(s))
  {}
};

using Variables = boost::program_options::variables_map;
using vBedMap = std::unordered_map<std::string, std::vector<posi>>;

Variables argv_init(int argc, char* argv[])
{
  boost::program_options::options_description allopts("Options");
  allopts.add_options()
    ("bed,b",  boost::program_options::value<std::string>(), "Bed file")
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

  for(auto &pair: mp) {
    int32_t nbed(pair.second.size());
    if(nbed == 1) continue;
    for (int32_t i=0; i<nbed; ++i) {
      for (int32_t j=i+1; j<nbed; ++j) {
	std::cout << pair.first << "\t"
		  << pair.second[i].chr   << "\t"
		  << pair.second[i].start << "\t"
		  << pair.second[j].chr   << "\t"
		  << pair.second[j].start << std::endl;
      }
    }
  }

  return 0;
}
