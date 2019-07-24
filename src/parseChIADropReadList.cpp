/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include <boost/program_options.hpp>
#include "SSP/common/BedFormat.hpp"

using namespace std;
using namespace boost::program_options;

variables_map argv_init(int argc, char* argv[])
{
  options_description allopts("Options");
  allopts.add_options()
    ("bed,b",      value<string>(), "Bed file")
    ("help,h", "print this message")
    ;

  variables_map values;

  if (argc==1) {
    cout << "\n" << allopts << endl;
    exit(0);
  }
  try {
    parsed_options parsed = parse_command_line(argc, argv, allopts);
    store(parsed, values);

    if (values.count("help")) {
      cout << "\n" << allopts << endl;
      exit(0);
    }
    for (auto x: {"bed"}) {
      if (!values.count(x)) {
	cerr << "specify --" << x << " option." << endl;
	exit(0);
      }
    }
    notify(values);

  } catch (exception &e) {
    cout << e.what() << endl;
    exit(0);
  }
  return values;
}

template <class T>
std::vector<T> parse_readlist(const std::string &fileName)
{
  std::vector<T> vbed;
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

    vbed.emplace_back(v);
  }

  return vbed;
}

int main(int argc, char* argv[])
{
  variables_map values = argv_init(argc, argv);

  string filename(values["bed"].as<string>());
  std::cout << filename << std::endl;
  std::vector<bed> vbed = parse_readlist<bed>(filename);

  return 0;
}
