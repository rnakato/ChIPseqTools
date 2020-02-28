#include <string>
#include <iostream>
#include <algorithm>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include "SSP/common/BedFormat.hpp"
#include "SSP/common/ReadAnnotation.hpp"
#include "gene_bed.h"
#include "SSP/common/inline.hpp"

using Variables = boost::program_options::variables_map;

Variables argv_init(int argc, char* argv[])
{
  boost::program_options::options_description allopts("Options");
  allopts.add_options()
    ("bed1",      boost::program_options::value<std::string>(), "1st bed file")
    ("bed2",      boost::program_options::value<std::string>(), "2nd bed file")
    ("loop",      boost::program_options::value<std::string>(), "Loop file")
    ("length,l",  boost::program_options::value<int32_t>()->default_value(0), "Extend length for overlap ")
    ("gt",        boost::program_options::value<std::string>(), "Genome table (tab-delimited file describing the name and length of each chromosome)")
    ("nobs", "do not output the overlapped loop list")
    ("hiccups",   "HICCups format as input (default: Mango)")
    ("help,h", "print this message")
    ;

  Variables values;

  if (argc==1) {
    std::cout << "\nUsage: compare_bed2loop [option] --bed1 <1st bed> --bed2 <2nd bed> --loop <loop file> -o <output> -gt <genome_table>\n"
	      << allopts
	      << std::endl;
    exit(0);
  }
  try {
    boost::program_options::parsed_options parsed = parse_command_line(argc, argv, allopts);
    store(parsed, values);

    if (values.count("help")) {
      std::cout << "\nUsage: compare_bed2loop [option] --bed1 <1st bed> --bed2 <2nd bed> --loop <loop file> -o <output> -gt <genome_table>\n"
		<< allopts
		<< std::endl;
      exit(0);
    }
    std::vector<std::string> opts = {};
    for (auto x: {"bed1", "bed2", "loop", "gt"}) {
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

int main(int argc, char* argv[])
{
  Variables values = argv_init(argc, argv);

  std::string bed1file(values["bed1"].as<std::string>());
  std::string bed2file(values["bed2"].as<std::string>());
  std::string loopfile(values["loop"].as<std::string>());

  auto vbed1 = parseBed<bed>(bed1file);
  auto vbed2 = parseBed<bed>(bed2file);

  std::string tool("mango");
  if(values.count("hiccups")) tool = "hiccups";
  InteractionSet interset(loopfile, "", tool);

  std::cout << "# Bed1: " << bed1file << std::endl;
  std::cout << "# Number: " << vbed1.size() << std::endl;
  std::cout << "# Bed2: " << bed2file << std::endl;
  std::cout << "# Number: " << vbed2.size() << std::endl;
  std::cout << "# Loop file: " << loopfile << std::endl;

  interset.compare_bed_loop(vbed1, vbed2, values.count("nobs"));

  return 0;
}
