#include <string>
#include <iostream>
#include <algorithm>
#include <boost/program_options.hpp>
#include "header_drompaplus.hpp"

using namespace std;
using namespace boost::program_options;

variables_map argv_init(int argc, char* argv[])
{
  options_description allopts("Options");
  allopts.add_options()
    ("gtf,g", value<string>(), "Gene annotation file (gtf format)")
//    ("name,n", "Output name instead of id")
    ("unique,u", "Output one reference transcript per one gene (default: all transcripts)")
    ("help,h", "Print this message")
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
    vector<string> opts = {};
    for (auto x: {"gtf"}) {
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

void printRefFlat(const HashOfGeneDataMap &mp, const int32_t unique)
{
  std::cout << "#geneName\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\tgene type\ttranscript type";
  if (unique) std::cout << "\t" << "reference transcript name\treference transcript id";
  else        std::cout << "\t" << "gene name\tgene id";
  std::cout << std::endl;

  for (auto pair: mp) {
    for (auto x: pair.second) {
      if (unique) std::cout << x.second.gname << "\t" << x.second.gid << "\t";
      else        std::cout << x.second.tname << "\t" << x.second.tid << "\t";

      std::cout << "chr"<< pair.first << "\t"
                << x.second.strand << "\t"
                << (x.second.txStart -1) << "\t"
                << x.second.txEnd << "\t";

      if (x.second.cdsStart) {
        std::cout << (x.second.cdsStart -1) << "\t" << x.second.cdsEnd << "\t";
      } else {
        std::cout << (x.second.txStart -1) << "\t" << x.second.txEnd  << "\t";
      }

      std::cout << x.second.exonCount << "\t";
      if (x.second.strand == "+") {
        for (int32_t i=0; i<x.second.exonCount; ++i) std::cout << (x.second.exon[i].start -1) << ",";
        std::cout << "\t";
        for (int32_t i=0; i<x.second.exonCount; ++i) std::cout << x.second.exon[i].end << ",";
      } else {
        for (int32_t i=x.second.exonCount-1; i>=0; --i) std::cout << (x.second.exon[i].start -1) << ",";
        std::cout << "\t";
        for (int32_t i=x.second.exonCount-1; i>=0; --i) std::cout << x.second.exon[i].end << ",";
      }

      std::cout << "\t" << x.second.gtype << "\t" << x.second.ttype;

      if (unique) std::cout << "\t" << x.second.tname << "\t" << x.second.tid;
      else        std::cout << "\t" << x.second.gname << "\t" << x.second.gid;

      std::cout << std::endl;
    }
  }
  return;
}

int main(int argc, char* argv[])
{
  variables_map values = argv_init(argc, argv);

  auto tmp = parseGtf(values["gtf"].as<string>()); // hash for transcripts
  auto gmp = construct_gmp(tmp);                   // hash for genes

  //  printMap(tmp);
  //printMap(gmp);

  if (values.count("unique")) printRefFlat(gmp, values.count("unique"));
  else                        printRefFlat(tmp, values.count("unique"));

  return 0;
}
