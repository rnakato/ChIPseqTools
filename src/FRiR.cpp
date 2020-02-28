/* Copyright(c) Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * All rights reserved.
 */
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <boost/filesystem.hpp>
#include "DROMPAplus/src/extendBedFormat.hpp"
#include "DROMPAplus/submodules/SSP/common/BoostOptions.hpp"
#include "DROMPAplus/submodules/SSP/common/inline.hpp"
#include "DROMPAplus/submodules/SSP/src/ParseMapfile.hpp"
#include "DROMPAplus/submodules/SSP/src/Mapfile.hpp"
#include "DROMPAplus/submodules/SSP/src/LibraryComplexity.hpp"

class globalopt{
  MyOpt::Opts opt;
  std::string repeatfilename;
  std::string repeattype;
  std::string output;

public:
  SeqStatsGenome genome;
  LibComp complexity;

  globalopt(): opt("Repeats",100), repeatfilename("") {
    opt.add_options()
      ("repeat,r",   boost::program_options::value<std::string>(), "Repeat file")
      ("repeattype", boost::program_options::value<std::string>()->default_value("class"), "Repeat type: [class|family|name]")
      ;
  }

  void setOpts(MyOpt::Opts &allopts) {
    genome.setOpts(allopts);
    allopts.add(opt);
  }
  void setValues(const MyOpt::Variables &values) {
    genome.setValues(values);
    repeatfilename = MyOpt::getVal<std::string>(values, "repeat");
    repeattype     = MyOpt::getVal<std::string>(values, "repeattype");
    output         = MyOpt::getVal<std::string>(values, "output");
  }
  const std::string & getrepeatfilename() const { return repeatfilename; }
  const std::string & getrepeattype()     const { return repeattype; }
  const std::string & getoutput()         const { return output; }
};

void getOpts(globalopt &p, int argc, char* argv[]);

void help_global()
{
  auto helpmsg = R"(
===============

Usage: FRiR [option] -r <repeatfile> -i <inputfile> -o <output> --gt <genome_table>)";

  std::cerr << "\nFRiR" << helpmsg << std::endl;
  return;
}

using mapRep = std::unordered_map<std::string, std::vector<bed>>;

mapRep read_Repeatfile(const std::string &filename, const std::string &type)
{
  std::ifstream in(filename);
  if(!in) {
    std::cerr << "Error: Repeat file does not exist." << std::endl;
    std::exit(1);
  }

  mapRep Repeat;
  std::string lineStr;
  std::vector<std::string> v;
  while (!in.eof()) {
    getline(in, lineStr);

    if(lineStr.empty() || lineStr[0] == '#') continue;
    boost::split(v, lineStr, boost::algorithm::is_any_of("\t"));

    std::string repName(v[10]);
    std::string repClass(v[11]);
    std::string repFamily(v[12]);
    if(type == "class")    Repeat[repClass].emplace_back(v[5], stoi(v[6]), stoi(v[7]));
    else if(type == "family") Repeat[repFamily].emplace_back(v[5], stoi(v[6]), stoi(v[7]));
    else if(type == "name")   Repeat[repName].emplace_back(v[5], stoi(v[6]), stoi(v[7]));
    else {
      std::cerr << "Error: invalid repeattype.\n";
      exit(0);
    }
  }

  return Repeat;
}

uint64_t countReadsChr(const std::vector<bed> &vbed, const SeqStats &chr, uint64_t &len)
{
  uint64_t nread(0);
  std::vector<BpStatus> array(chr.getlen(), BpStatus::MAPPABLE);
  OverrideBedToArray(array, chr.getname(), vbed);

  len += std::count(array.begin(), array.end(), BpStatus::INBED);

  for (auto strand: {Strand::FWD, Strand::REV}) {
    for (auto &x: chr.getvReadref(strand)) {
      //      if(x.duplicate) continue;
      int32_t s(std::min(x.F3, x.F5));
      int32_t e(std::max(x.F3, x.F5));
      for(int32_t i=s; i<=e; ++i) {
	if(array[i] == BpStatus::INBED) {
	  ++nread;
	  break;
	}
      }
    }
  }

  return nread;
}


uint64_t countReads(const std::vector<bed> &vbed, const SeqStatsGenome &genome, uint64_t &len)
{
  uint64_t nread(0);
  for(auto &chr: genome.chr) nread += countReadsChr(vbed, chr, len);
  return nread;
}

int main(int argc, char* argv[])
{
  globalopt p;
  getOpts(p, argc, argv);

  read_mapfile(p.genome);
  p.complexity.checkRedundantReads(p.genome);

  auto Repeat = read_Repeatfile(p.getrepeatfilename(), p.getrepeattype());

  std::ofstream out(p.getoutput());

  out << "Type\ttotal length\t% genome\tread number\t% genome\tlog10(read/length)" << std::endl;
  for(auto itr = Repeat.begin(); itr != Repeat.end(); ++itr) {
    uint64_t len(0);
    auto nread = countReads(itr->second, p.genome, len);
    double rread(getratio(nread, p.genome.getnread(Strand::BOTH)));
    double rlen(getratio(len, p.genome.getlen()));
    out << itr->first << "\t"
	<< len   << "\t" << rlen  << "\t"
	<< nread << "\t" << rread << "\t"
	<< log10(getratio(rread, rlen)) << std::endl;
	}

  return 0;
}

void getOpts(globalopt &p, int argc, char* argv[])
{
  DEBUGprint("setOpts...");

  MyOpt::Opts allopts("Options");
  MyOpt::setOptIO(allopts, "sspout");
  MyOpt::setOptPair(allopts);
  p.setOpts(allopts);
  MyOpt::setOptOther(allopts);

  DEBUGprint("getOpts...");

  MyOpt::Variables values;

  try {
    boost::program_options::parsed_options parsed = parse_command_line(argc, argv, allopts);
    store(parsed, values);
  }
  catch(const boost::program_options::error_with_option_name& e) {
    std::cout << e.what() << std::endl;
    exit(0);
  }
  if (argc ==1) {
    help_global();
    std::cerr << "Use --help option for more information on the other options\n\n";
    exit(0);
  }
  if (values.count("help")) {
    help_global();
    std::cout << "\n" << allopts << std::endl;
    exit(0);
  }
  std::vector<std::string> opts = {"input", "output", "gt"};
  for (auto x: opts) {
    if (!values.count(x)) PRINTERR_AND_EXIT("specify --" << x << " option.");
  }

  try {
    notify(values);
    p.setValues(values);

    //    init_dump(values, ssp);
  } catch(const boost::bad_any_cast& e) {
    std::cout << e.what() << std::endl;
    exit(0);
  }

  DEBUGprint("getOpts done.");
  return;
}
