/******************************************************************************************
  Jesse Waite
  CS516
  2/9/2016

  A NeedlemanWunsch and SmithWaterman sequence nucleotide alignment calculator.

  Usage:
      ./[executable] [local fasta file] [0 or 1] [optional .config file]

  The fasta file is expected to contain two gene sequences in the fasta file; so basically
  two fasta files concatenated. Select 0 to un global/NeedlemanWunsh, or 1 for local/SmithWaterman
  alignment.

  This code was written for correctness, not space efficiency or other code compactness.

  Use 'sh compile.sh' to compile (requires c++11, 02 level optimizations), or run g++ directly as:
    g++ -o seqcmp main.cpp Source.cpp -std=c++0x -O2
*******************************************************************************************/

#include "SequenceAlignment.hpp"

int main(int argc, char** argv)
{
    string fasta, config;
    Sequence seq1, seq2;
    bool isGlobal;
    Params params;
    Alignment alignment;

    if (argc < 3) {  // <executable name> <input file containing both s1 and s2> <0: global, 1: local> <optional: path to parameters config file>
        cout << "Usage: ./seqcmp [fasta file] [0 (global) or 1 (local)] [optional: path to parameter file]" << endl;
        return 1;
    }

    //get path to fasta file, containing two sequences
    fasta = argv[1];
    if(!fileExists(fasta)){
       cout << "ERROR file >" << fasta << "< not found or not accessible" << endl;
      return 1;
    }

    //get the flag for NeedlemanWunsch (global) or SmithWaterman (local)
    isGlobal = argv[2][0] == '0';
    if(argc == 4){
      config = argv[3];
      if(!fileExists(config)){
        cout << "ERROR file >" << config << "< not found or not accessible" << endl;
        return 1;
      }

      //populate the parameter object
      SequenceAlignment::ParseParamsFile(config, params);
    }
    else{
      config = "default";
      params.match = 1;
      params.mismatch = -2;
      params.g = -2;
      params.h = -5;
    }

    //read the sequences
    parse2FastaFile(fasta, seq1, seq2);

    SequenceAlignment* stringComp = new SequenceAlignment();

    if(isGlobal){
      stringComp->NeedlemanWunsch(seq1.seq,seq2.seq,params,alignment);
      stringComp->PrintResult(seq1, seq2, params, alignment);
    }
    else{
      stringComp->SmithWaterman(seq1.seq, seq2.seq, params, alignment);
      stringComp->PrintResult(seq1, seq2, params, alignment);
    } 
   
    delete stringComp;

    return 0;
}
