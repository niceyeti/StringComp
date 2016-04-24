#pragma once

#ifndef _INC_STDIO
#include <stdio.h>
#endif

#ifndef _STRING_
#include <string>
#endif

#ifndef _CSTRING_
#include <cstring>
#endif

#ifndef _IOSTREAM_
#include <iostream>
#endif

#ifndef _FSTREAM_
#include <fstream>
#endif

#ifndef _VECTOR_
#include <vector>
#endif

#ifndef _ALGORITHM_
#include <algorithm>
#endif

#ifndef _INC_LIMITS
#include <limits.h>
#endif

#include "..\..\mccreight\mccreight\Util.hpp"

using namespace std;

//consider using this to cut down on matrix size, using smallest native int possible, capable of holding max/min possible score.
//could also split score types since SW only needs unsigned positive values.
typedef int SCORE;
typedef enum { DEL, INS, SUB } CellFlag;
//typedef enum { GLOBAL, SEMIGLOBAL, LOCAL } AlignmentStrategy;

#define NEG_INF -10000000

typedef struct DpCell {
    SCORE deletionScore;
    SCORE insertionScore;
    SCORE substitutionScore;
    //CellFlag flag;

    //custom for prg3; these are terribly space inefficient for very large alignments, but fine for small ones
    //int matches;
    //int alignmentLen;

    ///... // add any other field(s) that you may need for the implementation
}Cell;

//For optimized implementations, there's no backtracking; the alignment matches/length are propagated in the forward step.
//The additional score data increases the table size, so there is a conscientious time/space tradeoff.
typedef struct optimizedScore {
    int score;
    int matches;
    int alignmentLength;
}OptimizedScore;

typedef struct optimizedCell {
    OptimizedScore deletionScore;
    OptimizedScore insertionScore;
    OptimizedScore substitutionScore;
}OptimizedCell;

//For time-optimized alignments, we're only intereted in numerical valus of #matches and overall alignment length. This
//will typically only be used for optimized SmithWaterman.
typedef struct optimizedAlignment {
    int length;
    int numMatches;
}OptimizedAlignment;

//Input parameters object to scoring metrics
typedef struct Parameters {
    int mismatch;
    int match;
    int h;
    int g;
}Params;


class Alignment {
public:

    Alignment()
    {
        Clear();
    }

    string method;
    int maxScore;
    int gaps;
    int matches;
    int mismatches;
    int openingGaps;
    string s1;
    string s2;
    string bridge;
    void Clear();
    void Print();
    void PrintValidity(const Parameters& params) const;
    int Length();
};

//TODO: factor optimized and non-optimized aligners into different objects, and make a common wrapper class
class SequenceAlignment
{
    private:
        bool _verbose;
        vector<vector<Cell> > _dpTable;
        vector<vector<OptimizedCell> > _optimizedTable;
        void _clearTable();
        void _resize(int rows, int cols);
        int _maxThree(int a, int b, int c);
        inline void _reportProgress(int row, int numRows);
        void _scoreAffine_NeedlemanWunsch(char a, char b, int row, int col, vector<vector<DpCell> >& dpTable, const Params& params);
        void _scoreAffine_SmithWaterman(char a, char b, int row, int col, vector<vector<DpCell> >& dpTable, const Params& params);
        void _scoreAffine_SmithWaterman_Opt(char a, char b, int row, int col, vector<vector<OptimizedCell> >& optimizedTable, const Params& params);
        void _printTable(const vector<vector<DpCell> >& dpTable);
        int _getMaxDirection(const Cell& predecessor);
        int _getAffineDeletionScore(const DpCell& predecessor, const Params& params);
        int _getAffineInsertionScore(const DpCell& predecessor, const Params& params);
        int _getAffineSubstitutionScore(bool isMatch, const DpCell& predecessor, const Params& params);
        //optimized scoring
        void _setAffineDeletionScore_Opt(const OptimizedCell& predecessor, OptimizedCell& cell, const Params& params);
        void _setAffineInsertionScore_Opt(const OptimizedCell& predecessor, OptimizedCell& cell, const Params& params);
        void _setAffineSubstitutionScore_Opt(bool isMatch, const OptimizedCell& predecessor, OptimizedCell& cell, const Params& params);
        void _updateAlignment(const int curState, const int prevState, Alignment& alignment, const bool isMatch, const int i, const int j, const string& seq1, const string& seq2);
        int _getPrevState(const int curState, const Cell& cell, const Params& params, const bool isMatch);
        bool _hasPositiveScore(const struct DpCell& cell);
        void _verboseUpdate(const int curState, const int prevState, Alignment& alignment, const bool isMatch, const int i, const int j, const string& seq1, const string& seq2);
        void _resizeInPlace_Optimized(int rows, int cols);
        void _clearOptimizedScore(OptimizedScore& optScore);
public:
        SequenceAlignment();
        ~SequenceAlignment();
        static void ParseParamsFile(const string& fname, Params& params);
        void ClearOptimizedTable();
        int GetOptimizedTableSize();
        void PrintResult(const Sequence& sequence1, const Sequence& sequence2, const Params& params, const Alignment& alignment);
        void NeedlemanWunsch(const string& seq1, const string& seq2, Params& params, Alignment& alignment, bool verbose=true);
        void SmithWaterman(const string& seq1, const string& seq2, Params& params, Alignment& alignment, bool verbose=true);
        void SmithWaterman_Optimized(int i1, int n1, const string& seq1, int i2, int n2, const string& seq2, Params& params, OptimizedAlignment& alignment);
        //int GetLevenshteinDist(string s1, string s2, );
        //string GetLCS(string s1, string s2);  <-- use a suffix tree, not an alignment
};

//void ParseFastaFile(const string fname, Sequence& s1, Sequence& s2);
//void ParseParamsFile(const string& fname, Params& params);
//bool fileExists(const string& path);
