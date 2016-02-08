#pragma once
#include <stdio.h>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <limits.h>

using namespace std;

typedef enum { DEL, INS, SUB } CellFlag;
//typedef enum { GLOBAL, SEMIGLOBAL, LOCAL } AlignmentStrategy;

#define NEG_INF -1000000000  //from limits.h

typedef struct DpCell {
    int deletionScore;
    int insertionScore;
    int substitutionScore;
    //CellFlag flag;
    ///... // add any other field(s) that you may need for the implementation
}Cell;

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
};

class SequenceComparer
{
    private:
        vector<vector<Cell> > _dpTable;
        void _clearTable();
        void _resize(int rows, int cols);
        int _maxThree(int a, int b, int c);
        void _reportProgress(int row, int numRows);
        void _scoreAffine_NeedlemanWunsch(char a, char b, int row, int col, vector<vector<DpCell> >& dpTable, const Params& params);
        void _scoreAffine_SmithWaterman(char a, char b, int row, int col, vector<vector<DpCell> >& dpTable, const Params& params);
        void _printTable(const vector<vector<DpCell> >& dpTable);
        int _getAffineDeletionScore(const DpCell& predecessor, const Params& params);
        int _getAffineInsertionScore(const DpCell& predecessor, const Params& params);
        int _getAffineSubstitutionScore(bool isMatch, const DpCell& predecessor, const Params& params);

    public:
        SequenceComparer();
        ~SequenceComparer();

        void PrintResult(const string& seq1, const string& seq2, const Params& params, const Alignment& alignment);
        void NeedlemanWunsch(const string& seq1, const string& seq2, Params& params, Alignment& alignment);
        void SmithWaterman(const string& seq1, const string& seq2, Params& params, Alignment& alignment);
        //needleman-wunsch
        //int GetLevenshteinDist(string s1, string s2, );
        //string GetLCS(string s1, string s2);
};

void ParseFastaFile(const string fname, string& s1, string& s2);
void ParseParamsFile(const string& fname, Params& params);












