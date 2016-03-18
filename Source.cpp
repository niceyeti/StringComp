#include "Header.hpp"

void Alignment::Print()
{

}

void Alignment::Clear()
{
    maxScore = 0;
    gaps = 0;
    matches = 0;
    mismatches = 0;
    openingGaps = 0;
    s1.clear();
    s2.clear();
    bridge.clear();
}

SequenceComparer::SequenceComparer()
{}

SequenceComparer::~SequenceComparer()
{
    _clearTable();
}

/*
Resizes the (square) matrix.
*/
void SequenceComparer::_resize(int rows, int cols)
{
    DpCell cell;
    cell.deletionScore = NEG_INF;
    cell.insertionScore = NEG_INF;
    cell.substitutionScore = NEG_INF;

    _clearTable();
    _dpTable.resize(rows);
    for (int i = 0; i < _dpTable.size(); i++) {
        _dpTable[i].resize(cols,cell);
    }
}

void SequenceComparer::PrintResult(const Sequence& sequence1, const Sequence& sequence2, const Params& params, const Alignment& alignment)
{
    int i, j, k, n, s1Ct, s2Ct;
    const string& seq1 = sequence1.seq;
    const string& seq2 = sequence2.seq;

    cout << "INPUT:\r\n******\r\n";

    cout << "\r\n" << sequence1.desc << ">\r\n";
    for (i = 0; i < seq1.size(); i++) {
        cout << seq1[i];
        if (i % 60 == 59)
            cout << "\r\n";
    }

    cout << "\r\n" << sequence2.desc << ">\r\n";
    for (i = 0; i < seq2.size(); i++) {
        cout << seq2[i];
        if (i % 60 == 59)
            cout << "\r\n";
    }
    cout << endl;

    //print the output
    cout << "\r\nOUTPUT:\r\n*******\r\n\r\n";
    cout << "scores:    match = " << params.match << ", mismatch = " << params.mismatch;
    cout << ", h = " << params.h << ", g = " << params.g << "\r\n\r\n";
    
    cout << "Sequence 1 = \"" << sequence1.desc << "\", length = " << seq1.length() << "\r\n";
    cout << "Sequence 2 = \"" << sequence2.desc << "\", length = " << seq2.length() << "\r\n" << endl;

    if (!((alignment.bridge.size() == alignment.s1.size()) && (alignment.bridge.size() == alignment.s2.size()) && alignment.s1.size() == alignment.s2.size())) {
        cout << "ERROR alignments should be of equal size. bridge=" << alignment.bridge.size();
        cout << " s1=" << alignment.s1.size() << " s2=" << alignment.s2.size() << endl;
        //return;
    }

    //print out the alignments, with a field width of 60 chars, and padded by character numbers
    i = 0; j = 0, k = 0; n = 0, s1Ct = 0; s2Ct = 0;
    while ((i < alignment.s1.length()) && (j < alignment.s2.length()) && (k < alignment.bridge.length())) {
        printf("s1: %5d ",(s1Ct+1));
        for (n = 0; n < 60 && i < alignment.s1.length(); n++, i++) {
            cout << alignment.s1[i];
            if (alignment.s1[i] != '-')
                s1Ct++;
        }
        cout << " " << s1Ct << "\r\n";

        //print 60 chars of the bridge sequence
        cout << "          ";
        for (n = 0; n < 60 && k < alignment.bridge.size(); n++, k++) {
            cout << alignment.bridge[k];
        }
        cout << "\r\n";

        //print 60 chars of s2
        printf("s2: %5d ",(s2Ct+1));
        for (n = 0; n < 60 && j < alignment.s2.length(); n++, j++) {
            cout << alignment.s2[j];
            if (alignment.s2[i] != '-')
                s2Ct++;
        }
        cout << " " << s2Ct << "\r\n" << endl;
    }

    cout << endl;

    //Reporting
    int identities = alignment.bridge.size();
    //cout << "alignment lengths (s1,s2,bridge): " << alignment.s1.length() << " " << alignment.s2.length() << " " << alignment.bridge.length() << endl;
    cout << "Report:\r\n\r\n" << alignment.method << " score = " << alignment.maxScore << endl;
    cout << "Number of:   matches = " << alignment.matches << ", mismatches = " << alignment.mismatches;
    cout << ", gaps = " << alignment.gaps << ", opening gaps = " << alignment.openingGaps << endl;
    cout << "Identities: " << alignment.matches << "/" << identities;
    cout << " (" << (int)(((float)alignment.matches / (float)identities) * 100) << "%)";
    cout << ", Gaps = " << alignment.gaps << "/" << identities << " (" << (int)(((float)alignment.gaps / (float)identities) * 100) << "%)";
    cout << endl;
}

/*
Zeroes the matrix.
*/
void SequenceComparer::_clearTable()
{
    for (int i = 0; i < _dpTable.size(); i++) {
        _dpTable[i].clear();
    }
    _dpTable.clear();
}

//Report progress every 50th row; clearly this ignores the cost/variation of row lengths
void SequenceComparer::_reportProgress(int row, int numRows)
{
    if (row % 50 == 0) {
        cout << "\r" << (((float)row / (float)numRows) * 100.0) << "% complete...      " << flush;
    }
}

/*
Implements NeedlemanWunsch algorithm using an affine scoring strategy. The algorithm may use either
local or global alignment.

This implementation follows the matrix formality of having seq1 represented row-wise, seq2 is column-wise.
*/
void SequenceComparer::NeedlemanWunsch(const string& seq1, const string& seq2, Params& params, Alignment& alignment)
{
    int i, j, state, prevState;
    bool isMatch;

    alignment.method = "NeedlemanWunsch";

    //init the matrix
    _resize(seq1.size() + 1, seq2.size() + 1);

    //set 0th row entries as follows: insertion-scores = [0,h+g,+g,+g,+g,+g ... +g]
    _dpTable[0][1].insertionScore = params.g + params.h;
    for (j = 2; j < _dpTable[0].size(); j++) {
        _dpTable[0][j].insertionScore = _dpTable[0][j - 1].insertionScore + params.g;
    }
    //init 0th row deletion and substitution scores to -INF
    for (j = 1; j < _dpTable[0].size(); j++) {
        _dpTable[0][j].substitutionScore = _dpTable[0][j].deletionScore = NEG_INF;
    }
    //init 0th column entries as follows: deletion-scores = [0,h+g,+g,+g,+g,+g ... +g]
    _dpTable[1][0].deletionScore = params.g + params.h;
    for (i = 2; i < _dpTable.size(); i++) {
        _dpTable[i][0].deletionScore = _dpTable[i-1][0].deletionScore + params.g;
    }
    //init 0th column insertion and substitution scores to -INF
    for (i = 1; i < _dpTable.size(); i++) {
        _dpTable[i][0].substitutionScore = _dpTable[i][0].insertionScore = NEG_INF;
    }
    //init the first cell scores to 0 
    _dpTable[0][0].insertionScore = 0;
    _dpTable[0][0].deletionScore = 0;
    _dpTable[0][0].substitutionScore = 0;

    cout << "Running NeedlemanWunsch forward algorithm..." << endl;
    //run forward algorithm
    for (i = 1; i < _dpTable.size(); i++) {
        for (j = 1; j < _dpTable[i].size(); j++) {
            //propagate the scores, by row
            _scoreAffine_NeedlemanWunsch(seq1[i - 1], seq2[j - 1], i, j, _dpTable, params);
        }
        //show progress
        _reportProgress(i, _dpTable.size());
    }

    //_printTable(_dpTable);

    cout << "\r\nBacktracking to find optimal alignment..." << endl;
    //global backtrack: from bottom-right cell to find optimal alignment, and track score
    i = _dpTable.size() - 1;
    j = _dpTable[0].size() - 1;
    alignment.Clear();
    alignment.maxScore = _maxThree(_dpTable[i][j].deletionScore, _dpTable[i][j].insertionScore, _dpTable[i][j].substitutionScore);
    cout << "MAx score IS: " << alignment.maxScore << endl;
    Cell& cell = _dpTable[i][j];
    //get the first maximal state from the bottom/right-most cell
    if (cell.deletionScore >= cell.insertionScore && cell.deletionScore >= cell.substitutionScore) {
      state = DEL;
    }
    else if (cell.insertionScore >= cell.deletionScore && cell.insertionScore >= cell.substitutionScore) {
      state = INS;
    }
    else{ //assume substitution for all other cases
      state = SUB;
    }

    while(i >= 0 && j >= 0){
      isMatch = seq1[i] == seq2[j];
      prevState = state;

      switch(state){
        case DEL:
          //given we're in the DEL state, get state of the max of predecessor scores using affine rules
          state = _getPrevState(DEL, _dpTable[i-1][j], params, isMatch);
          break;

        case INS:
          //given we're in the INS state, get state of the max of predecessor scores using affine rules
          state = _getPrevState(INS, _dpTable[i][j-1], params, isMatch);
          break;

        case SUB:
          //given we're in the SUB state, get state of the max of predecessor scores using affine rules
          state = _getPrevState(SUB, _dpTable[i-1][j-1], params, isMatch);
          break;

        default:
            cout << "UNKNOWN BACKTRACK STATE: " << state << endl;
          break;
      }

      _updateAlignment(state, prevState, alignment, isMatch, i, j, seq1, seq2);

      //update the table indices
      switch(state){
        case DEL: //traverse to previous row
          i--;
        break;
        case INS: //traverse to previous column
          j--;
        break;
        case SUB: //traverse to previous row and previous column
          j--; i--;
        break;
        default:
          cout << "UKNOWN BACKTRACK STATE" << endl;
        break;
      }
    }


/*    while (i > 0 && j > 0) {
        Cell& cell = _dpTable[i][j];

        //deletion score is max, so traverse up
        if (cell.deletionScore >= cell.insertionScore && cell.deletionScore >= cell.substitutionScore) {
            alignment.bridge = " " + alignment.bridge;
            alignment.s1 = seq1[i - 1] + alignment.s1;
            alignment.s2 = "-" + alignment.s2;
            alignment.gaps++;
            //look in upper cell: was this an opening deletion-gap?
            if (previousTraversal != DEL) {
                alignment.openingGaps++;
            }
            previousTraversal = DEL;
            i--;
        }
        //insertion score is max, so traverse left
        else if (cell.insertionScore >= cell.deletionScore && cell.insertionScore >= cell.substitutionScore) {
            alignment.bridge = " " + alignment.bridge;
            alignment.s1 = "-" + alignment.s1;
            alignment.s2 = seq2[j - 1] + alignment.s2;
            alignment.gaps++;
            //was this an opening insertion-gap?
            if (previousTraversal != INS) {
                alignment.openingGaps++;
            }
            previousTraversal = INS;
            j--;
        }
        else{ //assume substitution for all other cases, so traverse diagonally up and left
            alignment.s1 = seq1[i - 1] + alignment.s1;
            alignment.s2 = seq2[j - 1] + alignment.s2;
            if (seq1[i - 1] == seq2[j - 1]) {
                alignment.matches++;
                alignment.bridge = "|" + alignment.bridge;
            }
            else {
                alignment.mismatches++;
                alignment.bridge = " " + alignment.bridge;
            }
            previousTraversal = SUB;
            i--; j--;
        }
    }*/
    
    _clearTable();
}

/*
Given a current backtrack state (DEL,SUB,INS) this applies the affine rules in reverse to find the next direction
(backtrack state) in which to traverse.

curState: The method by which this cell was entered: DEL, INS, or SUB.
cell: The cell just entered.
params: The scoring parameters
isMatch: whether or not the characters in seq1 and seq2 match at this cell
*/
int SequenceComparer::_getPrevState(const int curState, const Cell& cell, const Params& params, const bool isMatch)
{
  int sub, del, ins;
  int previous;

  switch(curState){
    case SUB:
      sub = isMatch ? (cell.substitutionScore + params.match) : (cell.substitutionScore + params.mismatch);
      del = isMatch ? (cell.deletionScore + params.match) : (cell.deletionScore + params.mismatch);
      ins = isMatch ? (cell.insertionScore + params.match) : (cell.insertionScore + params.mismatch);
    break;

    case DEL:
      sub = cell.substitutionScore + params.h + params.g;
      ins = cell.insertionScore + params.h + params.g;
      del = cell.deletionScore + params.g;
    break;

    case INS:
      sub = cell.substitutionScore + params.h + params.g;
      ins = cell.insertionScore + params.g;
      del = cell.deletionScore + params.h + params.g;
    break;
  }

  if(ins >= del && ins >= sub){
    previous = INS;
  }
  else if(del >= ins && del >= sub){
    previous = DEL;
  }
  else{
    previous = SUB;
  }

  return previous;
}

/*
Updates the alignment and score, while backtracking.

curState: The method by which we entered the current max cell (INS, DEL, or SUB)
prevState: The pethod by which we entered the previous max cell (INS, DEL, or SUB); is used to detect opening gaps by edge-detecting wrt curState.
alignment: The alignment and score info
isMatch: Whether or not the current chars in the current max cell match one another
i: row of current cell, and ith index of seq1
j: col of current cell, and jth index of seq2
*/
void SequenceComparer::_updateAlignment(const int curState, const int prevState, Alignment& alignment, const bool isMatch, const int i, const int j, const string& seq1, const string& seq2)
{
  if(curState == DEL){  //deletion
    alignment.gaps++;
    if(prevState != DEL){
      alignment.openingGaps++;
    }
    alignment.s1 = seq1[i] + alignment.s1;
    alignment.s2 = "-" + alignment.s2;
    alignment.bridge = " " + alignment.bridge;
  }
  else if(curState == INS){  //insertion
    alignment.gaps++;
    if(prevState != INS){
      alignment.openingGaps++;
    }
    alignment.s1 = "-" + alignment.s1;
    alignment.s2 = seq2[j] + alignment.s2;
    alignment.bridge = " " + alignment.bridge;
  }
  else if(curState == SUB){  //substitution
    if(isMatch){
      alignment.matches++;
      alignment.s1 = seq1[i] + alignment.s1;
      alignment.s2 = seq2[j] + alignment.s2;
      alignment.bridge = "|" + alignment.bridge;
    }
    else{
      alignment.mismatches++;
      alignment.s1 = seq1[i] + alignment.s1;
      alignment.s2 = seq2[j] + alignment.s2;
      alignment.bridge = " " + alignment.bridge;
    }
  }
  else{
    cout << "ERROR unmapped state in updateAlignment" << endl;
  }
}

void SequenceComparer::SmithWaterman(const string& seq1, const string& seq2, Params& params, Alignment& alignment)
{
    int i, j;
    int maxScore;
    pair<int, int> maxIndices;

    alignment.method = "SmithWaterman";

    //init the matrix
    _resize(seq1.size() + 1, seq2.size() + 1);
    //init all values in first column to zero
    for (j = 0; j < _dpTable[0].size(); j++) {
        _dpTable[0][j].deletionScore = 0;
        _dpTable[0][j].insertionScore = 0;
        _dpTable[0][j].substitutionScore = 0;
    }
    //init all vals in first column to zero
    for (i = 0; i < _dpTable.size(); i++) {
        _dpTable[i][0].deletionScore = 0;
        _dpTable[i][0].insertionScore = 0;
        _dpTable[i][0].substitutionScore = 0;
    }

    cout << "Running SmithWaterman forward algorithm..." << endl;
    //run forward algorithm
    maxScore = NEG_INF;
    for (i = 1; i < _dpTable.size(); i++) {
        for (j = 1; j < _dpTable[i].size(); j++) {
            //propagate the scores, by row
            _scoreAffine_SmithWaterman(seq1[i - 1], seq2[j - 1], i, j, _dpTable, params);

            //update the max-score and its location, as needed
            Cell& cell = _dpTable[i][j];
            if (maxScore < cell.deletionScore) {
                maxScore = cell.deletionScore;
                maxIndices.first = i;
                maxIndices.second = j;
            }
            if (maxScore < cell.insertionScore) {
                maxScore = cell.insertionScore;
                maxIndices.first = i;
                maxIndices.second = j;
            }
            if (maxScore < cell.substitutionScore) {
                maxScore = cell.substitutionScore;
                maxIndices.first = i;
                maxIndices.second = j;
            }
        }
        //show progress
        _reportProgress(i, _dpTable.size());
    }

    //_printTable(_dpTable);

    cout << "\r\nBacktracking to find optimal alignment..." << endl;
    //global backtrack: from bottom-right cell to find optimal alignment, and track score
    int previousTraversal;
    i = maxIndices.first;
    j = maxIndices.second;
    alignment.Clear();
    alignment.maxScore = _maxThree(_dpTable[i][j].deletionScore, _dpTable[i][j].insertionScore, _dpTable[i][j].substitutionScore);
    //traverses the cells of the optimal path from cell(i,j) until no further positive path can be found
    while (i > 0 && j > 0 && (_maxThree(_dpTable[i][j].deletionScore,_dpTable[i][j].insertionScore,_dpTable[i][j].substitutionScore) > 0)) {
        Cell& cell = _dpTable[i][j];

        //deletion score is max, so traverse up
        if (cell.deletionScore >= cell.insertionScore && cell.deletionScore >= cell.substitutionScore) {
            alignment.bridge = " " + alignment.bridge;
            alignment.s1 = seq1[i - 1] + alignment.s1;
            alignment.s2 = "-" + alignment.s2;
            alignment.gaps++;
            //look in upper cell: was this an opening deletion-gap?
            if (previousTraversal != DEL) {
                alignment.openingGaps++;
            }
            previousTraversal = DEL;
            i--;
        }
        //insertion score is max, so traverse left
        else if (cell.insertionScore >= cell.deletionScore && cell.insertionScore >= cell.substitutionScore) {
            alignment.bridge = " " + alignment.bridge;
            alignment.s1 = "-" + alignment.s1;
            alignment.s2 = seq2[j - 1] + alignment.s2;
            alignment.gaps++;
            //was this an opening insertion-gap?
            if (previousTraversal != INS) {
                alignment.openingGaps++;
            }
            previousTraversal = INS;
            j--;
        }
        else { //assume substitution for all other cases, so traverse diagonally up and left
            alignment.s1 = seq1[i - 1] + alignment.s1;
            alignment.s2 = seq2[j - 1] + alignment.s2;
            if (seq1[i - 1] == seq2[j - 1]) {
                alignment.matches++;
                alignment.bridge = "|" + alignment.bridge;
            }
            else {
                alignment.mismatches++;
                alignment.bridge = " " + alignment.bridge;
            }
            previousTraversal = SUB;
            i--; j--;
        }
    }

    _clearTable();
}

void SequenceComparer::_printTable(const vector<vector<DpCell> >& dpTable)
{
    for (int i = 0; i < dpTable.size(); i++) {
        for (int j = 0; j < dpTable[i].size(); j++) {
            printf("%3d ",_maxThree(dpTable[i][j].deletionScore, dpTable[i][j].insertionScore, dpTable[i][j].substitutionScore));
        }
        cout << endl;
    }
    cout << endl;
}

void SequenceComparer::_scoreAffine_SmithWaterman(char a, char b, int row, int col, vector<vector<DpCell> >& dpTable, const Params& params)
{
    if (row <= 0 || col <= 0 || row >= dpTable.size() || col >= dpTable[0].size()) {
        cout << "BAIL; indices in minThree invalid: (row,col)=" << row << ":" << col << endl;
        exit(1);
    }

    //get the deletion score for this cell, per affine rules
    dpTable[row][col].deletionScore = _getAffineDeletionScore(dpTable[row - 1][col], params);
    //get the inserton score for this cell, per affine rules
    dpTable[row][col].insertionScore = _getAffineInsertionScore(dpTable[row][col - 1], params);
    //get the substitution score for this cell, per character match/mismatch
    dpTable[row][col].substitutionScore = _getAffineSubstitutionScore((a == b), dpTable[row - 1][col - 1], params);

    //reset any scores below zero (the smith waterman modification to needleman wunsch)
    if (dpTable[row][col].deletionScore < 0) {
        dpTable[row][col].deletionScore = 0;
    }
    if (dpTable[row][col].insertionScore < 0) {
        dpTable[row][col].insertionScore = 0;
    }
    if (dpTable[row][col].substitutionScore < 0) {
        dpTable[row][col].substitutionScore = 0;
    }
}

/*
Sets the deletion score for a given cell(i,j) according to affine rules.
For cell(i,j), the deletion score is:
    max {
        cell(i-1,j).deletionScore + g
        cell(i-1,j).insertionScore + h + g
        cell(i-1,j).substitutionScore + h + g
    }
*/
int SequenceComparer::_getAffineDeletionScore(const DpCell& predecessor, const Params& params)
{
    int del = predecessor.deletionScore + params.g;
    int sub = predecessor.substitutionScore + params.h + params.g;
    int ins = predecessor.insertionScore + params.h + params.g;
    
    return _maxThree(del,sub,ins);
}

int SequenceComparer::_getAffineSubstitutionScore(bool isMatch,const DpCell& predecessor, const Params& params)
{
    int max = _maxThree(predecessor.deletionScore,predecessor.substitutionScore,predecessor.insertionScore);

    if (isMatch) {
        max += params.match;
    }
    else {
        max += params.mismatch;
    }

    return max;
}

/*
Sets the insertion score for a given cell(i,j) according to affine rules.
For cell(i,j), the deletion score is:
max {
cell(i,j).insertionScore + g
cell(i,j - 1).deletionScore + h + g
cell(i,j - 1).substitutionScore + h + g
}
*/
int SequenceComparer::_getAffineInsertionScore(const DpCell& predecessor, const Params& params)
{
    int ins = predecessor.insertionScore + params.g;
    int sub = predecessor.substitutionScore + params.h + params.g;
    int del = predecessor.deletionScore + params.h + params.g;

    return _maxThree(del, sub, ins);
}

/*
    Scores a cell in the dp table according to left, right, upper predecessors, propagating these scores
    according to affine/Needleman-Wunsch rules. Scoring could be maintained here, but is performed in the
    backtracking step; this is the forward alg.
*/
void SequenceComparer::_scoreAffine_NeedlemanWunsch(char a, char b, int row, int col, vector<vector<DpCell> >& dpTable, const Params& params)
{
    if (row <= 0 || col <= 0 || row >= dpTable.size() || col >= dpTable[0].size()) {
        cout << "BAIL; indices in minThree invalid: (row,col)=" << row << ":" << col << endl;
        exit(1);
    }

    //get the deletion score for this cell, per affine rules
    dpTable[row][col].deletionScore = _getAffineDeletionScore(dpTable[row-1][col],params);
    //get the inserton score for this cell, per affine rules
    dpTable[row][col].insertionScore = _getAffineInsertionScore(dpTable[row][col-1],params);
    //get the substitution score for this cell, per character match/mismatch
    dpTable[row][col].substitutionScore = _getAffineSubstitutionScore( (a==b), dpTable[row-1][col-1], params);
}

int SequenceComparer::_maxThree(int a, int b, int c)
{
    int max = b;

    if (a > max) {
        max = a;
    }
    if (c > max) {
        max = c;
    }

    return max;
}

//returns a string filtered of all but characters in validChars
string filter(const string& s, const string& validChars)
{
    string filtered;

    for (int i = 0; i < s.length(); i++) {
        if (validChars.find( tolower(s[i]) ) != string::npos){
            filtered += s[i];
        }
    }

    return filtered;
}

/*

*/
void ParseFastaFile(const string fname, Sequence& seq1, Sequence& seq2)
{
    int seqnum;
    ifstream myReadFile;
    myReadFile.open(fname);
    string line;
    string validBases = "atcg";

    string& s1 = seq1.seq;
    string& s2 = seq2.seq;

    if (myReadFile.is_open()) {
        //init
        s1.clear();
        s2.clear();
        seqnum = 0;

        while (getline(myReadFile, line)) {
            //check if line begins with '>'; if so, advance state, and grab the description line
            if (line.length() > 0 && line[0] == '>') {
                if(seqnum==0) 
                    seq1.desc = line.substr(1,line.length()); 
                else
                    seq2.desc = line.substr(1,line.length()); 
                seqnum++;
            }
            else {
                switch (seqnum) {
                    //parsing sequence 1
                    case 1:
                        s1 += filter(line, validBases);
                    break;
                    //parsing sequence 2
                    case 2:
                        s2 += filter(line, validBases);
                    break;
                    default:
                        cout << "ERROR unknown seqnum: " << seqnum << endl;
                    break;
                }
            }
        }
    }
    myReadFile.close();

    //cout << "Parsed sequences:" << endl;;
    //cout << "s1 >>" << s1 << endl;
    //cout << "s2 >>" << s2 << endl;
}

void ParseParamsFile(const string& fname, Params& params)
{
    ifstream myReadFile;
    myReadFile.open(fname);
    string line, parameter, value;
    const char delim = ' ';
    int val;

    if (myReadFile.is_open()) {
        while (getline(myReadFile,line)){
            parameter = line.substr(0, line.find_first_of(delim));
            value = line.substr(line.rfind(delim), line.length());
            
            try {
                val = stoi(value);
            }
            catch (...){
                cout << "ERROR exception thrown for value: " << value << endl;
            }

            if (parameter == "match") {
                params.match = val;
            }
            else if (parameter == "mismatch") {
                params.mismatch = val;
            }
            else if (parameter == "g") {
                params.g = val;
            }
            else if (parameter == "h"){
                params.h = val;
            }
            else{
                cout << "ERROR param unmapped: " << parameter << endl;
            }
        }
    }
    myReadFile.close();

    cout << "Parsed params (mismatch,match,g,h): " << params.mismatch << "," << params.match << "," << params.g << "," << params.h << endl;
}

bool fileExists(const string& path)
{
  ifstream myStream(path);

  bool fileAccessible = !myStream.fail(); 
  if(fileAccessible){
    myStream.close();
  }

  return fileAccessible;
}



