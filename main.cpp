#include "Header.hpp"

/*
-56 	28/70(40%) 	0/70(0%) 	Plus/Plus

Test4 gives this, using parameters mismatch=-2, match=-1, h=-5 g=-2:
Query  1   GTGGCGCGAGCTTCTGAAACTAGGCGGCAGAGGCGGAGCCGCTGTGGCACTGCTGCGCCT  60
             |    | | | | |  | |  ||  | |  ||||   ||  | ||  | |    |||
Sbjct  1   CCGTGAGGCGTTGCCGTCAGTCAGCTACCGCTGCGGGAGCGGAGCGGGTCGGTGCGGCCG  60

Query  61  CTGCTGCGCC  70
               ||
Sbjct  61  GGTGTGGCGG  70


Test5.txt results (params same as above match=1,mismath=-2,g=-2,h=-5):
NW Score	Identities	Gaps	Strand
-252 	124/287(43%) 	14/287(4%) 	Plus/Plus

Query  1    GTGGCGCGAGCTTCTGAAACTAGGCGGCAGAGGCGGAGCCGCTGTGGCACTGCTGCGCCT  60
              |    | | | | |  | |  ||  | |  ||||   ||  | ||  | | |||| |
Sbjct  1    CCGTGAGGCGTTGCCGTCAGTCAGCTACCGCTGCGGGAGCGGAGCGGGTCGG-TGCGGC-  58

Query  61   CTGCTGCGCCTCGGGTGTCTTTTGCGGCGGTGGGTCGCCGCCGGGAGAAGCGTGAGGGGA  120
            | | || | |  | |||   |  | ||  ||| | |   |  ||||   || |||
Sbjct  59   CGGGTGTGGCGGGCGTGCGCTCCGGGGTCGTGAGGCCGTGA-GGGACGCGCCTGACACCC  117

Query  121  CAGATTTGTGACCGGCGCGGTTTTTGTCAGCTTACTCCGGCCAAAAAAGAACTGCACCTC  180
            | ||   | | | | ||   |  ||  | ||   |    ||           |||| | |
Sbjct  118  CGGAG--GAGCCAGTCGACCTCCTTCACGGCCACCCGGAGCAGCTGCCCCGGTGCAGCCC  175

Query  181  TGGAGCGGACTTATTTACCAAGCATTGGAGG--AATATCGTAGGTAAAAATGCCTAT-TG  237
              | |||  ||     |   | |    | |   |   ||   ||||      ||  | ||
Sbjct  176  GCG-GCGTCCTCCCAGAGGGATCCGGCGCGTCCAGAGTCCGCGGTAGCTGCCCCGTTCTG  234

Query  238  GA-TCCAAAGAGAGGCCAACATTTTTTGA---AATTTTTAAGACACG  280
             | ||    ||    | |||| ||  | |   |   |   ||| |
Sbjct  235  CAGTCGCCGGATTACCTAACACTTCCTTACCGAGCATCGGAGAAAT   280


Test1.txt:
NW Score	Identities	Gaps	Strand
-2 	6/10(60%) 	0/10(0%) 	Plus/Plus

Query  1  ATCGTTACGT  10
          ||| ||   |
Sbjct  1  ATCCTTGACT  10







*/

int main(int argc, char** argv)
{
    string fasta, config, seq1, seq2;
    bool isGlobal;
    Params params;
    Alignment alignment;

    /*
    if (argc != 4) {  // <executable name> <input file containing both s1 and s2> <0: global, 1: local> <optional: path to parameters config file>
        cout << "Usage: scmp [fname]" << endl;
        return 1;
    }

    fasta = argv[1];
     = argv[2][0] == '0' ? Alin
    config = argv[3];
    */

    //isGlobal = argv[2][0] == '0';

    fasta = "fasta.txt";
    //fasta = "test1.txt";
    //fasta = "test2.txt";
    //fasta = "test3.txt";
    //fasta = "test5.txt";
    //fasta = "ananthTest.txt";
    config = "parameters.config";
    //read the sequences
    ParseFastaFile(fasta, seq1, seq2);
    //parse the parameter object
    ParseParamsFile(config, params);

    SequenceComparer* stringComp = new SequenceComparer();

    stringComp->NeedlemanWunsch(seq1,seq2,params,alignment);
    stringComp->PrintResult(seq1, seq2, params, alignment);

    stringComp->SmithWaterman(seq1, seq2, params, alignment);
    stringComp->PrintResult(seq1, seq2, params, alignment);
    
    delete stringComp;

    return 0;
}


/*

SMithWatemrna, Test5.txt, run at
#=======================================
#
# Aligned_sequences: 2
# 1: EMBOSS_001
# 2: EMBOSS_001
# Matrix: EDNAFULL
# Gap_penalty: 5.0
# Extend_penalty: 1.0
#
# Length: 309
# Identity:     161/309 (52.1%)
# Similarity:   161/309 (52.1%)
# Gaps:         126/309 (40.8%)
# Score: 415.0
#
#
#=======================================

EMBOSS_001         1 GT--GGCG----CG--AGCTTCTGAAACTA--G--GCGGCAGAGGCGGAG     38
||  ||||    ||  ||  ||   |.|||  |  ||||  || ||||||
EMBOSS_001         3 GTGAGGCGTTGCCGTCAG--TC---AGCTACCGCTGCGG--GA-GCGGAG     44

EMBOSS_001        39 C----CGCTGTGGCAC---TGCT-GCGCCTCTGC-TGCGCCT-CGGGTGT     78
|    ||.||.||| |   || | |||    .|| |||| || |||| ||
EMBOSS_001        45 CGGGTCGGTGCGGC-CGGGTG-TGGCG----GGCGTGCG-CTCCGGG-GT     86

EMBOSS_001        79 CTTTTGCGGCGGT--GGGTCGCCGCC-G-------GGAGAAGC--GT-GA    115
|  .||.|||.||  |||.|| |||| |       ||||.|||  || ||
EMBOSS_001        87 C--GTGAGGCCGTGAGGGACG-CGCCTGACACCCCGGAGGAGCCAGTCGA    133

EMBOSS_001       116 ----------GG------GGA-CAGAT----TTGTG-A-CCGGCGCGGTT    142
||      ||| |||.|    ..||| | ||  |||||
EMBOSS_001       134 CCTCCTTCACGGCCACCCGGAGCAGCTGCCCCGGTGCAGCC--CGCGG--    179

EMBOSS_001       143 TTTGTCAGCTTAC-------TCCGG-----CCAAAAAAG--------AAC    172
.|||  ||..|       |||||     ||   |.||        |.|
EMBOSS_001       180 --CGTC--CTCCCAGAGGGATCCGGCGCGTCC---AGAGTCCGCGGTAGC    222

EMBOSS_001       173 TGCACC--TCTGGAG----CGGACTTA---------T--TTACCAAGCAT    205
|||.||  ||||.||    |||| |||         |  |||||.|||||
EMBOSS_001       223 TGCCCCGTTCTGCAGTCGCCGGA-TTACCTAACACTTCCTTACCGAGCAT    271

EMBOSS_001       206 TGGAGGAAT    214
.||||.|||
EMBOSS_001       272 CGGAGAAAT    280

Ananth test output:

Scores:    match = 1, mismatch = -2, h =-5, g = -2

Sequence 1 = "s1", length = 125 characters
Sequence 2 = "s2", length = 111 characters

s1  1    ACATGCTACACGTATCCGATACCCCGTAACCGATAACGATACACAGACCTCGTACGCTTG  60
|||||| ||||   ||||||||||||||||||||||||||||| ||||||||||||||||
s2  1    ACATGCGACACTACTCCGATACCCCGTAACCGATAACGATACAGAGACCTCGTACGCTTG  60

s1  61   CTACAACGTACTCTATAACCGAGAACGATTGACATGCCTCGTACACATGCTACACGTACT  120
|||           ||||||||||||||||||||| |||||||||   ||||||||||||
s2  61   CTA-----------ATAACCGAGAACGATTGACATTCCTCGTACA---GCTACACGTACT  106

s1  121  CCGAT  125
|||||
s2  107  CCGAT  111



Report:

Global optimal score = 55

Number of:  matches = 105, mismatches = 6, gaps = 14, opening gaps = 2

Identities = 105/125 (84%), Gaps = 14/125 (11%)


Fasta truncated results:
NW Score	Identities	Gaps	Strand
-949 	7608/11579(66%) 	1426/11579(12%) 	Plus/Plus

Ending with:
Query  11264  TCAACTAAAAATTCAAATACTTTAAATCAGAAGATTTCATAGTTAATTTATTTTTTTTTT  11323
||||
Sbjct  10338  -----------------------------GAAG---------------------------  10341

Query  11324  CAACAAAATGGTCATCCAAACTCAAACTTGAGAAAATATCTTGCTTTCAAATTGGCACT  11382
||
Sbjct  10342  ------------------------------------------------------GCGGG  10346





*/