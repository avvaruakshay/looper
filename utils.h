/*
    Utils: Auxiliary function for looper.cpp
    @file utils.h
    @author Akshay Avvaru
    @version 0.1 06/08/2020
*/
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <bitset>
#include <fstream>

using namespace std;


namespace utils {

    /*
     *  Check for length cutoff
     *  @param M Maximum motif size
     *  @param cutoff Cutoff length of repeat sequence
    */
    void length_cutoff_error(uint M, uint cutoff) {
        try { if (cutoff < 2*M) { throw 1; } }
        catch (int err) {
            cout << "Looper:" << endl;
            cout << endl << "\033[1m\033[31mLengthCutoffError: \033[0m"; 
            cout << "Length cutoff cannot be smaller than twice of ";
            cout << "maximum motif size" << '\n';
            exit (EXIT_FAILURE);
        }
    }
    
    /*
     *  Check for input file
     *  @param input Bool if file is good
     *  @param file_name Name of the input file
    */
    void input_file_error(bool input, string file_name) {
        try { if (!input) { throw 1; } }
        catch (int err) {
            cout << "Looper:" << endl;
            cout << "\033[1m\033[31mFileNotFoundError: \033[0m"; 
            cout << "File " << file_name << " doesn't exist \n";
            exit (EXIT_FAILURE);
        }
    }

    /*
     *  Check for input file
     *  @param input Bool if file is good
     *  @param file_name Name of the input file
    */
    void motif_range_error(uint m, uint M) {
        try { if (m > M) { throw 1; } }
        catch (int err) {
            cout << "Looper:" << endl;
            cout << endl << "\033[1m\033[31mMotifRangeError: \033[0m"; 
            cout << "Maximum motif size is smaller than minimum motif size.";
            exit (EXIT_FAILURE);
        }
    }

    
    /*
        Prints help message / usage of the program
    */
    void print_help() {
        cout << "usage: looper -i <file>"; 
        cout << " [-m <int>] [-M <int>] [-l <int>]";
        cout << " [-o <file>] " << endl << endl;

        cout << "Required arguments: " << endl;
        cout << "-i\t<file>\tInput fasta file" << endl << endl;
        cout << "Optional arguments: " << endl;
        cout << "-m\t<int>\tMinimum motif size. Default: 1" << endl;
        cout << "-M\t<int>\tMaximum motif size. Default: 6" << endl;
        cout << "-l\t<int>\tCutoff repeat length. Default: 2*M."<< endl;
        cout << " \t \tShould atleast be twice of maximum motif size." << endl;
        cout << "-o\t<file>\tOutput file name.";
        cout << "Default: Input file name + _looper.tsv"<< endl;
    }

    /*
        Parse command line arguments.
    */
    void parse_arguments(int argc, char* argv[], string &fin, string &fout,\
                        uint &m, uint &M, uint &cutoff) {
        for (int i=1; i < argc; ++i) {
            string arg = argv[i];
            if (arg == "-h") { utils::print_help(); exit (EXIT_SUCCESS);}
            else if (arg == "-i") {
                if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                    // Increment 'i' so we don't get the argument as the next argv[i].
                    fin = argv[i+1]; i++; 
                } else { // Uh-oh, there was no argument to the input file option.
                  cerr << "-i option requires one argument." << endl;
                } 
            } else if (arg == "-o") {
                if (i + 1 < argc) { fout = argv[i+1]; i++;  }
                else { cerr << "-o option requires one argument." << endl; } 
            } else if (arg == "-m") {
                if (i + 1 < argc) { m = atoi(argv[i+1]); i++;  }
                else { cerr << "-m option requires one argument." << endl; } 
            } else if (arg == "-M") {
                if (i + 1 < argc) { M = atoi(argv[i+1]); i++;  }
                else { cerr << "-M option requires one argument." << endl; } 
            } else if (arg == "-l") {
                if (i + 1 < argc) { cutoff = atoi(argv[i+1]); i++;  }
                else { cerr << "-l option requires one argument." << endl; } 
            }
        }
        if (m == 0) { m = 1; }
        if (M == 0) { M = 6; }
        if (fout == "") { fout = fin + "_looper.tsv"; }
        utils::motif_range_error(m, M);
        if (cutoff == 0) { cutoff = 2*M; }
    }


    /*
        Converts a 2-bit string to nucleotide sequence
        @param seq 64-bit integer representing 2-bit string of the sequence
        @param l length of the DNA sequence
        @return string of the nucleotide sequence
    */ 
    string bit2nuc(uint64_t seq, int l, int m) {
        string nuc = "";
        uint64_t fetch = 3ull << 2*(l-1);
        uint64_t c;
        int shift = 2*(l-1) ;
        for (int i=0; i<m; ++i) {
            c = (seq & fetch) >> shift;
            switch(c) {
                case 0: nuc+= "A"; break;
                case 1: nuc+= "C"; break;
                case 2: nuc+= "G"; break;
                case 3: nuc+= "T"; break;
                default: continue;
            }
            shift -= 2; fetch >>= 2;
        }
        return nuc;
    }


    /*
     *  Calculates the reverse complement of a DNA 2-bit string
     *  @param seq 64-bit integer representing 2-bit string of the sequence
     *  @param l length of the DNA sequence
     *  @return a 64-bit integer representing the reverse complement
    */
    uint64_t bit_reverse_complement(uint64_t seq, int l) {
        uint64_t rc = 0ull;
        uint64_t const NORM = ~(0ull) >> 2*(32-l);
        bitset<64> norm (NORM);
        seq = ~(seq); seq = seq & NORM;
        for (int i=0; i<l; i++) { 
            rc += (seq & 3ull) << 2*(l-1-i); seq = seq >> 2;
        }
        return rc;
    }


    /*
     *  Counts the number of sequences in the input fasta file
     *  @param filename name of the fasta file
     *  @return number of sequences in the file (int)
    */
    int count_seq(string filename) {
        ifstream file(filename);
        string fline;
        int seq_count = 0;
        while (getline(file, fline)) {
            if (fline[0] == '>') { seq_count += 1; }
        }
        file.close();
        return seq_count;
    }

    /*
     *  Calculates the repeat class of the sequence
     *  @param seq 64-bit integer representing 2-bit string of the sequence
     *  @param l length of the DNA sequence
     *  @param m length of the motif size
     *  @return string of repeat class motif with the strand orientation;
     *      example: "AGG+"
    */
    string get_repeat_class(uint64_t seq, int l, int m, unordered_map<string, string> &rClassMap) {
        string strand;
        // Throw error if length cutoff is smaller than 
        // twice the length of largest motif
        uint64_t expand = seq >> 2*(l-(2*m)); 
        uint64_t const NORM = ~(0ull) >> 2*(32-m);
        uint64_t min = ~(0ull);
        uint64_t cyc; uint64_t cyc_rc;
        int palindrome_check = 0;

        for (int i=0; i<m; i++) {
            cyc = expand & NORM;
            if (cyc < min) { min = cyc; strand = "+";}
            cyc_rc = utils::bit_reverse_complement(cyc, m);
            if (cyc_rc < min) { min = cyc_rc; strand = "-";}
            if (cyc == cyc_rc) { palindrome_check = 1;}
            expand = expand >> 2;
        }

        if (palindrome_check == 1) { strand = "+"; }
        string repeat_class = utils::bit2nuc(min, m, m);
        rClassMap[utils::bit2nuc(seq, l , m)] = repeat_class + strand;

        return repeat_class + strand;
    }

    /*
     *  Filters redundant motif sizes for division rule checks
     *  @param m minimum motif size
     *  @param M maximum motif size
     *  @return list of motif sizes to perform non-redundant checks
    */
    vector<uint> get_motif_sizes(uint m, uint M) {
        vector<uint> a = {M};
        int vsize = 0;
        for (int i=M-1; i >= m; --i) {
            bool check = false;
            for (int j=0; j < vsize; j++) {
                if (a[j] % i == 0) { check = true; break; }
            }
            if (!check) { a.push_back(i); vsize += 1;}
        }
        return a;
    }

    /*
     *  Calculates the atomicity of a motif
     *  @param seq 64-bit integer representing 2-bit string of the sequence
     *  @param l length of the DNA sequence
     *  @param m length of the motif size
     *  @return atomicity of the motif
    */
    uint check_atomicity(uint64_t seq, uint l, uint m) {
        seq = seq >> (2*(l-m));
        for (int i=1; i<m; i++) {
            if (m%i == 0) {
                uint64_t D = 0ull; uint d = m/i;
                for (int j=0; j<d; j++) { D = D << 2*i; D += 1; }
                if (seq%D == 0) { return i; }
            }
        }
        return m;
    }
}