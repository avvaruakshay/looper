/*
    Utils: Auxiliary functions for looper.cpp
    @file utils.h
    @author Akshay Avvaru
    @version 0.1 06/08/2020
*/
#include <cstdint>
#include <unordered_map>
#include <bitset>
#include <fstream>
#include <chrono>
#include <iomanip>

using namespace std;
using namespace chrono;


namespace utils {

    /*
     *  Check for length cutoff
     *  @param M Maximum motif size
     *  @param cutoff Cutoff length of repeat sequence
    */
    void length_cutoff_error(unsigned int M, unsigned int cutoff) {
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
    void motif_range_error(unsigned int m, unsigned int M) {
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
                        unsigned int &m, unsigned int &M, unsigned int &cutoff) {
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
    string bit2base(unsigned long long int seq, int l, int m) {
        string nuc = "";
        unsigned long long int fetch = 3ull << 2*(l-1);
        unsigned long long int c;
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
    unsigned long long int bit_reverse_complement(unsigned long long int seq, int l) {
        unsigned long long int rc = 0ull;
        unsigned long long int const NORM = ~(0ull) >> 2*(32-l);
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
    void count_seq(string filename, int &sequences, unsigned long long int &gsize, unsigned long long int &GC) {
        ifstream file(filename);
        string fline;
        while (getline(file, fline)) {
            if (fline[0] == '>') { sequences += 1; }
            else {
                for(const auto c: fline) {
                    gsize += 1;
                    switch(c) {
                        case 'c': case 'C': GC += 1; break;
                        case 'g': case 'G': GC += 1; break;
                        default: continue;
                    }
                }
            }
        }
        file.close();
    }

    /*
        Updating the progress bar
        @param start_time integer denoting start time of the program
        @param numseq number of sequences processed so far
        @param sequences total number of sequences in the fasta file
    */
    void update_progress_bar(unsigned long long int start_time, int numseq, int sequences) {
        const int BAR_WIDTH = 50;
        float progress = (((float) numseq) / ((float) sequences));
        unsigned long long int now = duration_cast<milliseconds>(
            system_clock::now().time_since_epoch()
        ).count();
        float total_time = now-start_time;
        int time_ps = int((total_time/float(numseq))*1000);
        float time_per_seq = float(time_ps)/1000.0;
        int total_secs = int(total_time);
        int h = total_secs/3600000;
        int m = (total_secs%3600000)/60000;
        int s = ((total_secs%3600000)%60000)/1000;
        int ms = ((total_secs%3600000)%60000)%1000;
        cout << setprecision(2) << h << ":" << m << ":" << s << ":" << ms << " ";
        cout << "[";
        int pos = BAR_WIDTH * progress;
        for (int i = 0; i < BAR_WIDTH; ++i) {
            if (i < pos) cout << "=";
            else if (i == pos) cout << ">";
            else cout << " ";
        }
        
        cout << "] " << "" << numseq << "/" << sequences << " seqs | ";
        cout << int(progress * 100.0) << "% | ";
        if (progress == 1) cout << time_per_seq << " sec/seq     " << endl;
        else cout << time_per_seq << " sec/seq" << "\r";

        cout.flush();
    }


    /*
     *  Calculates the repeat class of the sequence
     *  @param seq 64-bit integer representing 2-bit string of the sequence
     *  @param l length of the DNA sequence
     *  @param m length of the motif size
     *  @return string of repeat class motif with the strand orientation;
     *      example: "AGG+"
    */
    string get_repeat_class(unsigned long long int seq, int l, int m, unordered_map<string, string> &rClassMap) {
        string strand;
        // Throw error if length cutoff is smaller than 
        // twice the length of largest motif
        unsigned long long int expand = seq >> 2*(l-(2*m)); 
        unsigned long long int const NORM = ~(0ull) >> 2*(32-m);
        unsigned long long int min = ~(0ull);
        unsigned long long int cyc; unsigned long long int cyc_rc;
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
        string repeat_class = utils::bit2base(min, m, m);
        rClassMap[utils::bit2base(seq, l , m)] = repeat_class + strand;

        return repeat_class + strand;
    }

    /*
     *  Filters redundant motif sizes for division rule checks
     *  @param m minimum motif size
     *  @param M maximum motif size
     *  @return list of motif sizes to perform non-redundant checks
    */
    vector<unsigned int> get_motif_sizes(unsigned int m, unsigned int M) {
        vector<unsigned int> a = {M};
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
    static inline unsigned int check_atomicity(unsigned long long int seq, unsigned int l, unsigned int m) {
        seq = seq >> (2*(l-m));
        for (int i=1; i<m; i++) {
            if (m%i == 0) {
                unsigned long long int D = 0ull; unsigned int d = m/i;
                for (int j=0; j<d; j++) { D = D << 2*i; D += 1; }
                if (seq%D == 0) { return i; }
            }
        }
        return m;
    }
}