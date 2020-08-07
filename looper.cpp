#include <fstream>
#include <iostream>
#include <math.h>
#include <string.h>
#include <unordered_map>
#include <vector>
#include <chrono>
#include <assert.h>

#include "utils.h"

using namespace std;
using namespace std::chrono;

// global variable
unordered_map<string, string> rClassMap;

// Data structure to track 2-bit sequence of scanned window in the genome
// and it's location
struct bitSeqWindow {
    uint64_t seq, count;
    signed int cutoff;
    bitSeqWindow() { reset(); }
    void reset() { seq = count = cutoff = 0;}
};


/* Main function of LOOPER */
int main(int argc, char* argv[]) {
    ios_base::sync_with_stdio(false);

    string fin, fout;
    uint m = 0, M = 0, cutoff = 0;
    if (argc == 1) { utils::print_help(); exit (EXIT_FAILURE);}
    else if (argc > 1) { 
        utils::parse_arguments(argc, argv, fin, fout, m, M, cutoff);
        utils::length_cutoff_error(M, cutoff);
    }


    int sequences = utils::count_seq(fin); // total number of sequences
    ifstream ins(fin); // input fasta file
    utils::input_file_error(ins.good(), fin);
    ofstream out(fout); // output file
    string line;
    bitSeqWindow window;    

    cout << endl << "Searching for tandem repeats in " << argv[1] << endl;
    cout << "Min-motif: " << m << "\t Max-motif: " << M;
    cout << "\t Length-cutoff: " << cutoff <<  endl << endl;

    cout

    uint64_t start_time = duration_cast<milliseconds>(
        system_clock::now().time_since_epoch()
    ).count();

    // integer tracking the start of the repeat
    // -1 indicates no repeat is found 
    int start = -1; 
    uint atomicity; // atomicity of the repeat
    uint end;       // end of the repeat    
    uint rlen;      // length of the repeat
    int repeat_check;  // bool tracking if the window sequence is a repeat
    int numseq = 0; // current sequence number
    int const BAR_WIDTH = 60; // progress of the progress bar
    string motif, repeat_class, seq_name;

    // NORM is require to fetch the current window sequence
    uint64_t const NORM = ~(0ull) >> 2*(32-cutoff);
    // non-redundant list of motifs used for checks
    vector<uint> motif_checks = utils::get_motif_sizes(m, M);
    const uint N = motif_checks.size();
    uint64_t divisor[N]; // list of divisors
    uint rem_shift[N]; // list of remainder sizes
    for (int i=0; i<N; i++) {
        uint d = cutoff / motif_checks[i];
        uint r = cutoff % motif_checks[i];
        uint64_t D = 0ull;
        for (int j=0; j<d; j++) { D = D << (2*motif_checks[i]); D += 1; }
        D = D << (2*r);
        divisor[i] = D;
        rem_shift[i] = 2*(cutoff - r);
    }

    while(getline(ins, line)) {
        if (line[0] == '>') {
            float progress = ((float) numseq) / ((float) sequences);
            if (start != -1) {
                end = window.count; rlen = end - start;
                out << seq_name << "\t" << start << "\t" << end \
                << "\t" << repeat_class.substr(0, atomicity) << "\t" << rlen \
                << "\t" << repeat_class.substr(atomicity, 1) << "\t" \
                << rlen/atomicity << "\t" << motif << endl;
            }
            seq_name = line.substr(1);
            window.reset(); start = -1;
            
            uint64_t now = duration_cast<milliseconds>(
                system_clock::now().time_since_epoch()
            ).count();
            float total_time = float(now-start_time)/1000.0;
            float time_per_seq = total_time/float(numseq);

            cout << "Time elapsed: " << total_time << " secs\n";
            cout << "[";
            int pos = BAR_WIDTH * progress;
            for (int i = 0; i < BAR_WIDTH; ++i) {
                if (i < pos) cout << "=";
                else if (i == pos) cout << ">";
                else cout << " ";
            }
            
            cout << "] " << "" << numseq << "/" << sequences << " seqs | ";
            cout << int(progress * 100.0) << "% | ";
            cout << time_per_seq << " sec/seq" << "\x1b[A\r";

            cout.flush();
            numseq++;
        }
        else {
            for(const auto c: line) {
                switch(c) {
                    case 'a': case 'A': break;
                    case 'c': case 'C': window.seq |= 1ull; break;
                    case 'g': case 'G': window.seq |= 2ull; break;
                    case 't': case 'T': window.seq |= 3ull; break;
                    case 'N': case 'n': 
                        window.seq = 0; window.cutoff = -1;
                        if (start != -1) {
                            end = window.count; rlen = end - start;
                            out << seq_name << "\t" << start << "\t" << end \
                            << "\t" << repeat_class.substr(0, atomicity) << \
                            "\t" << rlen << "\t" << \
                            repeat_class.substr(atomicity, 1) << "\t" \
                            << rlen/atomicity << "\t" << motif << endl;
                        }
                        start = -1;
                        break;
                    default: continue;
                }
                window.count += 1;
                window.cutoff += 1;
                window.seq &= NORM;

                if (window.cutoff >= cutoff) { // To be optimized
                    repeat_check = 0;
                    for (int i=0; i<N; ++i){
                        if ( (window.seq % divisor[i]) == ( window.seq >> rem_shift[i]) ) {
                            if (start == -1) {
                                atomicity = utils::check_atomicity(window.seq, cutoff, motif_checks[i]);
                                
                                // atomicity should be greater than 
                                // minimum motif-size
                                if (atomicity >= m) {
                                    start = window.count - cutoff;
                                    motif = utils::bit2nuc(window.seq, cutoff, atomicity);
                                    if (rClassMap.find(motif) != rClassMap.end()) { 
                                        repeat_class = rClassMap[motif];
                                    } else {
                                        repeat_class = utils::get_repeat_class(window.seq, cutoff, atomicity);
                                    }
                                }
                            }
                            repeat_check = 1; break;
                        }
                    }
                    if (repeat_check == 0 & start != -1) {
                        end = window.count - 1; rlen = end - start;
                        out << seq_name << "\t" << start << "\t" << end \
                        << "\t" << repeat_class.substr(0, atomicity) << "\t" \
                        << rlen << "\t" << repeat_class.substr(atomicity, 1) \
                        << "\t" << rlen/atomicity << "\t" << motif << endl;
                        start = -1;
                    }
                }
                window.seq <<= 2;
            }
        }
    }
    if (start != -1) {
        end = window.count; rlen = end - start;
        out << seq_name << "\t" << start << "\t" << end \
        << "\t" << repeat_class.substr(0, atomicity) << "\t" << rlen \
        << "\t" << repeat_class.substr(atomicity, 1) << "\t" \
        << rlen/atomicity << "\t" << motif << endl;
    }

    uint64_t end_time = duration_cast<milliseconds>(
        system_clock::now().time_since_epoch()
    ).count();
    float total_time = float(end_time - start_time)/1000.0;
    float time_per_seq = (float(end_time-start_time)/1000.0)/float(numseq);
    
    cout << "Time elapsed: " << total_time << " secs" << endl;
    cout << "[";
    for (int i = 0; i < BAR_WIDTH; ++i) {
        if (i < BAR_WIDTH) cout << "=";
        else if (i == BAR_WIDTH) cout << ">";
        else cout << " ";
    }
    cout << "] " << numseq << "/" << sequences << " seqs | " << "100% | ";
    cout << time_per_seq << " sec/seq" << endl;
    cout << "Total time taken: " << total_time << " seconds" << endl;

    ins.close(); out.close();
    return EXIT_SUCCESS;
}