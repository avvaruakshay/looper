#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "fastq_gzip_reader.h"

int main(int argc, char* argv[]) {
    fastq::Input a = fastq::Input();

    std::cout << argv[1] << std::endl;
    uint line_count = 1;
    fastq::Read curr_read = a.fetch();
    while (curr_read.valid) {
        std::cout << curr_read.identifier << "\n";
        std::cout << curr_read.sequence << "\n";
        std::cout << curr_read.separator << "\n";
        std::cout << curr_read.baseQual << "\n";
        curr_read = a.fetch();
        std::cout << curr_read.valid << "\n";
    }

    std::cout << a.currentCount() << "\n";
    return 0;
}

