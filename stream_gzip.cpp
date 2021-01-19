#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "fastq_gzip_reader.h"

int main(int argc, char* argv[]) {
    fastq::Input a = fastq::Input();

    // char *inp = argv[1];
    std::cout << argv[1] << std::endl;
    if (argv[2]) {
        std::cout << argv[2] << std::endl;
    }

    fastq::Read read = a.fetch();
    std::cout << read.identifier << "\n";
    std::cout << read.sequence << "\n";
    std::cout << read.separator << "\n";
    std::cout << read.baseQual << "\n";
    return 0;
}

