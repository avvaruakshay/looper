#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

int main(int argc, char* argv[]) {
    
    // char *inp = argv[1];
    std::cout << argv[1] << std::endl;
    std::cout << argv[2] << std::endl;
    // char c;
    // std::cin >> c;
    // std::cout << c << std::endl;
    // char b;
    // while (std::cin >> b) {        
    //     std::cout << b << std::endl;
    // }

    // char c;
    // std::vector<int> cl;
    std::string line;
    uint line_count = 1;
    // std::getline(std::cin, line);
    // std::istringstream iss(line);
    for (int i=0; i< 5; i++) {
        while ( line_count % 5 != 0 && std::getline(std::cin, line) ) {
            std::cout << line << " " << line_count << "\n";
            line_count++;
        }
        line_count = 1;
    }
    return 0;
}

