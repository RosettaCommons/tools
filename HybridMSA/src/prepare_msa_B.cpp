#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <filesystem>
#include <iomanip>

namespace fs = std::filesystem;

std::vector<std::string> read_fasta(const std::string& filename) {
    std::vector<std::string> sequences;
    std::ifstream file(filename);
    std::string line;
    std::string current_seq;

    while (std::getline(file, line)) {
        if (line[0] == '>') {
            if (!current_seq.empty()) {
                sequences.push_back(current_seq);
                current_seq.clear();
            }
        } else {
            current_seq += line;
        }
    }
    if (!current_seq.empty()) {
        sequences.push_back(current_seq);
    }
    return sequences;
}

std::string repeat_char(char c, int n) {
    return std::string(n, c);
}

int main(int argc, char* argv[]) {
    if (argc < 4) { // Minimum arguments: program name, fasta file, MSA path, and cardinality
        std::cerr << "Usage: " << argv[0] << " <fasta_file> <a3m_file> <cardinality>" << std::endl;
        return 1;
    }

    // Store command-line arguments and get sequence length
    std::string fa = argv[1];
    std::string a3m = argv[2];
    std::string card = argv[3];
    std::vector<std::string> sequences = read_fasta(fa);
    int length = sequences[0].length();

    // Generate a3m files
    for (int i = 0; i < sequences.size(); ++i) {
        std::stringstream ss;
        ss << std::setw(2) << std::setfill('0') << i + 1;
        std::string i_formatted = ss.str();

        std::string output_filename = fa.substr(0, fa.length() - 3) + "_seq" + i_formatted + ".a3m";
        std::ofstream output_file(output_filename);

        std::string seq = sequences[i];

        if (card == "monomer") {
            output_file << "#" << length << "\t1" << std::endl;
        } else if (card == "homotrimer") {
            output_file << "#" << length << "\t3" << std::endl;
        }

        output_file << ">101" << std::endl;
        output_file << seq << std::endl;
        output_file << ">101" << std::endl;
        output_file << seq << std::endl;

        std::ifstream original_a3m(a3m);
        output_file << original_a3m.rdbuf();

        output_file.close();
    }

    return 0;
}
