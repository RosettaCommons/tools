#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <filesystem>
#include <map>
#include <regex>
#include <cstdlib>

namespace fs = std::filesystem;

std::string exec(const char* cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}

std::string get_sequence_from_file(const std::string& filename) {
    std::ifstream file(filename);
    std::string line;
    std::getline(file, line); // Skip first line
    std::getline(file, line); // Get second line (sequence)
    return line;
}

std::vector<std::string> split(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

int main(int argc, char* argv[]) {
    if (argc < 3) { // Minimum arguments: program name, fasta file, MSA path, and at least one contig
        std::cerr << "Usage: " << argv[0] << " <fasta_file> <msa path> <contig1> [contig2] ..." << std::endl;
        return 1;
    }

    // Get sequence and sequence length
    std::string fa = argv[1];
    std::string seq = get_sequence_from_file(fa);
    int length = seq.length();

    // Get MSA directory path
    std::string msa_path = argv[2];

    // Get contigs
    std::vector<std::string> substrings;
    for (int i = 3; i < argc; ++i) {
        substrings.push_back(argv[i]);
    }

    // Print contigs for debugging
    std::cout << "Fasta file: " << fa << std::endl;
    std::cout << "MSA directory: " << msa_path << std::endl;
    std::cout << "Contigs:" << std::endl;
    for (const auto& substring : substrings) {
        std::cout << " - " << substring << std::endl;
    }

    // for each contig get start and end postitions
    std::map<std::string, int> start_positions;
    std::vector<int> starts, ends;
    for (const auto& substring : substrings) {
        std::vector<size_t> positions;
        size_t pos = seq.find(substring);

        while (pos != std::string::npos) {
            positions.push_back(pos);
            pos = seq.find(substring, pos + 1); // Look for the next occurrence
        }

        if (positions.empty()) {
            std::cout << "Warning: No positions found for substring '" << substring << "'" << std::endl;
        } else if (positions.size() > 1) {
            std::cout << "Warning: Multiple positions found for substring '" << substring << "'" << std::endl;
            std::cout << "Positions: ";
            for (const auto& p : positions) {
                std::cout << p << " ";
            }
            std::cout << std::endl;

            starts.push_back(positions[0]);
            start_positions[substring] = positions[0];
            ends.push_back(positions[0] + substring.length());
        } else {
            starts.push_back(positions[0]);
            start_positions[substring] = positions[0];
            ends.push_back(positions[0] + substring.length());
        }
    }

    // Sort start and end postitions
    std::sort(starts.begin(), starts.end());
    std::sort(ends.begin(), ends.end());

    // Add N- and C-terminal linkers
    std::string Nlinker_dash_string = starts[0] == 0 ? "" : std::string(starts[0], '-');
    std::string Clinker_dash_string = ends.back() == length ? "" : std::string(length - ends.back(), '-');

    // Get linker lengths
    std::vector<int> linker_lengths;
    for (size_t i = 1; i < starts.size(); ++i) {
        linker_lengths.push_back(starts[i] - ends[i-1]);
    }

    // Get linkers
    std::vector<std::string> linkers;
    for (int len : linker_lengths) {
        linkers.push_back(std::string(len, '-'));
    }

    // Print linkers
    std::cout << "Linkers:" << std::endl;
    for (const auto& linker : linkers) {
        std::cout << linker << std::endl;
    }

    // Sort the substrings by starting position
    std::vector<std::string> sorted_substrings;
    for (int pos : starts) {
        for (const auto& pair : start_positions) {
            if (pair.second == pos) {
                sorted_substrings.push_back(pair.first);
                break;
            }
        }
    }

    // Read target names into arrays
    std::vector<std::string> target_names;
    std::ifstream target_file(msa_path + "/target_names.txt");
    std::string line;
    while (std::getline(target_file, line)) {
        target_names.push_back(line);
    }

    // Read all sorted contigs into nested array
    std::vector<std::vector<std::string>> contigs;
    for (const auto& substring : sorted_substrings) {
        std::vector<std::string> temp;
        std::ifstream contig_file(msa_path + "/" + substring + ".txt");
        while (std::getline(contig_file, line)) {
            temp.push_back(line);
        }
        contigs.push_back(temp);
    }

    // Create output files
    std::string base_name = fa.substr(0, fa.find_last_of('.'));
    std::ofstream output_file(base_name + ".a3m");

    // Populate output file
    for (size_t i = 0; i < target_names.size(); ++i) {
        std::string line = Nlinker_dash_string;

        for (size_t j = 0; j < sorted_substrings.size(); ++j) {
            if (j < sorted_substrings.size() - 1) {
                line += contigs[j][i] + linkers[j];
            } else {
                line += contigs[j][i];
            }
        }

        line += Clinker_dash_string;

        int checkSum = line.length() - std::count_if(line.begin(), line.end(), [](char c) { return std::islower(c); });

        if (checkSum != length) {
            std::cout << "Checksum failed" << std::endl;
        }

        std::string target = split(target_names[i], '\t')[0];

        output_file << target << "\t" << target.substr(1) << std::endl;
        output_file << line << std::endl;
    }

    return 0;
}