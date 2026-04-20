#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <filesystem>
#include <regex>
#include <cstdlib>
#include <sstream>

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
    if (argc < 4) { // Minimum arguments: program name, MSA path, and at least one substring
        std::cerr << "Usage: " << argv[0] << " <path_to_reformat.pl> <path_to_msas_directory> <contig1> [contig2] ..." << std::endl;
        return 1;
    }

    // Get MSA directory path
    std::string reformat = argv[1];

    // Get MSA directory path
    std::string path = argv[2];

    // Get contigs
    std::vector<std::string> substrings;
    for (int i = 3; i < argc; ++i) {
        substrings.push_back(argv[i]);
    }

    // Print contigs
    std::cout << "MSA directory: " << path << std::endl;
    std::cout << "Contigs:" << std::endl;
    for (const auto& substring : substrings) {
        std::cout << " - " << substring << std::endl;
    }

    // Reformat bfd_uniclust_hits.a3m
    std::string cmd = reformat + " a3m a3m " + path + "/bfd_uniclust_hits.a3m bfd_uniclust.a3m";
    system(cmd.c_str());

    // Reformat bfd_uniclust.a3m
    cmd = "awk '/^>/ {if (seq) print seq; printf \"%s\\n\", $0; seq = \"\"; next} {gsub(/[ \\t\\r\\n]/,\"\"); seq = seq $0} END {if (seq) print seq}' bfd_uniclust.a3m > bfd_uniclust_reformatted.a3m";
    system(cmd.c_str());

    // Concatenate reformatted files
    system("cat bfd_uniclust_reformatted.a3m > msa.a3m");

    std::string msa = "msa.a3m";
    std::string seq;

    // Read the second line of msa.a3m
    std::ifstream msa_file(msa);
    std::string line;
    std::getline(msa_file, line); // Skip first line
    std::getline(msa_file, seq);

    // Calculate start and end positions
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
            ends.push_back(positions[0] + substring.length() - 1);
        } else {
            starts.push_back(positions[0]);
            ends.push_back(positions[0] + substring.length() - 1);
        }
    }

    // Print starts and ends
    for (int start : starts) std::cout << start << " ";
    std::cout << std::endl;
    for (int end : ends) std::cout << end << " ";
    std::cout << std::endl;

    // Create output files
    std::ifstream input_msa(msa);
    std::ofstream target_names("target_names.txt");
    std::vector<std::ofstream> substring_files;
    for (const auto& substring : substrings) {
        substring_files.emplace_back(substring + ".txt");
    }

    // Populate output files
    int index = 0;
    std::getline(input_msa, line);
    std::getline(input_msa, line);
    while (std::getline(input_msa, line)) {
        if (index % 2 == 0) {
            target_names << line << std::endl;
        } else {
            int counter = 0;
            int saved_position = 0;
            int found = 0;

            for (int i = 0; i < line.length(); ++i) {
                char c = line[i];
                if (std::find(starts.begin(), starts.end(), counter) != starts.end() && 
                    (c == '-' || std::isupper(c))) {
                    saved_position = i;
                    int checkSum = i - std::count_if(line.begin(), line.begin() + i, 
                                                     [](char ch) { return std::islower(ch); });
                    if (starts[found] != checkSum) {
                        std::cout << "Element " << starts[found] << " is not equal to " << checkSum << std::endl;
                    }
                } else if (std::find(ends.begin(), ends.end(), counter) != ends.end() && 
                           (c == '-' || std::isupper(c))) {
                    substring_files[found] << line.substr(saved_position, i - saved_position + 1) << std::endl;
                    int checkSum = i - saved_position + 1 - 
                                   std::count_if(line.begin() + saved_position, line.begin() + i + 1, 
                                                 [](char ch) { return std::islower(ch); });
                    if (substrings[found].length() != checkSum) {
                        std::cout << "Element " << substrings[found].length() << " is not equal to " << checkSum << std::endl;
                    }
                    ++found;
                }
                if (c == '-' || std::isupper(c)) {
                    ++counter;
                }
            }
        }
        ++index;
    }

    // Remove lines where all contigs are blank
    std::vector<int> lines_to_remove;
    int line_count = 0;
    std::ifstream count_file("target_names.txt");
    while (std::getline(count_file, line)) ++line_count;

    for (int i = 1; i <= line_count; ++i) {
        bool non_dash = false;
        for (const auto& substring : substrings) {
            std::ifstream file(substring + ".txt");
            std::string line;
            for (int j = 1; j < i; ++j) std::getline(file, line);
            std::getline(file, line);
            if (std::regex_search(line, std::regex("[A-WYZ]"))) {
                non_dash = true;
                break;
            }
        }
        if (!non_dash) {
            lines_to_remove.push_back(i);
            std::cout << "rm " << i << std::endl;
        }
    }

    if (!lines_to_remove.empty()) {
        std::string sed_command = "sed -i '";
        for (int line : lines_to_remove) {
            sed_command += std::to_string(line) + "d;";
        }
        sed_command.pop_back();  
        sed_command += "' target_names.txt";
        system(sed_command.c_str());

        for (const auto& substring : substrings) {
            sed_command = "sed -i '";
            for (int line : lines_to_remove) {
                sed_command += std::to_string(line) + "d;";
            }
            sed_command.pop_back();
            sed_command += "' " + substring + ".txt";
            system(sed_command.c_str());
        }
    }

    return 0;
}
