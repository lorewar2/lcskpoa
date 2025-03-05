#include <iostream>
#include <chrono>
#include <fstream> 
#include <iostream> 
#include <string> 
#include "spoa/spoa.hpp"

int main(int argc, char** argv) {
    // Check if the filename is provided as a command-line argument
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
        return 1; // Return an error code
    }

    std::string filename = argv[1]; // Get the filename from the command-line argument
    std::ifstream file(filename); // Open the file for reading

    // Check if the file was successfully opened
    if (!file.is_open()) {
        std::cerr << "Unable to open file: " << filename << std::endl;
        return 1; // Return an error code
    }
    std::vector<std::vector<std::string>> sequences_vec;
    std::vector<std::string> sequences;
    std::string line; // Temporary string to hold each line
    int line_index = 0;
    // Read the file line by line
    while (std::getline(file, line)) {
        if (line_index % 2 == 1) {
            sequences.push_back(line); // Add each line to the vector
            if (sequences.size() >= 3) {
                sequences_vec.push_back(sequences);
                sequences.clear();
            }
        }
        line_index++;
    }
    file.close(); // Close the file
    for(int i = 0; i < sequences_vec.size(); i++) {
        auto alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW, 2, -2, -2);  // linear gaps
        spoa::Graph graph{};

        int run_index = 0;
        std::chrono::steady_clock::time_point begin;
        for (const auto& it : sequences_vec[i]) {
            auto alignment = alignment_engine->Align(it, graph);
            graph.AddAlignment(alignment, it);
            if (run_index == 1) { //start time measurement
                begin = std::chrono::steady_clock::now();
            }
            run_index++;
        }

        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cout << i << "iter, time elasped = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
    }
    return 0;
}