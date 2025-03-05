#include <iostream>
#include <chrono>
#include <fstream> 
#include <iostream> 
#include <string> 
#include <vector>
#include "include/abpoa.h"

unsigned char nt4_table[256] = {
    0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4 /*'-'*/, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

// 65,97=>A, 67,99=>C, 71,103=>G, 84,85,116,117=>T, else=>N
const char nt256_table[256] = {
    'A', 'C', 'G', 'T',  'N', '-', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', '-',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'A', 'N', 'C',  'N', 'N', 'N', 'G',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'T', 'T', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'A', 'N', 'C',  'N', 'N', 'N', 'G',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'T', 'T', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N'
};

int main(int argc, char** argv) {
    std::chrono::steady_clock::time_point begin;
    
    int i, j, n_seqs = 3;
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
    for(int index = 0; index < sequences_vec.size(); index++) {
        char seqs[3][90000] = {'0'}; // Initialize array with null characters
        for (size_t j = 0; j < sequences_vec[index].size() && j < 3; ++j) {
            strncpy(seqs[j], sequences_vec[index][j].c_str(), 90000 - 1);
            seqs[j][90000 - 1] = '\0'; // Ensure null termination
        }
        // initialize variables
        
        begin = std::chrono::steady_clock::now();
        // initialize variables
        abpoa_t *ab = abpoa_init();
        abpoa_para_t *abpt = abpoa_init_para();

        // alignment parameters
        abpt->align_mode = 0; // 0:global 1:local, 2:extension
        // abpt->mat_fn = strdup("HOXD70.mtx"); abpt->use_score_matrix = 1; // score matrix instead of constant match/mismatch score
        abpt->match = 2;      // match score
        abpt->mismatch = -2;   // mismatch penalty
        // abpt->gap_mode = ABPOA_CONVEX_GAP; // gap penalty mode
        abpt->gap_open1 = 2;  // gap open penalty #1
        abpt->gap_ext1 = 2;   // gap extension penalty #1
        abpt->gap_open2 = 2; // gap open penalty #2
        abpt->gap_ext2 = 2;   // gap extension penalty #2
                                // gap_penalty = min{gap_open1 + gap_len * gap_ext1, gap_open2 + gap_len * gap_ext2}
        // abpt->bw = 10;        // extra band used in adaptive banded DP
        // abpt->bf = 0.01; 
        
        // output options
        abpt->out_msa = 0; // generate Row-Column multiple sequence alignment(RC-MSA), set 0 to disable
        abpt->out_cons = 1; // generate consensus sequence, set 0 to disable
        abpt->w = 6, abpt->k = 9; abpt->min_w = 10; // minimizer-based seeding and partition
        abpt->progressive_poa = 0;
        abpt->max_n_cons = 0; // to generate 2 consensus sequences
        // abpt->sub_aln = 1;

        abpoa_post_set_para(abpt);

        // collect sequence length, trasform ACGT to 0123
        int *seq_lens = (int*)malloc(sizeof(int) * n_seqs);
        uint8_t **bseqs = (uint8_t**)malloc(sizeof(uint8_t*) * n_seqs);
        int **weights = (int**)malloc(sizeof(int*) * n_seqs);
        for (i = 0; i < n_seqs; ++i) {
            seq_lens[i] = strlen(seqs[i]);
            bseqs[i] = (uint8_t*)malloc(sizeof(uint8_t) * seq_lens[i]);
            weights[i] = (int*)malloc(sizeof(int) * seq_lens[i]);
            for (j = 0; j < seq_lens[i]; ++j) {
                bseqs[i][j] = nt4_table[(int)seqs[i][j]];
                if (j >= 12) weights[i][j] = 2;
                else weights[i][j] = 0;
            }
        }
        abpt->use_qv = 1;
        abpoa_msa(ab, abpt, n_seqs, NULL, seq_lens, bseqs, weights, NULL);

        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::cout << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << std::endl;
        for (i = 0; i < n_seqs; ++i) { free(bseqs[i]); free(weights[i]); }
        free(bseqs); free(seq_lens); free(weights);

        // free abpoa-related variables
        abpoa_free(ab); abpoa_free_para(abpt); 
    }
    return 0;
}