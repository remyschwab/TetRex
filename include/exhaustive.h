#ifndef UNTITLED_EXHAUSTIVE_H
#define UNTITLED_EXHAUSTIVE_H
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include "heuristic1.h"

// HELPER FUNCTIONS
// ---------------------------------------------------------------------------------------
std::map <char, int> create_alphabet_dict(std::string& alphabet){
    // This is a temporary function required for the funciton below
    // that turns an alphabet {'A', 'C', 'G', 'T'}; into a dictionary {'A': 0, 'C': 1, 'G': 2, 'T': 3};
    std::map <char, int>  alphabet_dict;
    for (int i = 0; i < alphabet.size(); i++) {
    alphabet_dict[alphabet[i]] = i;
    }
    return alphabet_dict;
}

std::map <char, int> alphabet_dict;
int nucl_to_rank(char& nucl){
    // This is a temporary function to obtain the rank of a letter in the alphabet.
    // This should later be replaced when using SeqAn alphabet.
    return alphabet_dict[nucl];
}


std::vector<std::vector<float>>  matrix_operation(std::vector<std::vector<float>> & m, std::size_t &n, std::size_t &s,
                                                    float (*func)(float)){
    // A small helper function that can be used to do operations, such as taking logarithms, on matrices of s x n in size.
    for(int i=0; i<s; i++) {
        for(int j=0; j<n; j++){
            m[i][j] = func(m[i][j]);
        }
    }
    return m;
}

// ---------------------------------------------------------------------------------------
void enumerate_strings(std::vector<std::vector<float>>  & profile, std::size_t &n, std::size_t &s, std::string alphabet, float T,
                        std::vector<std::string>& strings, std::vector<float>& scores, float& best_score,
                        int ik, std::vector<float>&score_cashe, std::string& string){
    // This enumerates the strings on demand. It starts at the first position and iterates recursively, depth-first through
    // all the letters in the alphabet. E.g. for n=3, it first creates AAA. If this string reaches the threshold,
    // it will be added to the strings vector.
    for (int a = 0; a < s; a++) {
        float score = score_cashe[ik] + profile[ik][a];
        if (score > T) {
            string[ik] = alphabet[a];
            if (ik == n - 1) {
                strings.push_back(string);
                scores.push_back(score+  best_score);
            } else {
                score_cashe[ik + 1] = score;
                enumerate_strings(profile, n,s, alphabet, T, strings, scores, best_score, (ik + 1), score_cashe, string);
            }
        }
    }
}

auto exhaustive(std::vector<std::vector<float>>  profile, std::size_t &n, std::size_t &s, std::string alphabet, float T){
    // This function does the exhaustive search, and is mainly a wrapper for the recursive function that enumerates
    // the strings.
    std::string best_string(n, ' '); float best_score=0;
    get_best_string (profile,n,s, best_string, best_score, alphabet, 1);
    T -= best_score;

    std::vector<std::string> strings;
    std::vector<float> scores={};
    std::vector<float> score_cashe(n,0); //init with n elements.
    std::string string = std::string(n, ' ');
    enumerate_strings(profile, n,s, alphabet,  T, strings, scores, best_score, 0, score_cashe, string);
    return std::make_tuple(std::move(strings), std::move(scores));
}







#endif //UNTITLED_EXHAUSTIVE_H
