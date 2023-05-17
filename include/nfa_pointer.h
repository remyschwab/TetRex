#pragma once

#include <string>
#include <vector>
#include "index.h"


enum
{
    Match = 256,
    Split = 257
};


/*
 * Represents an NFA state plus zero or one or two arrows exiting.
 * if c == Match, no arrows out; matching state.
 * If c == Split, unlabeled arrows to out and out1 (if != NULL).
 * If c < 256, labeled arrow with character c to out.
*/

struct State
{
    int c_;
    State *out1_ = nullptr;
    State *out2_ = nullptr;
    int lastlist_ = 0;
};


/*
 * A partially built NFA without the matching state filled in.
 * Frag.start points at the start state.
 * Frag.out is a list of places that need to be set to the
 * next state for this fragment.
*/
struct Frag
{
    State *start;
    std::vector<State *> out;
};

std::vector<State *> appendVec(const std::vector<State *>& vec1, const std::vector<State *>& vec2);

void patchVec(std::vector<State *>& in, State *s);

std::vector<State *> getVec(State *input);

std::string getRandomWord(State* startptr);

void add(std::vector<State* >& a, State* b);

void deleteGraph(State* startptr);

State* post2nfaE(const std::string& postfix);

void deleteGraph(State* startptr);
