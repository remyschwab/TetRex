#pragma once

#include <string>
#include <vector>
#include <stack>
#include "nfa_pointer.h"
#include "robin_hood.h"
#include "utils.h"
#include "index.h"


struct kState
{
    std::string qGram_;
    std::vector<kState *> outs_ = {};
    int marked_ = 0;
    bool start_ = 0;
};

struct keyState
{
    std::string qGramFrag_{};
    State *positionNFA_= nullptr;
    kState *home_ = nullptr;
};

struct Path
{
    uint16_t qPath_;
    kState* position_;
};

keyState* key(const std::string& qGramFrag, State* positionNFA_, kState* home);

Path* findPath(kState* position);

void oneStep(std::stack<keyState *>& stack, State* it_ptr, kState* kptr, std::string& qGram);

void firstPhase(State *it_ptr, std::vector<keyState *>& output, const size_t& q);

int linSearch(const std::vector<keyState *>& liste, keyState* obj);

void nextStep(std::stack<keyState *>& stack, keyState* input);

void nextKeys(std::vector<keyState *>& liste, keyState* input, kState* match);

std::vector<kState *> nfa2knfa(State* nfa_ptr, const int& q);

/*
 * Depth first search, generates the matrix with the possible paths
 */
template <typename Alphabet>
void dfs(kState* input, std::vector<std::vector<std::string>>& matrix, robin_hood::unordered_map<uint64_t, bitvector> &hash_to_bits, auto agent)
{
    std::vector<std::string> line{};
    std::stack<Path*> stack{};

    uint8_t qlength = input->qGram_.length();
    auto hash_adaptor = seqan3::views::kmer_hash(seqan3::ungapped{qlength});

    Path* p = findPath(input);
    stack.push(p);

    while(!stack.empty())
    {
        p = stack.top();

        if(p->position_->marked_ == 0)
        {
            auto acid_vec = convertStringToAlphabet<Alphabet>(p->position_->qGram_);
            auto digest = acid_vec | hash_adaptor;
            if(!hash_to_bits.count(digest[0]))
            {
                hash_to_bits[digest[0]] = agent.bulk_contains(digest[0]);
            }
            line.push_back(p->position_->qGram_);
            p->position_->marked_ = 1;
        }
        if(p->qPath_ < p->position_->outs_.size())
        {
            if(p->position_->outs_[p->qPath_]->qGram_ == "$")
            {
                matrix.push_back(line);
                p->qPath_++;
            }
            else
            {
                if(p->position_->outs_[p->qPath_]->marked_ == 0)
                {
                    stack.push(findPath(p->position_->outs_[p->qPath_]));
                }
                p->qPath_++;
            }
        }
        else
        {
            line.pop_back();
            p->position_->marked_ = 0;
            stack.pop();
            delete p;
        }
    }
}
