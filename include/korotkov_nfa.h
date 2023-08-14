#pragma once

#include <string>
#include <vector>
#include <stack>
#include "nfa_pointer.h"
#include "robin_hood.h"
#include "utils.h"
#include "index.h"
#include "robin_hood.h"


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

std::vector<kState *> nfa2knfa(State* nfa_ptr, const int& q);

Path* findPath(kState* position);

static void dfs_old(kState* input, std::vector<std::vector<std::string>>& matrix)
{
  std::vector<std::string> line{};
  std::stack<Path*> stack{};

  Path* p = findPath(input);
  stack.push(p);

  while(!stack.empty())
  {
    p = stack.top();

    if(p->position_->marked_ == 0)
    {
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

/*
 * Depth first search, generates the matrix with the possible paths
 */
void dfs(kState* input, std::vector<std::vector<uint64_t>>& matrix, 
            robin_hood::unordered_map<uint64_t, bitvector> &hash_to_bits,
            auto &agent, IndexStructure &ibf)
{
    std::vector<uint64_t> line{};
    std::stack<Path*> stack{};

    uint8_t qlength = input->qGram_.length();
    std::vector<unsigned char> digit_vec(qlength);
    uint64_t digest;
    uint64_t rev_digest; // Only used for nucleotides
    Path* p = findPath(input);
    stack.push(p);

    while(!stack.empty())
    {
        p = stack.top();
        if(p->position_->marked_ == 0)
        {
            if(ibf.molecule_ == "aa")
            {
                convertStringToVec(p->position_->qGram_, ibf, digit_vec);
                decompose_peptide_record(digit_vec, 0, ibf);
                digest = ibf.forward_store_;
            }
            else
            {
                digest = encode_dna(p->position_->qGram_);
                rev_digest = revComplement(digest, ibf.k_);
                digest = digest <= rev_digest ? digest : rev_digest;
            }
            if(!hash_to_bits.count(digest)) hash_to_bits[digest] = agent.bulk_contains(digest);
            line.push_back(digest);
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
