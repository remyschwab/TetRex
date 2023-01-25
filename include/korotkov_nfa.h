#ifndef KOROTKOV_NFA_H
#define KOROTKOV_NFA_H

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
  ~keyState(){}
};

struct Path
{
  uint16_t qPath_;
  kState* position_;
  ~Path(){};
};

kState* kstate(const std::string& qGram);

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
template <typename T>
void dfs(kState* input, std::vector<std::vector<std::string>>& matrix, uint32_t &vector_idx,
          robin_hood::unordered_map<uint64_t, uint32_t> &hash_to_idx,
          std::vector<bitvector> &kmer_bitvex, T agent)
{
  std::vector<std::string> line{};
  std::stack<Path*> stack{};
  
  size_t qlength = input->qGram_.length();
  auto hash_adaptor = seqan3::views::kmer_hash(seqan3::ungapped{qlength});

  Path* p = findPath(input);
  stack.push(p);

  while(!stack.empty())
  {
    p = stack.top();

    if(p->position_->marked_ == 0)
    {
      auto acid_vec = convertStringToDNA(p->position_->qGram_);
      auto digest = acid_vec | hash_adaptor;
      hash_to_idx[digest[0]] = vector_idx;
      kmer_bitvex.push_back(agent.bulk_contains(digest[0]));
      vector_idx++;
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

#endif
