#include <algorithm>
#include <iostream>
#include <math.h>
#include <queue>
#include <set>
#include <stack>
#include <string>
#include <vector>
#include "nfa_pointer.h"
#include "korotkov_nfa.h"
#include "graphMaker.h"



kState* kstate(const std::string& qGram)
{
  kState* ptr = new kState;

  ptr-> qGram_ = qGram;

  return ptr;
}

keyState* key(const std::string& qGramFrag, State *positionNFA, kState * home)
{
  keyState* ptr = new keyState;

  ptr-> qGramFrag_ = qGramFrag;
  ptr-> positionNFA_ = positionNFA;
  ptr-> home_ = home;

  return ptr;
}

/*
 * Helpfunction of firstPhase
 * represents one Step in the Automaton 
 */

void oneStep(std::stack<keyState *>& stack, State* itptr, kState* kptr, std::string& qGram)
{
  keyState *e1, *e2;
  int c = itptr->c_;
  switch(c)
  {
    default:
        qGram += c;
        e1 = key(qGram, itptr->out1_, kptr);
        stack.push(e1);
        break;
    case Split:
        e1 = key(qGram, itptr->out1_, kptr);
        e2 = key(qGram, itptr->out2_, kptr);
        stack.push(e2);
        stack.push(e1);
        break;
    case Match:
        throw int();
        break;
  }
}

/*
 * represents the first phase in the construction of the Automaton
 * creates the keys that are accessible from the start node 
 * Throws an error, if the q-gram length longer than the shortest possible q-gram
 */

void firstPhase(State* it_ptr, std::vector<keyState *>& output, const size_t& q)
{
  std::stack<keyState *> stack;
  keyState* k;
  std::string qGram = "";

  oneStep(stack,it_ptr,nullptr,qGram);

  while(!stack.empty())
  {
    k = stack.top();
    stack.pop();
    if(k->positionNFA_->c_ == Match)
    {
      delete k;
      throw int();
      //std::cerr<<"Q-Gram zu lang gewählt"<<"\n";
    }
    if(k-> qGramFrag_.size() == q-1 && k->positionNFA_->c_ != Split)
    {
      //std::cout<<k->qGramFrag_<<"\n";
      output.push_back(k);
    }
    else
    {
      oneStep(stack, k->positionNFA_, nullptr, k->qGramFrag_);
      delete k;
    }
  }
}

/*
 * linSearch and linSearchK are simple linear search for kStates and key objects
 * can be optimized in example with an hash table
 */

int linSearchK(const std::vector<kState *>& liste, std::string obj)
{
  for(size_t i = 0; i<liste.size(); i++)
  {
    if(liste[i]->qGram_ == obj)
    {
      return i;
    }
  }
  return -1;
}

int linSearch(const std::vector<keyState *>& liste, keyState* obj)
{
  for(size_t i = 0; i<liste.size(); i++)
  {
    if(liste[i]->qGramFrag_ == obj->qGramFrag_ && liste[i]-> positionNFA_ == obj->positionNFA_)
    {
      return i;
    }
  }
  return -1;
}


/*
 * Helpfunction for nextKeys and is part of the sec. Phase
 * is similar to oneStap 
 * Match State is not accesable
 */

void nextStep(std::stack<keyState *>& stack, keyState* input)
{
  keyState *e1, *e2;
  State* itptr = input->positionNFA_;
  int c = itptr->c_;
  switch(c)
  {
    default:
        stack.push(input);
        break;
    case Split:
        e1 = key(input->qGramFrag_, itptr->out1_, nullptr);
        e2 = key(input->qGramFrag_, itptr->out2_, nullptr);
        stack.push(e2);
        stack.push(e1);
        break;
  }
}

/*
 * creates the other keys and nodes for the korotkov Automaton
 * get the vec of key from firstPhase
 */

void nextKeys(std::vector<keyState *>& liste, keyState* input, kState* match)
{
  std::stack<keyState *> stack;
  keyState* k;
  kState* e;
  std::string qGramFrag = input->home_->qGram_;
  qGramFrag = qGramFrag.substr(1);
  std::string qGram = qGramFrag;

  k = key(qGramFrag, input->positionNFA_->out1_, nullptr);
  nextStep(stack, k);

  while(!stack.empty())
  {
    k = stack.top();
    stack.pop();
    if(k->positionNFA_->c_== Match)
    {
      input->home_->outs_.push_back(match);
    }
    else if(k->positionNFA_->c_ == Split)
    {
      nextStep(stack, k);
    }
    else
    {
      int i = linSearch(liste, k);
      if(i == -1)
      {
        qGram += k->positionNFA_->c_;
        e = kstate(qGram);
        input->home_->outs_.push_back(e);
        k->home_ = e;
        liste.push_back(k);
        qGram = qGramFrag;
      }
      else
      {
        delete k;
        k = liste[i];
        if(k->home_ != input->home_ && k->home_->start_ == 0)
        {
          input->home_->outs_.push_back(k->home_);
        }
      }
    }
  }


}

std::vector<kState *> nfa2knfa(State* nfa_ptr, const int& q)
{
  /*
  Phase 1:  Schlange befüllen mit keyStates
            zuerst alle vom start durchgehen, mittels stack;
  */
  std::vector<kState *> output{}; //fungiert als start
  std::vector<keyState *> queue{};
  kState* match = kstate("$");

  kState* e;

  //---------------------------------------------------
  State *it_ptr = nfa_ptr;
  try
  {
    firstPhase(it_ptr, queue, q);
  }
  catch(const int &Exception)
  {
    std::cerr<<"QGram zu lang gewählt"<<"\n";
  }

  //Phase 2
  //erstellen der start states und umschreiben der pointer der keys
  std::string edge;
  output.reserve(queue.size());

  for(size_t i = 0; i < queue.size(); i++)
  {
    edge = queue[i]->qGramFrag_;
    edge += queue[i]->positionNFA_->c_;
    int l = linSearchK(output, edge);
    if(l != -1)
    {
        queue[i]->home_ = output[l];
    }
    else
    {
      e = kstate(edge);
      e->start_ = 1;
      queue[i]->home_ = e;
      output.push_back(e);
    }
  }

  //neue keys erstellen und in queue eintragen, sowie kstate erstellen und verknüpfen

  for(size_t i = 0; i < queue.size(); i++)
  {
    nextKeys(queue, queue[i], match);
  }
  //delete the keys
  for(auto b:queue)
  {
    delete b;
  }
  return output;
}

Path* findPath(kState* position)
{
  Path* p = new Path;
  p->qPath_ = 0;
  p->position_ = position;
  return p;
}

/*
 * depht first search, generates the matrix with the possible paths
 */

void dfs(kState* input, std::vector<std::vector<std::string>>& matrix, uint32_t &vector_idx,
          robin_hood::unordered_map<uint64_t, uint32_t> &hash_to_idx,
          std::vector<bitvector> &kmer_bitvex, IndexStructure &ibf)
{
  std::vector<std::string> line{};
  std::stack<Path*> stack{};
  
  size_t qlength = input->qGram_.length();
  auto hash_adaptor = seqan3::views::kmer_hash(seqan3::ungapped{qlength});
  auto && ibf_ref = ibf.getIBF();
  auto agent = ibf_ref.membership_agent();

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

// int main()
// {

//   std::string regex = "ab+c+|.d.";
//   int qlength = 3;
//   State* nfa = post2nfaE(regex);
//   std::vector<kState *> knfa = nfa2knfa(nfa, qlength);
//   std::vector<std::vector<std::string>> matrix{};
//   for(auto i : knfa)
//   {
//     dfs(i,matrix);
//   }

//   for(auto i : matrix)
//   {
//     for(auto j : i)
//     {
//       std::cout<<j<<" ";
//     }
//     std::cout<<"\n";
//   }
//   printGraph(knfa,"out.dot");
//   return 0;
// }


// std::set<std::string> makeLine(std::vector<kState* > qPath)
// {
//   std::set<std::string> line{};
//   std::string s;
//   for(auto e : qPath)
//   {
//     s = e->qGram_;
//     line.insert(s);
//   }
//   //std::sort(line.begin(), line.end());
//   return line;

// }

// std::set<std::set<std::string>> getMatrix(std::vector<kState* > input)
// {
//   std::set<std::set<std::string>> matrix{};

//   std::stack<Path> stack;
//   for(uint i = 0; i < input.size(); i++)
//   {
//     std::vector<kState* > vec{input[i]};
//     stack.push(path(vec, input[i]));
//   }
//   Path p;
//   while(!stack.empty())
//   {
//     p = stack.top();
//     stack.pop();
//     for(auto e : p.position_->outs_)
//     {
//       if(e->qGram_ == "$")
//       {
//         std::set<std::string> line = makeLine(p.qPath_);
//         matrix.insert(line);
//       }
//       else if(std::find(p.qPath_.begin(), p.qPath_.end(), e) == p.qPath_.end())
//       {
//         std::vector<kState* > vec = p.qPath_;
//         vec.push_back(e);
//         stack.push(path(vec, e));
//       }
//     }
//   }
//   return matrix;
// }
