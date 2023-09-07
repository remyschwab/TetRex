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



// Anonymous namespace (avoids naming conflicts with other translation units)
namespace {

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
        e1 = new keyState{qGram, itptr->out1_, kptr};
        stack.push(e1);
        break;
    case SplitA:
        e1 = new keyState{qGram, itptr->out1_, kptr};
        e2 = new keyState{qGram, itptr->out2_, kptr};
        stack.push(e2);
        stack.push(e1);
        break;
    case SplitB:
        e1 = new keyState{qGram, itptr->out1_, kptr};
        e2 = new keyState{qGram, itptr->out2_, kptr};
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
        if(k->positionNFA_->c_ == Match) // The match state of the NFA can only be reached by Phase1 by error
        {
            delete k;
            throw int();
        }
        if(k-> qGramFrag_.size() == q-1 && k->positionNFA_->c_ != SplitA && k->positionNFA_->c_ != SplitB) // Phase1 only collects k-1mers
        {
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
int linSearchK(const std::vector<kState *>& liste, std::string obj) // Looks for a kState by kmer
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

int linSearch(const std::vector<keyState *>& liste, keyState* obj) // Searches for keyState by kmer AND NFA pointer
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
 * is similar to oneStep
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
    case SplitA:
            e1 = new keyState{input->qGramFrag_, itptr->out1_, nullptr};
            e2 = new keyState{input->qGramFrag_, itptr->out2_, nullptr};
            stack.push(e2);
            stack.push(e1);
            break;
    case SplitB:
            e1 = new keyState{input->qGramFrag_, itptr->out1_, nullptr};
            e2 = new keyState{input->qGramFrag_, itptr->out2_, nullptr};
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
    std::string qGramFrag = input->home_->qGram_; // Get the kmer from the kState (and not the NFA)
    qGramFrag = qGramFrag.substr(1);
    std::string qGram = qGramFrag;

    k = new keyState{qGramFrag, input->positionNFA_->out1_, nullptr};
    nextStep(stack, k); // I think this will always just push k on the stack... See condition at line 67

    while(!stack.empty())
    {
        k = stack.top();
        stack.pop();
        if(k->positionNFA_->c_== Match)
        {
            input->home_->outs_.push_back(match); // Completed a path through kNFA
        }
        else if(k->positionNFA_->c_ == SplitA || k->positionNFA_->c_ == SplitB)
        {
            nextStep(stack, k); // Unpack a split
        }
        else
        {
            int i = linSearch(liste, k);
            if(i == -1)
            {
                qGram += k->positionNFA_->c_;
                e = new kState{qGram};
                input->home_->outs_.push_back(e);
                k->home_ = e;
                liste.push_back(k); // The array being modified here is also being iterated over outside this scope
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

}


std::vector<kState *> nfa2knfa(State* nfa_ptr, const int& q)
{
    /*
    Phase 1:
        - Fill queue with keyStates
        - First fill stack with keystates reachable from start
    */
    std::vector<kState *> output{}; //fungiert als start
    std::vector<keyState *> queue{};
    kState* match = new kState{"$"};

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
    //create the start states and rewrite the pointers of the keys
    // Initialize a list of starting kStates and set the kStates of the keyStates
    std::string edge;
    output.reserve(queue.size());

    for(size_t i = 0; i < queue.size(); i++)
    {
        edge = queue[i]->qGramFrag_;
        edge += queue[i]->positionNFA_->c_;
        int l = linSearchK(output, edge);
        if(l != -1) // If there is already a kState that has this kmer
        {
            queue[i]->home_ = output[l]; // Set the kState component of the keyStates in the queue
        }
        else // If there's no kState with this qgram in the output, then add it as the start of a path
        {
            e = new kState{edge};
            e->start_ = 1;
            queue[i]->home_ = e;
            output.push_back(e);
        }
    }
    //neue keys erstellen und in queue eintragen, sowie kstate erstellen und verknüpfen
    for(size_t i = 0; i < queue.size(); i++) // This is very confusing: see line 168
        nextKeys(queue, queue[i], match);
    //delete the keys
    for(auto b:queue)
        delete b;
    return output;
}

Path* findPath(kState* position)
{
    Path* p = new Path;
    p->qPath_ = 0;
    p->position_ = position;
    return p;
}
