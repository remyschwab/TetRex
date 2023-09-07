#include <algorithm>
#include <stack>
#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>
#include "korotkov_nfa.h"
#include "graphMaker.h"

void printState(kState *input)
{
    std::cout<<"----------------------"<<"\n";
    std::cout<<"ID: "<<input<<"\n";
    std::cout<<"q-Gram: "<<input->qGram_<<"\n";
    for(auto e : input->outs_)
    {
        std::cout<<e<<"\n";
    }
    std::cout<<"----------------------"<<"\n";
}

void print(const std::vector<kState *>& input)
{
    std::stack<kState *> stack;
    std::cout<<"Start: "<<"\n";

    for(size_t i = 0; i < input.size(); i++)
    {
        std::cout<<input[i]<<"\n";
        stack.push(input[i]);
    }
    std::cout<<"#################"<<"\n";
    kState* k;
    while(!stack.empty())
    {
        k = stack.top();
        stack.pop();
        if(k->marked_ == 0)
        {
            printState(k);
            k->marked_ = 1;
            for(auto v : k->outs_)
            {
                stack.push(v);
            }
        }
    }
}

std::string printNode(kState *input)
{
    std::stringstream inputS;
    inputS << input;
    std::string line = "\"";
    line += inputS.str();
    if(input->start_)
    {
        line+= "\" [style=filled fillcolor=red label=\"";
    }
    else if(input->qGram_ == "$")
    {
        line+= "\" [shape=doublecircle label=\"";
    }
    else
    {
        line += "\" [label=\"";
    }

    line += input->qGram_;
    line += "\"]\n\"";
    line += inputS.str();
    line += "\"->{";
    for(auto e : input->outs_)
    {
        std::stringstream s;
        s << e;
        line += "\"";
        line += s.str();
        line += "\" ";

    }
    line += "}";
    return line;
}

void printGraph(const std::vector<kState *>& input, const std::string& filename)
{
    std::fstream f;
    f.open(filename, std::ios::out);
    f << "digraph graphname"<<"\n";
    f << "{"<<"\n";

    std::string newline;
    std::stack<kState *> stack;

    for(size_t i = 0; i < input.size(); i++)
    {
        stack.push(input[i]);
    }
    kState* k;
    while(!stack.empty())
    {
        k = stack.top();
        stack.pop();
        if(k->marked_ == 0)
        {
            newline = printNode(k);
            f << newline <<"\n";
            k->marked_ = 1;
            for(auto v : k->outs_)
            {
                stack.push(v);
            }
        }
    }
    f << "}";
    f.close();
}

void uMatrix(std::vector<std::vector<std::string>>& matrix)
{
    for(size_t i = 0; i < matrix.size(); ++i)
    {
        std::sort(matrix[i].begin(), matrix[i].end());
        auto last = std::unique(matrix[i].begin(), matrix[i].end());
        matrix[i].erase(last, matrix[i].end());
    }
    std::sort(matrix.begin(), matrix.end());
    auto last = std::unique(matrix.begin(), matrix.end());
    matrix.erase(last, matrix.end());
}
