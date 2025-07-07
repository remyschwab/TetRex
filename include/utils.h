#pragma once

#include <iostream>
#include <type_traits>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include <fstream>
#include <stack>
#include <filesystem>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <cstring>
#include <cctype>
#include <zlib.h>

#include "hibf/interleaved_bloom_filter.hpp"
#include "hibf/hierarchical_interleaved_bloom_filter.hpp"
#include "hibf/misc/bit_vector.hpp"
#include <seqan3/core/debug_stream.hpp>
#include "kseq.h"


KSEQ_INIT(gzFile, gzread)

#define DBG(x) seqan3::debug_stream << x << std::endl;
#define FOR_SIZE(container) for(size_t i = 0; i < (container).size(); ++i)
#define FLOOP(limit) for(size_t i = 0; i < limit; ++i)

/////////////// Type Declarations ///////////////
using bitvector = seqan::hibf::bit_vector;
/////////////// ****** END ****** ///////////////

// Token types for our lexer
enum class TokenType {
    CHAR,           // Regular character
    DOT,            // . (any character)
    STAR,           // *
    PLUS,           // +
    QUESTION,       // ?
    PIPE,           // |
    LPAREN,         // (
    RPAREN,         // )
    LBRACE,         // {
    RBRACE,         // }
    COMMA,          // ,
    NUMBER,         // Digit sequence
    MINMAX_OP,      // Our custom {m,n} operator
    EXACT_OP,       // Our custom {m} operator
    CHAR_CLASS,     // Character class [abc] or [^abc]
    CONCAT,         // Explicit concatenation operator
    END_OF_INPUT
};

constexpr std::array<char, 20> amino_acids = {
    'A', // Alanine
    'C', // Cysteine
    'D', // Aspartic acid
    'E', // Glutamic acid
    'F', // Phenylalanine
    'G', // Glycine
    'H', // Histidine
    'I', // Isoleucine
    'K', // Lysine
    'L', // Leucine
    'M', // Methionine
    'N', // Asparagine
    'P', // Proline
    'Q', // Glutamine
    'R', // Arginine
    'S', // Serine
    'T', // Threonine
    'V', // Valine
    'W', // Tryptophan
    'Y'  // Tyrosine
};

struct Token {
    TokenType type;
    std::string value;
    int min_count = 0;  // For quantifier tokens
    int max_count = 0;  // For quantifier tokens
    bool is_negated = false;  // For character class tokens
    // std::vector<std::pair<char, char>> char_ranges;  // For character class ranges
    std::vector<char> char_ranges;
    
    Token(TokenType t, const std::string& v = "") : type(t), value(v) {}
    Token(TokenType t, int min_val, int max_val) : type(t), min_count(min_val), max_count(max_val) {}
    Token(TokenType t, bool negated, const std::vector<char>& ranges) 
        : type(t), is_negated(negated), char_ranges(ranges) {}
};

class RegexLexer {
private:
    std::string input;
    size_t pos;
    
    // Parse a number from current position
    int parseNumber() {
        int num = 0;
        while (pos < input.length() && std::isdigit(input[pos])) {
            num = num * 10 + (input[pos] - '0');
            pos++;
        }
        return num;
    }
    
    // Parse character class [abc] or [^abc]
    Token parseCharacterClass() {
        pos++; // Skip '['
        
        if (pos >= input.length()) {
            throw std::runtime_error("Invalid character class: unexpected end of input");
        }
        
        bool is_negated = false;
        if (input[pos] == '^') {
            is_negated = true;
            pos++;
        }
        
        // std::vector<std::pair<char, char>> ranges;
        std::vector<char> ranges;
        // std::string chars;
        
        while (pos < input.length() && input[pos] != ']') {
            char current = input[pos];
            
            if (current == '\\') {
                // Handle escaped character
                pos++;
                if (pos >= input.length()) {
                    throw std::runtime_error("Invalid escape in character class");
                }
                char escaped = input[pos];
                
                // Handle common escape sequences
                switch (escaped) {
                    case 'n': escaped = '\n'; break;
                    case 't': escaped = '\t'; break;
                    case 'r': escaped = '\r'; break;
                    case '\\': escaped = '\\'; break;
                    case ']': escaped = ']'; break;
                    case '^': escaped = '^'; break;
                    case '-': escaped = '-'; break;
                    // Add more escape sequences as needed
                }
                
                ranges.push_back(escaped);
                pos++;
            }
            // else if (pos + 2 < input.length() && input[pos + 1] == '-' && input[pos + 2] != ']') {
            //     // Handle character range like a-z
            //     char start = current;
            //     char end = input[pos + 2];
                
            //     if (start > end) {
            //         throw std::runtime_error("Invalid character range: start > end");
            //     }
                
            //     ranges.push_back({start, end});
            //     pos += 3; // Skip start, '-', and end
            // }
            else {
                // Regular character
                // chars += current;
                ranges.push_back(current);
                pos++;
            }
        }
        
        if (pos >= input.length() || input[pos] != ']')
        {
            throw std::runtime_error("Invalid character class: missing closing ']'");
        }
        
        pos++; // Skip ']'
        
        // Convert individual characters to single-character ranges
        // for (char c : chars)
        // {
        //     ranges.push_back(c);
        // }
        
        if (ranges.empty())
        {
            throw std::runtime_error("Empty character class");
        }
        
        return Token(TokenType::CHAR_CLASS, is_negated, ranges);
    }
    Token parseQuantifier() {
        pos++; // Skip '{'
        
        if (pos >= input.length() || !std::isdigit(input[pos])) {
            throw std::runtime_error("Invalid quantifier: expected number after '{'");
        }
        
        int min_val = parseNumber();
        
        if (pos >= input.length()) {
            throw std::runtime_error("Invalid quantifier: unexpected end of input");
        }
        
        if (input[pos] == '}') {
            // {m} - exact quantifier
            pos++; // Skip '}'
            return Token(TokenType::EXACT_OP, min_val, min_val);
        } else if (input[pos] == ',') {
            pos++; // Skip ','
            
            if (pos >= input.length()) {
                throw std::runtime_error("Invalid quantifier: unexpected end after ','");
            }
            
            if (input[pos] == '}') {
                // {m,} - min quantifier (not handling this case for now)
                throw std::runtime_error("Open-ended quantifiers {m,} not supported");
            }
            
            if (!std::isdigit(input[pos])) {
                throw std::runtime_error("Invalid quantifier: expected number after ','");
            }
            
            int max_val = parseNumber();
            
            if (pos >= input.length() || input[pos] != '}') {
                throw std::runtime_error("Invalid quantifier: expected '}' after max value");
            }
            
            pos++; // Skip '}'
            
            if (min_val > max_val) {
                throw std::runtime_error("Invalid quantifier: min > max");
            }
            
            return Token(TokenType::MINMAX_OP, min_val, max_val);
        } else {
            throw std::runtime_error("Invalid quantifier: expected ',' or '}' after min value");
        }
    }

public:
    RegexLexer(const std::string& regex) : input(regex), pos(0) {}
    
    std::vector<Token> tokenize() {
        std::vector<Token> tokens;
        
        while (pos < input.length()) {
            char c = input[pos];
            
            switch (c) {
                case '.':
                    tokens.push_back(Token(TokenType::DOT, "."));
                    pos++;
                    break;
                case '*':
                    tokens.push_back(Token(TokenType::STAR, "*"));
                    pos++;
                    break;
                case '+':
                    tokens.push_back(Token(TokenType::PLUS, "+"));
                    pos++;
                    break;
                case '?':
                    tokens.push_back(Token(TokenType::QUESTION, "?"));
                    pos++;
                    break;
                case '|':
                    tokens.push_back(Token(TokenType::PIPE, "|"));
                    pos++;
                    break;
                case '(':
                    tokens.push_back(Token(TokenType::LPAREN, "("));
                    pos++;
                    break;
                case ')':
                    tokens.push_back(Token(TokenType::RPAREN, ")"));
                    pos++;
                    break;
                case '[':
                    tokens.push_back(parseCharacterClass());
                    break;
                case '{':
                    tokens.push_back(parseQuantifier());
                    break;
                case '\\':
                    // Handle escaped characters
                    pos++; // Skip backslash
                    if (pos >= input.length()) {
                        throw std::runtime_error("Invalid escape: end of input after '\\'");
                    }
                    tokens.push_back(Token(TokenType::CHAR, std::string(1, input[pos])));
                    pos++;
                    break;
                default:
                    // Regular character
                    tokens.push_back(Token(TokenType::CHAR, std::string(1, c)));
                    pos++;
                    break;
            }
        }
        
        tokens.push_back(Token(TokenType::END_OF_INPUT));
        return tokens;
    }
};

class PostfixConverter {
private:
    // Check if we need to insert concatenation between two tokens
    bool needsConcatenation(const Token& current, const Token& previous) {
        // After these tokens, we might need concatenation
        bool after_operand = (previous.type == TokenType::CHAR || 
                             previous.type == TokenType::DOT ||
                             previous.type == TokenType::CHAR_CLASS ||
                             previous.type == TokenType::RPAREN);
        bool after_quantifier = (previous.type == TokenType::STAR ||
                               previous.type == TokenType::PLUS ||
                               previous.type == TokenType::QUESTION ||
                               previous.type == TokenType::MINMAX_OP ||
                               previous.type == TokenType::EXACT_OP);
        
        // Before these tokens, we might need concatenation
        bool before_operand = (current.type == TokenType::CHAR ||
                              current.type == TokenType::DOT ||
                              current.type == TokenType::CHAR_CLASS ||
                              current.type == TokenType::LPAREN);
        
        return (after_operand || after_quantifier) && before_operand;
    }
    
    // Get operator precedence
    int getPrecedence(TokenType type) {
        switch (type) {
            case TokenType::PIPE: return 1;  // Alternation (lowest)
            case TokenType::CONCAT: return 2;  // Concatenation
            case TokenType::STAR:
            case TokenType::PLUS:
            case TokenType::QUESTION:
            case TokenType::MINMAX_OP:
            case TokenType::EXACT_OP:
                return 3;  // Repetition (highest)
            default: return 0;
        }
    }
    
    // Check if token is an operator
    bool isOperator(TokenType type) {
        return type == TokenType::PIPE || type == TokenType::CONCAT ||
               type == TokenType::STAR || type == TokenType::PLUS ||
               type == TokenType::QUESTION || type == TokenType::MINMAX_OP ||
               type == TokenType::EXACT_OP;
    }
    
    // Convert token to postfix string representation
    std::string tokenToPostfix(const Token& token) {
        switch (token.type) {
            case TokenType::CHAR:
                return token.value;
            case TokenType::DOT:
                return "FQ|L|T|K|P|A|Y|R|N|H|G|E|C|I|V|D|W|S|M|"; // Wildcards need to be expanded at some point I guess
            case TokenType::CHAR_CLASS: {
                std::string result = "";
                if(token.is_negated)
                {
                    Token& mutable_t = const_cast<Token&>(token); // Careful here!
                    std::sort(mutable_t.char_ranges.begin(), mutable_t.char_ranges.end());
                    std::vector<char> difference;
                    std::set_difference(amino_acids.begin(), amino_acids.end(), mutable_t.char_ranges.begin(), mutable_t.char_ranges.end(), std::back_inserter(difference));
                    result += difference[0];
                    for(size_t i = 1; i < difference.size(); ++i)
                    {
                        result += difference[i];
                        result += '|';
                    }
                    return result;
                }
                result += token.char_ranges[0];
                for(size_t i = 1; i < token.char_ranges.size(); ++i)
                {
                    result += token.char_ranges[i];
                    result += '|';
                }
                return result;
            }
            case TokenType::STAR:
                return "*";
            case TokenType::PLUS:
                return "+";
            case TokenType::QUESTION:
                return "?";
            case TokenType::PIPE:
                return "|";
            case TokenType::CONCAT:
                return "-";
            case TokenType::EXACT_OP:
                return "{" + std::to_string(token.min_count) + "}";
            case TokenType::MINMAX_OP:
                return "{" + std::to_string(token.min_count) + "," + 
                       std::to_string(token.max_count) + "}";
            default:
                return "";
        }
    }

public:
    // Convert infix regex to postfix with custom quantifier operators
    std::string infixToPostfix(const std::string& regex) {
        RegexLexer lexer(regex);
        std::vector<Token> tokens = lexer.tokenize();
        
        // Insert explicit concatenation operators
        std::vector<Token> withConcat;
        for (size_t i = 0; i < tokens.size(); i++) {
            if (i > 0 && needsConcatenation(tokens[i], tokens[i-1])) {
                withConcat.push_back(Token(TokenType::CONCAT));
            }
            withConcat.push_back(tokens[i]);
        }
        
        // Convert to postfix using Shunting Yard algorithm
        std::string result;
        std::stack<Token> operators;
        
        for (const Token& token : withConcat) {
            if (token.type == TokenType::CHAR || token.type == TokenType::DOT || 
                token.type == TokenType::CHAR_CLASS) {
                result += tokenToPostfix(token);
            }
            else if (token.type == TokenType::LPAREN) {
                operators.push(token);
            }
            else if (token.type == TokenType::RPAREN) {
                while (!operators.empty() && operators.top().type != TokenType::LPAREN) {
                    result += tokenToPostfix(operators.top());
                    operators.pop();
                }
                if (!operators.empty()) {
                    operators.pop();  // Remove the '('
                }
            }
            else if (isOperator(token.type)) {
                while (!operators.empty() && 
                       operators.top().type != TokenType::LPAREN &&
                       getPrecedence(operators.top().type) >= getPrecedence(token.type)) {
                    result += tokenToPostfix(operators.top());
                    operators.pop();
                }
                operators.push(token);
            }
            else if (token.type == TokenType::END_OF_INPUT) {
                break;
            }
        }
        
        // Pop remaining operators
        while (!operators.empty()) {
            result += tokenToPostfix(operators.top());
            operators.pop();
        }
        
        return result;
    }
};

std::string translate(const std::string& pattern);
