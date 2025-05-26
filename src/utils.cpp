#include "utils.h"

std::string translate(const std::string& pattern)
{
    PostfixConverter converter;
    std::string postfix;
    try
    {
        postfix = converter.infixToPostfix(pattern);
    }
    catch (const std::exception& e)
    {
        std::cout << "Error: " << e.what() << "\n---\n";
    }
    return postfix;
}