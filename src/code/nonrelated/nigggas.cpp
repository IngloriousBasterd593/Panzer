#include <iostream>
#include <string>
#include <vector>

int main()
{
    std::cout << "oh hel yeah" << std::endl;

    std::vector<std::string> vec;

    vec.push_back("oh");
    vec.push_back("hel");
    vec.push_back("yeah");

    for(auto i : vec)
    {
        std::cout << i << std::endl;
    }



    return 0;
}