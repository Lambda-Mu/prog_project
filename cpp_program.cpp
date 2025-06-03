#include <iostream>
#include <fstream>
#include <string>

int main(){
    std::ifstream in("input.txt");
    std::ofstream out("output.txt");

    if(!in or !out){
        std::cerr << "Error opening file!\n";
        return 1;
    }

    std::string data;
    std::getline(in, data);

    for(char& c: data){
        c = toupper(c);
    }
    
    in.close();
    out.close();
}
