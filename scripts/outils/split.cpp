#include <iostream>
#include <sstream>
#include <vector>
#include <string>

using namespace std;

std::vector<std::string> split(const std::string &s, char delimiter) {
    std::vector<std::string> tokens;
    std::istringstream stream(s);
    std::string token;

    while (std::getline(stream, token, delimiter)) {
        tokens.push_back(token);
    }

    return tokens;
}

int main() {
    std::string input = "1000Genomes:0.9399,0.06012,.|ALSPAC:0.9995,0.0005189,.|GnomAD:0.9376,0.06242,.|GoNL:0.999,0.001002,.|HapMap:0.855,0.145,.|Qatari:0.9213,0.0787,.|SGDP_PRJ:0.425,0.575,.|TOPMED:0.9351,0.06487,.|TWINSUK:0.9992,0.0008091,.|dbGaP_PopFreq:0.9312,0.06877,0";
    char delimiter = '|'; // Split each source and its frequencies from the others

    std::vector<std::string> tokens = split(input, delimiter);

    for (int i=0; i < tokens.size(); i++){
        cout << tokens[i] << endl;
        std::vector<std::string> subtokens = split(tokens[i], ':');
        string sourcename = subtokens[0];
        string allele_frequencies = subtokens[1];
        cout << sourcename << endl;
        cout << '\t' << allele_frequencies << endl;
        std::vector<string> frequencies = split(allele_frequencies, ',');
        for(int j=0; j < frequencies.size(); j++){
            //string test = frequencies[j];
            if (frequencies[j] == "."){
                frequencies[j] = (float)'0.0';
            }
            cout << '\t' << frequencies[j] << endl;
        }
        
    }

    return 0;
}
