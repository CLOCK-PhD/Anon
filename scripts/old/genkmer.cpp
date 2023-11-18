#include <iostream>
#include <string>

bool isDegenerate(char base) {
    // Define a function to check if a base is degenerate
    return (base == 'R' || base == 'Y' || base == 'S' || base == 'W' || base == 'K' ||
            base == 'M' || base == 'B' || base == 'D' || base == 'H' || base == 'V' || base == 'N');
}

int main() {
    std::string sequence = "ATCGRYTAGCTYAAAAATTCGGRASTTAYGCGTY"; // Replace this with your sequence
    int k = 3; // Set the desired k-mer length

    int sequenceLength = sequence.length();

    std::cout << "Référence" << std::endl;
    for (int i = 0; i <= sequenceLength - k; i++) {
        // Check if the current k-mer contains any degenerate nucleotides
        bool containsDegenerate = false;
        for (int j = i; j < i + k; j++) {
            if (isDegenerate(sequence[j])) {
                containsDegenerate = true;
                break;
            }
        }

        if (containsDegenerate) {
            // Skip this k-mer if it contains a degenerate nucleotide
            continue;
        }

        std::string kmer = sequence.substr(i, k);

        // Convert characters to uppercase before outputting
        for (char& c : kmer) {
            c = std::toupper(c);
        }

        std::cout << "K-mer " << i + 1 << ": " << kmer << std::endl;
    }

    std::cout << "TEST" << std::endl;
    int i = 0;
    while(i <= sequenceLength - k){
        bool containsDegenerate = false;
        for (int j = i; j < i + k; j++) {
            if (isDegenerate(sequence[j])) {
                containsDegenerate = true;
                i = j+1;
                break;
            }
        }
        if (containsDegenerate) {
            // Skip this k-mer if it contains a degenerate nucleotide
            continue;
        }
        std::string kmer = sequence.substr(i, k);
        // Convert characters to uppercase before outputting
        for (char& c : kmer) {
            c = std::toupper(c);
        }
        std::cout << "K-mer " << i + 1 << ": " << kmer << std::endl;
        i++;
    }

    return 0;
}
