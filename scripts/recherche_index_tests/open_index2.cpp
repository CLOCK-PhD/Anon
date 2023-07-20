#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <filesystem>

namespace fs = std::filesystem;

// Recherche binaire dans le vecteur sur le premier mot
bool binarySearch(const std::vector<std::string>& lines, const std::string& word) {
    int left = 0;
    int right = lines.size() - 1;

    while (left <= right) {
        int mid = left + (right - left) / 2;
        std::istringstream lineStream(lines[mid]);
        std::string firstWord;
        lineStream >> firstWord;

        if (firstWord == word) {
            return true; // Found the word in the first position of a line
        } else if (firstWord < word) {
            left = mid + 1;
        } else {
            right = mid - 1;
        }
    }

    return false; // Word not found
}

// Charger tous les fichiers dans la map :
// VERSION 1 - PAS DE COMPTAGE DE LIGNE
std::vector<std::string> loadFileContents(const std::string& filePath) {
    std::vector<std::string> lines;
    std::ifstream file(filePath);
    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {
            lines.push_back(line);
        }
        file.close();
    }
    return lines;
}

// Version 2 - Optimisation de la taille des vecteurs
/*std::vector<std::string> loadFileContents(const std::string& filePath) {
    std::ifstream file(filePath);
    if (!file.is_open()) {
        std::cout << "Failed to open file: " << filePath << std::endl;
        return {};
    }

    // Count the number of lines in the file
    int lineCount = 0;
    std::string line;
    while (std::getline(file, line)) {
        lineCount++;
    }

    // Resize the vector to the exact number of lines
    std::vector<std::string> lines(lineCount);

    // Read the file again to populate the vector
    file.clear();
    file.seekg(0, std::ios::beg);
    int index = 0;
    while (std::getline(file, line)) {
        lines[index++] = line;
    }

    file.close();

    return lines;
}*/

std::map<std::string, std::vector<std::string>> loadFilesInDirectory(const std::string& directoryPath) {
    std::map<std::string, std::vector<std::string>> fileMap;

    int fileCount = 0;
    
    for (const auto& entry : fs::directory_iterator(directoryPath)) {
        if (entry.is_regular_file()) {
            std::string filePath = entry.path().string();
            std::vector<std::string> lines = loadFileContents(filePath);
            fileMap[entry.path().filename().string()] = lines;

            std::cout << "Loaded file: " << filePath << " (" << fileCount << ")" <<  std::endl;
            fileCount++;
        }
    }
    
    std::cout << "Total files loaded: " << fileCount << std::endl;

    return fileMap;
}

// Pour extraire d'autres informations de la ligne
std::string extractSecondWord(const std::string& line) {
    std::istringstream lineStream(line);
    std::string firstWord;
    std::string secondWord;
    lineStream >> firstWord >> secondWord;
    return secondWord;
}

// Recherche dichotomique sur l'index
std::string searchFirstWord(const std::vector<std::string>& lines, const std::string& word) {
    int left = 0;
    int right = lines.size() - 1;

    while (left <= right) {
        int mid = left + (right - left) / 2;
        std::istringstream lineStream(lines[mid]);
        std::string firstWord;
        lineStream >> firstWord;
        if (firstWord == word) {
            return extractSecondWord(lines[mid]); // Found the word, extract the second word
        } else if (firstWord < word) {
            left = mid + 1;
        } else {
            right = mid - 1;
        }
    }

    return ""; // Word not found
}

// TEST - Faire la recherche avec une fonction - OK !
void searchWordInFile(const std::map<std::string, std::vector<std::string>> &fileMap, const std::string &key, const std::string &word) {
    auto it = fileMap.find(key);
    if (it != fileMap.end()) {
        const std::vector<std::string>& lines = it->second;
        std::string secondWord = searchFirstWord(lines, word);
        if (!secondWord.empty()) {
            std::cout << "The word '" << word << "' was found in the first position of a line in the file '" << key << "'." << std::endl;
            std::cout << "The second word in the matching line is: " << secondWord << std::endl;
        } else {
            std::cout << "The word '" << word << "' was not found in the first position of a line in the file '" << key << "'." << std::endl;
        }
    } else {
        std::cout << "The specified key '" << key << "' does not exist in the map." << std::endl;
    }
}

int main() {

    /* TEST BASIC POUR L'INDEX EN MAP
    // Créer l'index en tant que map
    // Clé (string): préfix
    // Valeur (vecteur) : chaque ligne du fichier
    std::map<std::string, std::vector<std::string>> fileMap;
    
    std::cout << "Loading files..." << std::endl;
    // Load file contents into vector
    //std::ifstream file1("file1.txt");
    std::ifstream file1("indexDirLink/AAAAA");
    std::vector<std::string> lines1;
    std::string line1;
    while (std::getline(file1, line1)) {
        lines1.push_back(line1);
    }
    fileMap["AAAAA"] = lines1;
    file1.close();

    // Load file contents into vector
    std::ifstream file2("indexDirLink/AAAAC");
    std::vector<std::string> lines2;
    std::string line2;
    while (std::getline(file2, line2)) {
        lines2.push_back(line2);
    }
    fileMap["AAAAC"] = lines2;
    file2.close();

    // Load file contents into vector
    std::ifstream file3("indexDirLink/AAAAG");
    std::vector<std::string> lines3;
    std::string line3;
    while (std::getline(file3, line3)) {
        lines3.push_back(line3);
    }
    fileMap["AAAAG"] = lines3;
    file3.close();
    std::cout << "\tDone." << std::endl;*/


    // TEST - CHARGER TOUS LES FICHIERS DE L'INDEX - OK !!!
    std::string directoryPath = "indexDirLink"; // Specify the directory path
    
    std::map<std::string, std::vector<std::string>> fileMap = loadFilesInDirectory(directoryPath);

    std::string key = "AAAAA"; // Key where you want to search for the word
    std::string word = "TTTTTTTTTTTTTTTTTTTTTTTGCT"; // le mot à chercher

    // TEST AVEC FONCTION BINARY SEARCH - OK
    /*std::cout << "PREMIER TEST" << std::endl;
    auto it = fileMap.find(key);
    if (it != fileMap.end()) {
        const std::vector<std::string>& lines = it->second;
        if (binarySearch(lines, word)) {
            std::cout << "The word '" << word << "' was found in the first position of a line in the file '" << key << "'." << std::endl;
        } else {
            std::cout << "The word '" << word << "' was not found in the first position of a line in the file '" << key << "'." << std::endl;
        }
    } else {
        std::cout << "The specified key '" << key << "' does not exist in the map." << std::endl;
    }*/

    // TEST AVEC SEARCHFIRSTWORD+SECONDWORD - OK !!
    /*std::cout << "DEUXIEME TEST" << std::endl;
    auto it = fileMap.find(key);
    if (it != fileMap.end()) {
        const std::vector<std::string>& lines = it->second;
        std::string secondWord = searchFirstWord(lines, word);
        if (!secondWord.empty()) {
            std::cout << "The word '" << word << "' was found in the first position of a line in the file '" << key << "'." << std::endl;
            std::cout << "The second word in the matching line is: " << secondWord << std::endl;
        } else {
            std::cout << "The word '" << word << "' was not found in the first position of a line in the file '" << key << "'." << std::endl;
        }
    } else {
        std::cout << "The specified key '" << key << "' does not exist in the map." << std::endl;
    }*/

    // TEST - Même chose mais c'est maintenant une fonction - OK !!!
    searchWordInFile(fileMap, key, word);

    return 0;
}
