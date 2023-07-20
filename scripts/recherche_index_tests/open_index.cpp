#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <filesystem>
#include <sstream>
#include <vector>
#include <algorithm>

// TEST POUR OUVRIR ET CHARGER L'INDEX EN MEMOIRE
// Compilation : C++17 pour filesystem

// Fonction pour charger un fichier en tant que valeur d'une clé de map
std::string readFile(const std::string &filename){
    std::ifstream file(filename);
    if (!file){
        std::cerr << "Failed to open the file: " << filename << std::endl;
        return "";
    }

    std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    file.close();

    return content;
}

// TEST RECHERCHE
bool searchFirstWord(const std::string &content, const std::string &word){
    std::istringstream iss(content);
    std::string line;
    while(std::getline(iss, line)){
        std::istringstream lineStream(line);
        std::string firstWord;
        if (lineStream >> firstWord && firstWord == word){
            return true;
        }
    }
    return false;
}

// TEST RECHERCHE 2
bool searchFirstWord2(const std::string &content, const std::string &word){
    std::istringstream iss(content);
    std::string line;
    std::vector<std::string> lines;

    while (std::getline(iss, line)){
        lines.push_back(line);
    }

    auto it = std::lower_bound(lines.begin(), lines.end(), word);

    if (it != lines.end() && *it == word){
        return true;
    }

    return false;
}



// TEST RECHERCHE 4
bool binarySearch(const std::string& line, const std::string& word) {
    std::istringstream lineStream(line);
    std::string firstWord;
    lineStream >> firstWord;

    return (firstWord == word);
}

bool searchFirstWord4(const std::string& content, const std::string& word) {
    std::istringstream iss(content);
    std::string line;

    while (std::getline(iss, line)) {
        if (binarySearch(line, word)) {
            return true; // Found the word in the first position of a line
        }
    }

    return false; // Word not found
}

int main(){

    std::string directoryPath = "indexDirLink"; // Répertoire contenant l'index

    std::map<std::string, std::string> fileMap; // Map - Clé : nom du fichier ; Valeur : contenu du fichier


    // TEST 1 - Ajout des fichiers manuellement - OK
    // TEST 3 - On utilise ça pour faire les tests pour la recherche, c'est plus court.
    
    std::string file1 = "indexDirLink/AAAAA";
    //std::string file2 = "indexDirLink/AAAAT";
    //std::string file3 = "indexDirLink/AAAAG";
    // Charer en mémoire le contenu des fichiers dans les clés des indexes
    std::cout << "Chargement des fichiers..." << std::endl;
    fileMap["AAAAA"] = readFile(file1);
    //fileMap["AAAAT"] = readFile(file2);
    //fileMap["AAAAG"] = readFile(file3);
    std::cout << "\tFin." << std::endl;

    // TEST 2 - Ajout de tous les fichiers depuis le répertoire - OK
    // Prends un peu de temps, on utilise le code d'avant pour la partie recherche
    /*int loadingCounter = 0;
    for (const auto& entry : std::filesystem::directory_iterator(directoryPath)){
        if(entry.is_regular_file()){
            std::string filename = entry.path().filename().string();
            std::string filepath = entry.path().string();
            std::cout << "Loading file: " << filename << "(" << loadingCounter <<")" << std::endl;
            fileMap[filename] = readFile(filepath);
        }
    }*/

    // Accéder au contenu d'un fichier - pour voir si ça marche bien - OK TEST 1 & 2
    //std::cout << "Contenu du fichier AAAAA:" << std::endl;
    //std::cout << fileMap["AAAAA"] << std::endl;

    // TEST - Recherche d'un mot dans une clé spécifique
    std::cout << "Début de la recherche" << std::endl;
    std::string key = "AAAAA";  // Préfixe
    std::string word = "TTTTTTTTTTTTTTTTTTTTTTTTCT"; // Suffixe

    auto it = fileMap.find(key);
    if (it != fileMap.end()) {
        const std::string &content = it->second;
        // if (searchFirstWord(content, word)){
        if (searchFirstWord4(content, word)){
            std::cout << key << word << " - trouvé " << std::endl;
        } else {
            std::cout << "Pas trouvé" << std::endl;
        }
    } else {
        std::cout << "La clé spécifiée '" << key << "' n'existe pas." << std::endl;
    }

    // Test binary search
    /*std::string content = "word1\tbli\nword2\tbla\nword3\tmocos\nword4\tblu\n";
    std::string word3 = "word3";

    if (binarySearch(content, word3)) {
        std::cout << "The word '" << word3 << "' was found in the first position of a line in the content." << std::endl;
    } else {
        std::cout << "The word '" << word3 << "' was not found in the first position of a line in the content." << std::endl;
    }*/

    return 0;
}