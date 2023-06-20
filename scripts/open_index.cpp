#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <filesystem>

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

int main(){

    std::string directoryPath = "indexDirLink"; // Répertoire contenant l'index

    std::map<std::string, std::string> fileMap; // Map - Clé : nom du fichier ; Valeur : contenu du fichier


    // TEST 1 - Ajout des fichiers manuellement
    /*
    std::string file1 = "indexDirLink/AAAAA";
    std::string file2 = "indexDirLink/AAAAT";
    std::string file3 = "indexDirLink/AAAAG";
    // Charer en mémoire le contenu des fichiers dans les clés des indexes
    fileMap["AAAAA"] = readFile(file1);
    fileMap["AAAAT"] = readFile(file2);
    fileMap["AAAAG"] = readFile(file3);*/

    // TEST 2 - Ajout de tous les fichiers depuis le répertoire
    int loadingCounter = 0;
    for (const auto& entry : std::filesystem::directory_iterator(directoryPath)){
        if(entry.is_regular_file()){
            std::string filename = entry.path().filename().string();
            std::string filepath = entry.path().string();
            std::cout << "Loading file: " << filename << "(" << loadingCounter <<")" << std::endl;
            fileMap[filename] = readFile(filepath);
        }
    }
    

    // Accéder au contenu d'un fichier - pour voir si ça marche bien
    std::cout << "Contenu du fichier AAAAA:" << std::endl;
    std::cout << fileMap["AAAAA"] << std::endl;

    return 0;
}