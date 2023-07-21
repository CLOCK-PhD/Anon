#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <filesystem>

/* Programme de test pour la recherche des k-mers dans l'index

Réalisé pour faire des tests et éviter de faire n'importe quoi dans KIM.
Pour les tests, les accès aux fichiers sont en dur dans le code.

Pas mal d'optimisations à envisager:
- utiliser des arrays à la place de vecteur pour stocker les lignes dans le map qui contient l'index
(pour gagner en mémoire et en vitesse de recherche ?)
(implique de parcourir les fichiers de l'index une première fois - penser à supprimer les arrays à la fin du programme)
- otpimiser la lecture du fichier fastq (pour l'instant, c'est plus ou moins le projet de M1 adapté)
- voir s'il n'y a pas une méthode plus efficace et rapide pour exclure les k-mers qui contiennent un N
- envisager une extraction des informations  plus rapide sur les lignes où un k-mer est trouvé
- Renommer les fonctions avec des noms plus explicites
- inclure tdqm pour faire des barres de progressions
*/

// Copie de search_index.cpp, les commentaires et tests en moins.
// Ajout du programme open_index2.cpp pour charger l'index en mémoire et faire les opérations.
// Les opérations se déroulent dans la partie "IN TRANSITION", quand on a récupéré toutes les infos du reads.

using namespace std;
namespace fs = std::filesystem;

enum FileReadState {
    IN_HEADER,
    IN_SEQUENCE,
    IN_SEPARATOR,
    IN_QUALITY,
    IN_TRANSITION
};

// CHARGEMENT DE L'INDEX
// Charger un fichier dans un vecteur et le mettre dans l'index
std::vector<std::string> loadFileContents(const std::string &filePath) {
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

// Créer la map qui contient l'index
std::map<std::string, std::vector<std::string>> loadFilesInDirectory(const std::string &directoryPath) {
    std::map<std::string, std::vector<std::string>> fileMap;

    int fileCount = 0;
    
    for (const auto &entry : fs::directory_iterator(directoryPath)) {
        if (entry.is_regular_file()) {
            std::string filePath = entry.path().string();
            std::vector<std::string> lines = loadFileContents(filePath);
            fileMap[entry.path().filename().string()] = lines;
            // Message d'indication pour le chargement
            std::cout << "Loaded file: " << filePath << " (" << fileCount << ")" <<  std::endl;
            fileCount++;
        }
    }
    // Récap du nombre de fichiers chargés
    std::cout << "Total files loaded: " << fileCount << std::endl;

    return fileMap;
}

// RECHERCHE DANS L'INDEX
// Pour extraire d'autres informations de la ligne qui a trouvé le k-mer dans l'index
std::string extractSecondWord(const std::string& line) {
    std::istringstream lineStream(line);
    std::string firstWord;
    std::string secondWord;
    lineStream >> firstWord >> secondWord;
    return secondWord;
}

// Recherche dichotomique dans le vecteur sur le premier mot
std::string searchFirstWord(const std::vector<std::string>& lines, const std::string& word) {
    int left = 0;
    int right = lines.size() - 1;

    while (left <= right) {
        int mid = left + (right - left) / 2;
        std::istringstream lineStream(lines[mid]);
        std::string firstWord;
        lineStream >> firstWord;
        if (firstWord == word) {
            return extractSecondWord(lines[mid]); // Une fois trouvé, on extrait le deuxième mot sur la ligne (rsid)
        } else if (firstWord < word) {
            left = mid + 1;
        } else {
            right = mid - 1;
        }
    }

    return ""; // Mot non trouvé
}


// Fonction réalisant la recherche et retournant le rsid
std::string searchWordInFile(const std::map<std::string, std::vector<std::string>> &fileMap, const std::string &key, const std::string &word) {
    auto it = fileMap.find(key);
    if (it != fileMap.end()) {
        const std::vector<std::string>& lines = it->second;
        std::string secondWord = searchFirstWord(lines, word);
        if (!secondWord.empty()) {
            // C'est ici qu'on fait les choses quand on trouve le k-mer
            //std::cout << "The word '" << word << "' was found in the first position of a line in the file '" << key << "'." << std::endl;
            //std::cout << "The second word in the matching line is: " << secondWord << std::endl;
            //std::cout << key << word << '\t' << secondWord << std::endl;
            return secondWord;
        }/* else {
            std::cout << "The word '" << word << "' was not found in the first position of a line in the file '" << key << "'." << std::endl;
        }*/
    } /*else {
        std::cout << "The specified key '" << key << "' does not exist in the map." << std::endl;
    }*/
}

// Vérifier la présence d'un N dans la séquence
bool containsCharacter(const std::string &str, char targetChar){
    for(char c : str){
        if (c == targetChar){
            return true;
        }
    }
    return false;
}

int main(){

    // Récupérer les résultats de la recherche
    std::map<std::string, std::vector<string>> resultMap;

    /////////////////////////////////////////
    //               INDEX                 //
    /////////////////////////////////////////

    // Charger tous les fichiers de l'index
    //std::string directoryPath = "indexDirLink"; // Chemin d'accès du dossier, version longue
    std::string directoryPath = "../data/inde_lite_redux";
    std::map<std::string, std::vector<std::string>> fileMap = loadFilesInDirectory(directoryPath);

    /////////////////////////////////////////
    //               FASTQ                 //
    /////////////////////////////////////////

    //ifstream file("../data/raw_seq_data/ERR020236_1.fastq");
    ifstream file("fq"); // Chemin d'accès du fichier fastq
    string filename = "fq"; // Un peu inutile, j'avais fait ça quand je n'avais pas fait de lien dans mon dossier où est le script; à virer

    // Vérifier si le fichier s'est bien ouvert
    if (!file) {
        cerr << "Unable to open file!" << endl;
        return 1;
    } else {
        cout << "Fichier -" << filename << "- ouvert avec succès." << endl;
    }

    // Déclaration de variables pour les opérations sur les kmers
    int kmerSize = 31;
    int prefixSize = 5;
    int suffixSize = kmerSize - prefixSize;

    // Déclaration de variables pour la lecture du fichier fastq
    string header;
    string seq;
    string separator;
    string qual;
    char last_c = '\n';
    size_t beg = 0; // beginning
    size_t seqSize = 0; // Taille de la séquence
    size_t qual_size = 0;   // taille de la séquence qualité
    size_t i = 0;   // Notre itérateur préféré
    size_t charNb = 0;  // nombre de caractères dans le buffer
    size_t num_ligne = 1; // Numéro de la ligne ; pour le message d'erreur
    char buffer[1024];
    FileReadState state = IN_HEADER;

    // tant que le nombre de char dans le buffer n'atteint pas 0, on reste dans la boucle
    do {
        if (i == charNb){
            // Suppression du contenu et remplissage
            i = 0;
            beg = file.tellg();
            file.read(buffer, 1024);
            charNb = file.gcount();
        }
        while (i < charNb){
            switch(state){
                case IN_HEADER: {
                    if (buffer[i] == '\n'){
                        ++ num_ligne;
                        state = IN_SEQUENCE;
                    } else {
                        header += buffer[i];
                    }
                    break;
                }
                case IN_SEQUENCE: {
                    ++seqSize;

                    if(isspace(buffer[i])){
                        num_ligne += (buffer[i] == '\n');
                    } else {
                        if (last_c == '\n' && buffer[i] == '+'){
                            state = IN_SEPARATOR;
                            separator += buffer[i];
                        } else {
                            //cerr << "Erreur : " << filename << " - ligne " << num_ligne << endl;
                            //cerr << '\t' << buffer[i] << " : pas un nucléotide" << endl;
                            seq += buffer[i];
                        }
                    }
                    break;
                }
                case IN_SEPARATOR: {
                    if (buffer[i] == '\n'){
                        ++num_ligne;
                        separator.clear();
                        state = IN_QUALITY;
                    } else {
                        separator += buffer[i];
                    }
                    break;
                }
                case IN_QUALITY: {
                    qual_size++;
                    if (buffer[i] == '\n'){
                        ++num_ligne;
                        state = IN_TRANSITION;
                        if (seqSize -1 != qual_size){ // avant : seqSize != qual_size -1
                            cout << "Erreur : " << filename << " - " << header << endl;
                            cout << "\tLa qualité et la séquence ne correspondent pas" << endl;
                        }
                        qual_size = 0;
                    }
                    break;
                }
                case IN_TRANSITION:{
                    //cout << header << endl;
                    //cout << seq << endl;
                    std::string rsid;

                    ////////////////////////////
                    // OPERATIONS SUR LE READ //
                    ////////////////////////////
                    int nbKmers = seq.length() - kmerSize;
                    for (int i = 0; i <= nbKmers ; i++){
                        //string kmer = seq.substr(i, kmerSize);
                        //cout << kmer << endl;
                        string prefix = seq.substr(i, prefixSize);
                        string suffix = seq.substr(i + prefixSize, suffixSize);
                        //cout << prefix << " - " << suffix << endl;
                        
                        // Test sans vérification que le suffix contient un N
                        //searchWordInFile(fileMap, prefix, suffix);

                        // Si le suffixe ne contient pas de N, on fait la recherche
                        if (!containsCharacter(suffix, 'N')){
                            rsid = searchWordInFile(fileMap, prefix, suffix);
                            resultMap.emplace(rsid, header);
                        }
                    }
                    ////////////////////////////////
                    // FIN OPERATIONS SUR LE READ //
                    ////////////////////////////////
                    
                    // Nettoyage pour le prochain read
                    header.clear();
                    seq.clear();
                    seqSize = 0;
                    // Préparation pour la lecture du prochain read
                    if (last_c == '\n' && buffer[i] == '@'){
                        state = IN_HEADER;
                        // On écrit dans le header
                        header += buffer[i];
                    } else if (buffer[i] == '\n') {
                        ++num_ligne;
                    }
                }
            }
            last_c = buffer[i++];
        }

    }while(charNb);


    if (!header.empty()){
        header.clear();
    }

    // Close the file
    file.close();

    cout << "le chat" << endl;

    return 0;
}