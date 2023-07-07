#include <iostream>
#include <fstream>
#include <string>
//#include "fastq.h"

using namespace std;

enum FileReadState {
    IN_HEADER,
    IN_SEQUENCE,
    IN_SEPARATOR,
    IN_QUALITY,
    IN_TRANSITION
};

int main(){
    // Open the fastq file
    //ifstream file("../data/raw_seq_data/ERR020236_1.fastq");
    ifstream file("fq");
    string filename = "fq";

    // Check if the file was successfully opened
    if (!file) {
        cerr << "Unable to open file!" << endl;
        return 1;
    } else {
        cout << "Fichier -" << filename << "- ouvert avec succès." << endl;
    }

    // TEST 3

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
        //cout << i << charNb << endl;
        if (i == charNb){
            //cout << "Condition 1 de la boucle" << endl; // test
            // Suppression du contenu et remplissage
            i = 0;
            beg = file.tellg();
            file.read(buffer, 1024);
            charNb = file.gcount();
        }
        //cout << beg << endl;
        //cout << buffer[i] << endl;
        //cout << charNb << endl;
        //break ;
        while (i < charNb){
            switch(state){
                case IN_HEADER: {
                    //cout << "IN HEADER" << endl;
                    if (buffer[i] == '\n'){
                        ++ num_ligne;
                        //cout << header << endl;
                        //header.clear();
                        state = IN_SEQUENCE;
                    } else {
                        header += buffer[i];
                        //cout << header << endl; // test
                    }
                    break;
                }
                case IN_SEQUENCE: {
                    //cout << "IN SEQUENCE" << endl;
                    ++seqSize;

                    if(isspace(buffer[i])){
                        num_ligne += (buffer[i] == '\n');
                        //cout << "BLA" << endl; // test
                    } else {
                        if (last_c == '\n' && buffer[i] == '+'){
                            state = IN_SEPARATOR;
                            separator += buffer[i];
                            //cout << seq << endl;
                            //cout << seq.length() << endl;
                            //cout << seqSize << endl;
                            //seq.clear();
                        } else {
                            //cerr << "Erreur : " << filename << " - ligne " << num_ligne << endl;
                            //cerr << '\t' << buffer[i] << " : pas un nucléotide" << endl;
                            seq += buffer[i];
                        }
                    }
                    break;
                }
                case IN_SEPARATOR: {
                    //cout << "IN SEPARATOR" << endl;
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
                    //cout << "IN QUALITY" << endl;
                    qual_size++;
                    if (buffer[i] == '\n'){
                        ++num_ligne;
                        state = IN_TRANSITION;
                        if (seqSize -1 != qual_size){ // avant : seqSiver != qual_size -1
                            cout << "Erreur : " << filename << " - " << header << endl;
                            cout << "\tLa qualité et la séquence ne correspondent pas" << endl;
                        }
                        qual_size = 0;
                    }
                    break;
                }
                case IN_TRANSITION:{
                    cout << header << endl;
                    //cout << seq << endl;

                    // OPERATIONS SUR LE READ
                    int nbKmers = seq.length() - kmerSize;
                    for (int i = 0; i <= nbKmers ; i++){
                        //string kmer = seq.substr(i, kmerSize);
                        //cout << kmer << endl;
                        string prefix = seq.substr(i, prefixSize);
                        string suffix = seq.substr(i + prefixSize, suffixSize);
                        cout << prefix << " - " << suffix << endl;
                    }
                    
                    // FIN OPERATIONS SUR LE READ
                    
                    header.clear();
                    seq.clear();
                    seqSize = 0;
                    //cout << "IN TRANSITION" << endl;
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