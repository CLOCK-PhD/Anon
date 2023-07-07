#include <fstream>
#include <iostream>
#include "fastq.h"

enum FileReadState {
            IN_HEADER,
            IN_SEQUENCE,
            IN_SEPARATOR,
            IN_QUALITY,
            IN_TRANSITION
        };

// CONSTRUCTEURS
Fastq::Fastq():
    header("Inconnu"),
    seq_pos(),
    seq_size(){}

Fastq::Fastq(std::string header, size_t seq_pos, size_t seq_size, size_t qual_pos, size_t qual_size):
    header(header),
    seq_pos(seq_pos),
    seq_size(seq_size),
    qual_pos(qual_pos),
    qual_size(qual_size){}

// GETTERS
std::string Fastq::getHeader() const {
    return header;
}

size_t Fastq::getSeq_pos() const {
    return seq_pos;
}

size_t Fastq::getSeq_size() const {
    return seq_size;
}

size_t Fastq::getQual_pos() const {
    return qual_pos;
}

size_t Fastq::getQual_size() const {
    return qual_size;
}

// SETTERS
void Fastq::setHeader(std::string header){
    this->header = header;
}

void Fastq::setSeq_pos(size_t seq_pos){
    this->seq_pos = seq_pos;
}

void Fastq::setSeq_size(size_t seq_size){
    this->seq_size = seq_size;
}

void Fastq::setQual_pos(size_t qual_pos){
    this->qual_pos = qual_pos;
}

void Fastq::setQual_size(size_t qual_size){
    this->qual_size = qual_size;
}

// Autres
bool Fastq::isBlank(char c){
    return ((c == ' ') || (c == '\t') || (c == '\n'));
}

// Je pense que c'est pas bon du tout et que c'est pas du tout ce que je veux faire
std::string Fastq::getSeq(const std::string &the_file){
    std::ifstream file(the_file.c_str(), std::ios::in);
    char seq_buffer[1024];
    size_t charNb;
    size_t length = getSeq_size();
    file.seekg(getSeq_pos());
    std::string seq(length, '?');
    size_t p = 0;
    size_t i;

    do{
        if(i == charNb){
            i = 0;
            file.read(seq_buffer, 1024);
            charNb = file.gcount();
        }

        file.read(seq_buffer, 1024);
        charNb = file.gcount();
        for(size_t i = 0; (i < charNb) && (p < length); i++){
            seq[p++] = seq_buffer[i];
        }
    }while(charNb);

    file.close();
    return seq;
}

static void lecture(std::string &the_file){
    std::ifstream file(the_file.c_str(), std::ios::in);
    Fastq fq;
    std::string header;
    std::string separator;
    std::string qual;
    char last_c = '\n';
    size_t beg = 0;
    size_t seqSize = 0;
    size_t qualSize = 0;
    size_t i = 0;
    size_t charNb = 0;
    size_t num_ligne = 1;
    char buffer[1024];
    FileReadState state = IN_HEADER;

    // Tant que le nombre de caractères dans le buffer n'atteing pas 0, on est dans la boucle
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
                        fq.setHeader(header);
                        std::cout << header << std::endl;
                        header.clear();
                        state = IN_SEQUENCE;
                        fq.setSeq_pos(beg + i + 1);
                    } else {
                        header += buffer[i];
                    }
                    break;
                }

                case IN_SEQUENCE: {
                    ++seqSize;
                    if(isblank(buffer[i])){
                        num_ligne += (buffer[i] == '\n');
                    }else{
                        if(last_c == '\n' && buffer[i] == '+'){
                            fq.setSeq_size(seqSize - 2);
                            seqSize = 0;
                            state = IN_SEPARATOR;
                            separator += buffer[i];
                        } else {
                            std::cerr << "Erreur : " << the_file << " - Ligne " << num_ligne << std::endl;
                            std::cerr << "\t" << buffer[i] << " : pas un nucléotide" << std::endl;
                        }
                    }
                    break;
                }

                case IN_SEPARATOR: {
                    if (buffer[i] == '\n') {
                        ++num_ligne;
                        separator.clear();
                        state = IN_QUALITY;
                        fq.setQual_pos(beg + i + 1);
                    } else {
                        separator += buffer[i];
                    }
                    break;
                }
                case IN_QUALITY: {
                    qualSize++;
                    if (buffer[i] == '\n'){
                        ++num_ligne;
                        state = IN_TRANSITION;
                        fq.setQual_size(qualSize-1);
                        if (fq.getSeq_size() != qualSize -1){
                            std::cout << "Erreur : " << the_file << " - " << fq.getHeader() << std::endl;
                            std::cout << '\t' << "La qualité et la séquence ne correspondent pas" << std::endl;
                        }
                        qualSize = 0;
                        //std::cout << fq.getSeq(the_file) << std::endl;
                    }
                    break;
                }
                case IN_TRANSITION:{
                    if (last_c == '\n' && buffer[i] == '@'){
                        state = IN_HEADER;
                        header += buffer[i];
                    } else if (buffer[i] == '\n'){
                        ++num_ligne;
                    }
                }
            }
            last_c = buffer[i++];
        }

    }while(charNb);

    if (!header.empty()){
        fq.setHeader(header);
        header.clear();
    }

    fq.setSeq_size(seqSize);
    fq.setQual_size(qualSize);

    file.close();
}
