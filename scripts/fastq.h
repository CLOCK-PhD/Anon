#ifndef __FASTQ_H__
#define __FASTQ_H__

#include <string>
#include <vector>

class Fastq{
    protected :
        std::string header;
        std::string type;
        size_t seq_pos;
        size_t seq_size;

    private :
        size_t qual_pos;
        size_t qual_size;

    public:
        // CONSTRUCTEUR
        Fastq();
        Fastq(std::string header, size_t seq_pos, size_t seq_size, size_t qual_pos, size_t qual_size);

        // GETTERS
        std::string getHeader() const;
        size_t getSeq_pos() const;
        size_t getSeq_size() const;
        size_t getQual_pos() const;
        size_t getQual_size() const;
        // A voir
        std::string getSeq(const std::string &the_file);        // Récupérer la séquence
        std::string getQual(const std::string &the_file);       // Récupérer la qualité

        // SETTERS
        void setHeader(std::string header);
        void setSeq_pos(size_t seq_pos);
        void setSeq_size(size_t seq_size);
        void setQual_pos(size_t qual_pos);
        void setQual_size(size_t qual_size);

        // Autres
        static bool isBlank(char c);
        static bool isNucl(char c, bool degeneratted);

        // Test du coup
        static void lecture(std::string &the_file);
        //static void lecture(std::ifstream the_file);

        enum FileReadState {
            IN_HEADER,
            IN_SEQUENCE,
            IN_SEPARATOR,
            IN_QUALITY,
            IN_TRANSITION
        };

};

#endif