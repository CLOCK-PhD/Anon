#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <string>

using namespace std;

struct snp_kmer{
    string kmer;
    string id;
    string chromosome;
    size_t snp_pos;
    size_t kmer_pos;
};

int main(){
    
    fstream myFile;
    
    string line;
    //char buffer[1024];
    size_t charNb = 0;
    
    //myFile.open("../download/example.fas", ios::in); //read
    myFile.open("../download/grch38p7_prim_assembly.fasta", ios::in); //read
    if(myFile.is_open()){
        while(getline(myFile, line)){
            cout << line << endl;
        }
        myFile.close();
    }
    else cout << "Unable to open file";

    
    /*ifstream is("../download/grch38p7_prim_assembly.fasta", ifstream::binary);
    if(is){
        // get length of file
        is.seekg(0, is.end);
        int length = is.tellg();
        is.seekg(0, is.beg);

        // allocate memory
        char * buffer = new char [length];

        // read data as a block:
        is.read(buffer, length);

        is.close();

        // print content:
        cout.write(buffer, length);

        delete[] buffer;
    }*/

    return 0;

}