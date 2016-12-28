
#ifndef DEIMMUNIZATION_H
#define DEIMMUNIZATION_H


#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <utility>


using namespace std;


class Deimmunization {
    
    public:
        
        //Constructor
        Deimmunization();
        
        //Deconstructor
        virtual ~Deimmunization();
        
        //// Load the human 9mer sequence database
        vector<string> load_humseq_database(string & fileName);
        
        //Given a sequence, identify all the possible 9mers sequences
        vector<string> identify_9mers_sequences(string & sequence);
        
        // Reads sequence from pdb and returns in single amino acid code as a string
        string extract_pdbseq(string & pdbName);
        
        // Create a amino acid three-letter and one-letter name map
        void create_aamap();
        
        // Create 20 amino acids one-letter name arrays
        void create_aaArray();
         
        // Count the mutations between two seqs
        int count_mutations(string & seq1, string & seq2);
        
        //Count the minimum mutations between a target seq and a seq database
        int count_minimum_mutations(string & seq, vector<string> & database);
        
        //Count the total mutations given a series of sequences and human sequence database
        int count_total_mutations(vector<string> & sequences, vector<string> & database); 
        
        // Determine the begin position
        int determine_begin_position(int & position, int & gap);
        
        // Determine the end position
        int determine_end_position(int & length, int & position, int & gap);

        //Given a sequence, begin and end position, gap, extract all the possible (gap+1)-mer sequences 
        vector<string> extract_all_sequences(string & seq, int & begin, int & end, int & gap);
        // Make the integer cuts
        void make_humanization_cuts(string & sequence, vector<int> & positions);
    
    private:

         map<string, string>  aaMap; 
         string aaArray [20];
};

#endif
