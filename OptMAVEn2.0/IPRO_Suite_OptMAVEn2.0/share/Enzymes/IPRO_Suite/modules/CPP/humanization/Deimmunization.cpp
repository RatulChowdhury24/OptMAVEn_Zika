
#include <fstream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <initializer_list>

#include "Deimmunization.h"
//Constructor
Deimmunization::Deimmunization() 
{
   Deimmunization::create_aamap(); 
   Deimmunization::create_aaArray(); 
}

//Deconstructor
Deimmunization::~Deimmunization() {}

// Create a amino acid three-letter and one-letter name map
void Deimmunization::create_aamap() 
{
    aaMap["ALA"] = "A";
    aaMap["A"] = "ALA";
    aaMap["CYS"] = "C";
    aaMap["C"] = "CYS";
    aaMap["ASP"] = "D";
    aaMap["D"] = "ASP";
    aaMap["GLU"] = "E";
    aaMap["E"] = "GLU";
    aaMap["PHE"] = "F";
    aaMap["F"] = "PHE";
    aaMap["GLY"] = "G";
    aaMap["G"] = "GLY";
    aaMap["HIS"] = "H";
    aaMap["H"] = "HIS";
    aaMap["HSD"] = "H";
    aaMap["ILE"] = "I";
    aaMap["I"] = "ILE";
    aaMap["LYS"] = "K";
    aaMap["K"] = "LYS";
    aaMap["LEU"] = "L";
    aaMap["L"] = "LEU";
    aaMap["MET"] = "M";
    aaMap["M"] = "MET";
    aaMap["ASN"] = "N";
    aaMap["N"] = "ASN";
    aaMap["PRO"] = "P";
    aaMap["P"] = "PRO";
    aaMap["GLN"] = "Q";
    aaMap["Q"] = "GLN";
    aaMap["ARG"] = "R";
    aaMap["R"] = "ARG";
    aaMap["SER"] = "S";
    aaMap["S"] = "SER";
    aaMap["THR"] = "T";
    aaMap["T"] = "THR";
    aaMap["VAL"] = "V";
    aaMap["V"] = "VAL";
    aaMap["TRP"] = "W";
    aaMap["W"] = "TRP";
    aaMap["TYR"] = "Y";
    aaMap["Y"] = "TYR";
}

void Deimmunization::create_aaArray()
{
    aaArray[0] = "A";
    aaArray[1] = "R";
    aaArray[2] = "N";
    aaArray[3] = "D";
    aaArray[4] = "C";
    aaArray[5] = "Q";
    aaArray[6] = "E";
    aaArray[7] = "G";
    aaArray[8] = "H";
    aaArray[9] = "I";
    aaArray[10] = "L";
    aaArray[11] = "K";
    aaArray[12] = "M";
    aaArray[13] = "F";
    aaArray[14] = "P";
    aaArray[15] = "S";
    aaArray[16] = "T";
    aaArray[17] = "W";
    aaArray[18] = "Y";
    aaArray[19] = "V";
}
// Load the human 9mer sequence database
vector<string> Deimmunization::load_humseq_database(string & fileName) 
{

    vector<string> seqList;
    ifstream infile;
    infile.open(fileName.c_str());
    string line;
    while(!infile.eof()) 
    {
        infile >> line;
        seqList.push_back(line);
    }

    infile.close();

    return seqList;
        
}

//Given a sequence, identify all the possible 9mers sequences
vector<string> Deimmunization::identify_9mers_sequences(string & seq) 
{
    
    const int minLength = 9;
    int seqLength = seq.length();
    vector<string> seqList;
    if (seqLength < minLength) 
    {
        cout << "Sequence length is less than 9 and no enough amino acids to identify" << endl;
        exit(EXIT_FAILURE);
    }

    else 
    {
        int gap = seqLength - minLength;
        for (int i = 0; i <= gap ; i++) 
        {
            string str = "";
            for (int j = i ; j < (i + minLength) ; j++)
            {
                str += seq[j];
            }
            seqList.push_back(str);

        }

    }

    return seqList;
}


// Read sequence from pdb and return amino acid sequence as a string
string Deimmunization::extract_pdbseq(string & pdbName) 
{

    string sequence;
    ifstream infile;
    infile.open(pdbName.c_str());
    
    string line;
    string prefix = "ATOM";
    while (!infile.eof())
    {
        getline(infile, line);
        //cout << line << endl;
        //cout << line[13] << " " << line[14] << " " << line[17] << endl;
        if (line.find(prefix) == 0 && line[13] == 'C' && line[14] == 'A')
        {
            
            string aa = "";
            aa += line[17];
            aa += line[18];
            aa += line[19];
            sequence += aaMap[aa];
            
        }
    }

    return sequence;
}



// Count the mutations between two seqs amd the two seqs have the same the length of amino acids
int Deimmunization::count_mutations(string & seq1, string & seq2)
{
    int len1 = seq1.length();
    int len2 = seq2.length();
    int mutations = 0;
    // If the two seqs have no any AA or the same number of AA, exit
    if (len1 == 0 && len1 != len2) 
    {
        cout << "The two sequences both have no any amino acid or have amino acids with different lengths" << endl;
        exit(EXIT_FAILURE);
    }
    else 
    {
        for (int i = 0; i < len1 ; i++)
        {
            if ( seq1[i] != seq2[i] ) 
            {
                mutations++;
            }
        }
    }

    return mutations;

}

//Count the minimum mutations between a target seq and a seq database
// Return a pair with the first element is a vector including all the seqs with the same minimu tations
// The second element is the minimum mutations
int Deimmunization::count_minimum_mutations(string & seq, vector<string> & database)
{
    // Define a variable to store the minumum muations
    int minMutations = 10000;
    // Store each sequence in the database and mutations between this seq and the target seq
    //map<string, int> seqMutMap;
    // All the sequences in the database having the minimum mutations with the target seq
    // Int is the minMutations
    //pair<vector<string>, int> seqMutPair ;
    // All the sequences in the database having the minimum mutations
    //vector<string> miniSeqs;
    //Calulate the minimum mutations between a target seq and the seq database
    for(vector<string>::iterator it = database.begin(); it != database.end(); ++it)
    {
        int mutations = count_mutations(seq, *it);
        //cout << mutations << endl;
        //seqMutMap[*it] = mutations;

        if (mutations < minMutations) 
        {
            minMutations = mutations;
        }
    }
    
    return minMutations;
    // Store all the sequences with the same minimum mutations with the target seq
    //for(map<string, int>::iterator it = seqMutMap.begin(); it != seqMutMap.end(); ++it)
    //{
    //    if (it->second == minMutations)
    //    {
    //        miniSeqs.push_back(it->first);
    //    }

    //}
    // Make the seqs and minMutations pair
    //seqMutPair = make_pair(miniSeqs, minMutations);
    
    //return seqMutPair;

}
            
//Count the total mutations given a series of sequences and human sequence database
int Deimmunization::count_total_mutations(vector<string> & sequences, vector<string> & database)
{
    int total = 0;
    for(vector<string>::iterator it = sequences.begin(); it != sequences.end(); ++it)
    {
        
        //pair<vector<string>, int> seqMutPair = count_minimum_mutations(*it, database);
        int minMutations = count_minimum_mutations(*it, database);
        total += minMutations;
        //total += seqMutPair.second;
    }

    return total;
}

//Given a position, gap, determine the begin position for extracting all the possible (gap+1)-mer sequence

int Deimmunization::determine_begin_position(int & position, int & gap)
{
    int begin = (position - gap) > 0 ? position - gap : 0;
    
    return begin;
}

//Given a sequence length, position, gap, determine the end position for extracting all the possible (gap+1)-mer sequence
int Deimmunization::determine_end_position(int & length, int & position, int & gap)
{
    int end = (position + gap) < length ? position + gap : length;
    
    return end;

}

//Given a sequence, begin and end position, gap, extract all the possible (gap+1)-mer sequences 
vector<string> Deimmunization::extract_all_sequences(string & seq, int & begin, int & end, int & gap)
{
   vector<string> sequences; 

   for (int index = begin; index <= end - gap; index++)
   {
        sequences.push_back(seq.substr(index, index + gap + 1));
   }

   return sequences;
}

//Given a protein amino acid sequence and the perturbation positions, return the cuts
void Deimmunization::make_humanization_cuts(string & sequence, vector<int> & positions)
{
    //cout << "sequence: " << sequence << endl;
    int lenSequence = sequence.size();
    int lenPositions = positions.size(); 
    if ((lenPositions > 3) || (lenPositions == 0))
    {
        cout << "The perturbation positions are not correct, please check" << endl;
        exit(EXIT_FAILURE);
    }
    //cout << lenSequence << " " << lenPositions << endl;
    string humanSeqs = "Human_9mer_Sequences.txt";
    vector<string> database = load_humseq_database(humanSeqs);
    int gap = 8;    
    vector< map<int, string> > cuts; 
    
    if (lenPositions == 1)
    {
        int p = positions[0];
        int begin = determine_begin_position(p, gap);
        int end = determine_end_position(lenSequence, p, gap);
        //cout << "begin " << begin << endl;
        //cout << "end " << end << endl;
        vector<string> seqs = extract_all_sequences(sequence, begin, end, gap);
        int count = count_total_mutations(seqs, database);
        //cout << "count " << count << endl;
        int iteration = 0;
        for(int index = 0; index < 20 ; index++)
        {
            string sequence_mut(sequence);
            sequence_mut.replace(p, 1, aaArray[index]);
            //cout << "sequence " << sequence << endl; 
            //cout << "sequence_mut " << sequence_mut << endl; 
            vector<string> seqs_mut = extract_all_sequences(sequence_mut, begin, end, gap);
            int count_mut = count_total_mutations(seqs_mut, database);
            //cout << "count_mut " << count_mut << endl;
            iteration++;
            //cout << iteration << endl;
            if (count_mut > count)
            {
                map<int, string> solutions;
                solutions[p] = aaArray[index];
                cuts.push_back(solutions);
            }

        }

    }

    if (lenPositions == 2)
    {
        int p1 = positions[0];
        int p2 = positions[1];
        int begin = determine_begin_position(p1, gap);
        int end = determine_end_position(lenSequence, p2, gap);
        //cout << "begin " << begin << endl;
        //cout << "end " << end << endl;
        vector<string> seqs = extract_all_sequences(sequence, begin, end, gap);
        int count = count_total_mutations(seqs, database);
        //cout << "count " << count << endl;
        int iteration = 0;
        for(int index1 = 0; index1 < 20 ; index1++)
        {
            string sequence_mut1(sequence);
            sequence_mut1.replace(p1, 1, aaArray[index1]);
            for(int index2 = 0; index2 < 20 ; index2++)
            {
                string sequence_mut2(sequence_mut1);
                sequence_mut2.replace(p2, 1, aaArray[index2]);
                //cout << "sequence " << sequence << endl; 
                //cout << "sequenc2 " << sequence_mut2 << endl; 
                
                vector<string> seqs_mut2 = extract_all_sequences(sequence_mut2, begin, end, gap);
                int count_mut = count_total_mutations(seqs_mut2, database);
                //cout << "count_mut " << count_mut << endl;
                iteration++;
                //cout << iteration << endl;
                if (count_mut > count)
                {
                    map<int, string> solutions;
                    solutions[p1] = aaArray[index1];
                    solutions[p2] = aaArray[index2];
                    cuts.push_back(solutions);
                }
            }
        }
    }
    
    if (lenPositions == 3)
    {
        int p1 = positions[0];
        int p2 = positions[1];
        int p3 = positions[2];
        int begin = determine_begin_position(p1, gap);
        int end = determine_end_position(lenSequence, p3, gap);
        //cout << "begin: " << begin << " " << "end: " << end << endl;
        vector<string> seqs = extract_all_sequences(sequence, begin, end, gap);
        int count = count_total_mutations(seqs, database);
        //cout << "count: " << count << endl;
        int iteration = 0;
        for(int index1 = 0; index1 < 20 ; index1++)
        {
            string sequence_mut1(sequence);
            sequence_mut1.replace(p1, 1, aaArray[index1]);
            for(int index2 = 0; index2 < 20 ; index2++)
            {
                string sequence_mut2(sequence_mut1);
                sequence_mut2.replace(p2, 1, aaArray[index2]);
                for(int index3 = 0; index3 < 20 ; index3++)
                {
                    string sequence_mut3(sequence_mut2);
                    sequence_mut3.replace(p3, 1, aaArray[index3]);
                    //cout << "sequence " << sequence << endl; 
                    //cout << "sequence_mut3 " << sequence_mut3 << endl; 
                    vector<string> seqs_mut3 = extract_all_sequences(sequence_mut3, begin, end, gap);
                    int count_mut = count_total_mutations(seqs_mut3, database);
                    //cout << "count_mut " << count_mut << endl;
                    iteration++;
                    //cout << iteration << endl;
                    if (count_mut > count)
                    {
                        map<int, string> solutions;
                        solutions[p1] = aaArray[index1];
                        solutions[p2] = aaArray[index2];
                        solutions[p3] = aaArray[index3];
                        cuts.push_back(solutions);
                    }
                }
            }
        }
    }

    
    ofstream out;
    out.open("cuts.txt", ios::out);
    //cout << cuts.size() << endl;
    for(vector< map<int, string> >::iterator it = cuts.begin(); it != cuts.end(); it++)
    {
        //cout << (*it).first << " " << (*it).second << endl; 
        for(map<int, string>::iterator mt = (*it).begin(); mt != (*it).end(); mt++)
        {
            out << (*mt).first << ":" << aaMap[(*mt).second] << " ; " ;
        }
        out << endl;
    }
    
    out.close();

}

















