
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "InitialAntigen.h"
#include "MOLECULES.h"

using namespace std;

int main(int argc, char * argv[]) 
{
    //Command line example: initialantigen.out 1ACY epitope.txt
    // A intance of InitialAntigen class
    InitialAntigen ia = InitialAntigen();
    // The pdb file prefix. 
    // For example 1ACY.pdb  prefix-->1ACY
    string prefix = argv[1];
    // Parameter file for the corresponding pdb file
    string antigenName = prefix + ".txt";
    // Read the epitope residue index from a text file
    // In the text file, each row stores a residue index
    string epitopeFile = argv[2];
    // Load the parameter file
    Molecule antigen = ia.load_antigen(antigenName);
    // Store the epitope residue index
    // Each index is substracted by 1 because of the residue array starting from 0
    vector<int> epitopeResiduesNames;
    ifstream in1(epitopeFile.c_str());
    string line;
    in1 >> line;
    while (!in1.eof())
    {
        int index = atof(line.c_str());
        epitopeResiduesNames.push_back(index - 1);
        in1 >> line;
    }
    
    //for(vector<int>::iterator it = epitopeResiduesNames.begin(); it != epitopeResiduesNames.end(); it++)
    //{
    //    cout << *it << endl;
    //}
    // Initial antigen positions
    // During this progress, the antigen is rotated to make the epitope having the most negative z coordinates
    // and head towareds the negative z direction. Then the antigen is translated along z direction (min->3.75, max->16.25, interval->1.25), x direction(min->-10, max->5, interval->2.5), y direction->(min->-5, max->10, interval->2.5) 
    // After that, the antigen is rotated along the z axis for 360 degree and 60 degree each time 
    // In total there are 11*7*7*6 = 3234 conformations
    vector<Molecule> antigen_news =ia.initial_antigen_position(antigen, epitopeResiduesNames);
    // Save each conformation to a pdb file
    int count = 1;
    for(vector<Molecule>::iterator it = antigen_news.begin(); it != antigen_news.end(); it++)
    {
        string name = to_string(count);
        string outFile = "antigen" + name + ".pdb";
        ia.output_molecule_pdb(*it, outFile);
        count += 1;
    }
     
}

