
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "InitialAntigen.h"
#include "MOLECULES.h"

using namespace std;

int main(int argc, char * argv[]) 
{
    InitialAntigen ia = InitialAntigen();
    string prefix = argv[1];
    string antigenName = prefix + ".txt";
    Molecule antigen = ia.load_antigen(antigenName);
    vector<int> epitopeResiduesNames = {164, 165, 166, 226, 227, 228, 229, 232, 310, 314};
    //vector<int> epitopeResiduesNames = {149, 150};
    //cout << antigen.residues[265].kind << endl;
    //cout << antigen.residues[266].kind << endl;
    //cout << antigen.residues[267].kind << endl;
    Molecule antigen_new = ia.move_antigen_epitope_znegative(antigen, epitopeResiduesNames);
    
    string outFile = "gp120_test.pdb";
    ia.output_molecule_pdb(antigen_new, outFile);
     
}

