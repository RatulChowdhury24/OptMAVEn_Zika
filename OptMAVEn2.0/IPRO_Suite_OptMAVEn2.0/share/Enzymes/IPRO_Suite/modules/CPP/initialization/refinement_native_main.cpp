
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "InitialAntigen.h"
#include "MOLECULES.h"

using namespace std;

int main(int argc, char * argv[]) 
{
    // refinement.out   
    string antigenPrefix = argv[1];
    string antigenName = antigenPrefix + ".txt";
    string antibodyHPrefix = argv[2];
    string antibodyHName = antibodyHPrefix + ".txt";
    string antibodyLPrefix = argv[3];
    string antibodyLName = antibodyLPrefix + ".txt";
    InitialAntigen ia = InitialAntigen();
    Molecule antigen = ia.load_antigen(antigenName);
    Molecule antibodyH = ia.load_antigen(antibodyHName);
    Molecule antibodyL = ia.load_antigen(antibodyLName);
    
    Molecule antigen_refined = ia.antigen_position_refinement(antigen, antibodyH, antibodyL);
    
    string outFile = antigenPrefix +  "_refined.pdb";
    ia.output_molecule_pdb(antigen_refined, outFile);
    
    return 0; 
}

