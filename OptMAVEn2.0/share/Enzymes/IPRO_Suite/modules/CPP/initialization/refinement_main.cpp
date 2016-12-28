
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
    string prefix = argv[1];
    string antigenName = prefix + ".txt";
    InitialAntigen ia = InitialAntigen();
    Molecule antigen = ia.load_antigen(antigenName);
    string domainsName = "domains.txt";
    string summary = prefix + "_summary.txt";
    vector<string> partsBest = ia.antigen_position_refinement(antigen, domainsName, summary);
    
    string outFile = prefix +  "_refined.pdb";
    ia.output_molecule_pdb(antigen, outFile);
    string finalOptimizedParts = prefix + "_final_parts";
    ofstream out(finalOptimizedParts.c_str(), std::ios_base::out);
    
    for (vector<string>::iterator it = partsBest.begin(); it != partsBest.end(); it++)
    {
        out << *it << " ";
    }
    out << antigenName << endl;
    
    return 0; 
}

