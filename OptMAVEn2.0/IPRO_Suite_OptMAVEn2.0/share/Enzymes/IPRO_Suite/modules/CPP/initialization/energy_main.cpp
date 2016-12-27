
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "InitialAntigen.h"
#include "MOLECULES.h"

using namespace std;

int main(int argc, char * argv[]) 
{
    string antigenPrefix = argv[1];
    //cout << antigenPrefix << endl;
    string domainsName = "domains.txt";
    InitialAntigen ia = InitialAntigen();
    string antigenName = antigenPrefix + ".txt";
    //cout << antigenName << endl;
    Molecule antigen = ia.load_antigen(antigenName);
    map<string, Molecule> parts = ia.load_parts(domainsName);
    string energyFile = antigenPrefix + "_maps_energies.txt";
    bool solvation = false;
    ia.output_maps_energies(parts, antigen, energyFile, solvation);
    return 0; 
}

