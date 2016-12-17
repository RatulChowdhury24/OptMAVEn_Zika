
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
    InitialAntigen ia = InitialAntigen();
    string antigenName = antigenPrefix + ".txt";
    Molecule antigen = ia.load_antigen(antigenName);
    string parts [] = {"HV.txt", "LV.txt", "HCDR3.txt", "LCDR3.txt", "HJ.txt", "LJ.txt"};
    float totalEnergy = 0.0;
    //for (vector<string>::iterator it= parts.begin(); it != parts.end(); it++)
    for (int n = 0; n < 6; ++n)
    {
        Molecule part = ia.load_antigen(parts[n]);
        bool solvation = false;
        float energy = ia.calculate_interaction_energy(part, antigen, solvation);
        totalEnergy += energy;
    }

    string energyFile = antigenPrefix + "_xray_cpp_energy.txt";
    ofstream of;
    of.open(energyFile.c_str(), ios::out);
    of << totalEnergy << endl;
    
    return 0; 
}

