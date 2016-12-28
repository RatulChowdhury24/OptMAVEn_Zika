
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include "InitialAntigen.h"
#include "MOLECULES.h"
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>

using namespace std;

int main(int argc, char * argv[]) 
{
    // Given the name of the conformations file  "conformations.txt"
    // which includes the antigen name prefix such antigen1
    string conformations = argv[1];
    string domainsName = "domains.txt";
    InitialAntigen ia = InitialAntigen();
    map<string, Molecule> parts = ia.load_parts(domainsName);
    vector<string> antigenNames;
    //cout << "1" << endl;
    ifstream in1(conformations.c_str());
    string line;
    in1 >> line;
    while (!in1.eof())
    {
        boost::trim(line);
        antigenNames.push_back(line);
        in1 >> line;
    }
    in1.close();
    //cout << "2" << endl;
    for(vector<string>::iterator it = antigenNames.begin(); it != antigenNames.end(); it++)
    {
        string antigenName = *it + ".txt";
        string energyFile = *it + "_maps_energies.txt";
        ifstream ifile(energyFile.c_str());
        if(!ifile) 
        { 
            Molecule antigen = ia.load_antigen(antigenName);
            if (!ifile)
            {
                ifile.close();
                bool solvation = false;
                ia.output_maps_energies(parts, antigen, energyFile, solvation);
            }
        }
    }
    return 0; 
}

