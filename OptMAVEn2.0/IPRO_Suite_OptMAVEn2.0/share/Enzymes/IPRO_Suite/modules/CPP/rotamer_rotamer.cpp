/* Written in 2013 by Robert Pantazes of the Costas Maranas Lab in the Chemical
 * Engineering Department of the Pennsylvania State University
 
 * This file calculates the energies between rotamers */

// Include needed C++ files
# include <iostream>
# include <fstream>
# include <sstream>
# include <cstdlib>
# include <vector>
# include <string>
// Include the MOLECULES and ENERGY_FUNCTIONS header files
# include "MOLECULES.h"
# include "ENERGY_FUNCTIONS.h"

// Use the standard namespace
using namespace std;

// Start the program
int main(){
// Include the Residue and Atom classes
using MOLECULES::Residue; using MOLECULES::Atom;
// Include the IE_calculation function
using ENERGY_FUNCTIONS::IE_calculation;

// Load the rotamers. Prepare by opening the file and creating all needed
// variables
vector<Atom> residue; vector<Residue> rotamers; string line;
ifstream in1("ROTAMERS.txt");
// Read in the first line, get the Atom, and get its labelling information
getline(in1, line); Atom initial; initial.load(line);
// Do this differently depending on the file format
if (initial.fileFormat == "PDB"){
    int resNum = initial.residueNumber; int rotNum = initial.rotamerNumber;
    // Loop through the lines of the rotamers file
    while (!in1.eof()){
        // Create an Atom
        Atom atom; atom.load(line);
        // If this Atom belongs in a new rotamer
        if ((atom.residueNumber != resNum) || (atom.rotamerNumber != rotNum)){
            // Make a Residue of the Atoms
            Residue res; res.load(residue);
            // Store the rotamer
            rotamers.push_back(res); residue.clear();
            // Update the numbering information
            resNum = atom.residueNumber; rotNum = atom.rotamerNumber;}
        // Store the Atom and get the next line
        residue.push_back(atom); getline(in1, line);}
    // Store the last rotamer
    Residue res; res.load(residue); rotamers.push_back(res);
    // End the PDB section
    }
// If the file format is not supported
else {string temp = "The rotamer_rotamer.cpp program does not support the "
    + initial.fileFormat + " file format.";
    cout << temp << endl; exit(EXIT_FAILURE);}

// Calculate the rotamer-rotamer energies differently for each supported force
// field
if (initial.forceField == "CHARMM"){
    // Put the calculated results in this stream
    stringstream out;
    // Loop through the rotamers
    for(int i=0; i<rotamers.size()-1; i++){
        // Loop through them again
        for(int j=i+1; j<rotamers.size(); j++){
            // If these rotamers are at the same position, don't do the energy
            // calculation
            if (rotamers[j].number == rotamers[i].number){continue;}
            // CHARMM calculates three energy values. Store them in this array
            float energy [3] = {0, 0, 0};
            // Loop through the Atoms of the first rotamer
            for(int k=0; k<rotamers[i].atoms.size(); k++){
                Atom atom1 = rotamers[i].atoms[k];
                // Loop through the Atoms of the second rotamer
                for(int l=0; l<rotamers[j].atoms.size(); l++){
                    Atom atom2 = rotamers[j].atoms[l];
                    // Calculate the energies between the two Atoms
                    IE_calculation(atom1, atom2, energy, 0);
                    // End the for loops
                    }}
            // Put the calculated energy into the output stream
            float total = energy[0] + energy[1] + energy[2];
            out << rotamers[i].number << " " << rotamers[i].rotamerNumber << " "
                << rotamers[j].number << " " << rotamers[j].rotamerNumber << " "
                << total << "\n";
            // End the for loops searching through the rotamers
            }}
    // Write the results to a file
    ofstream outFile("RR_ENERGIES.txt");
    int count = 0;
    while (!out.eof()){
        if (count == 5){outFile << "\n"; count = 0;}
        count += 1; string text; out >> text; outFile << text; outFile << " ";}
    outFile.close();
    // End the CHARMM section of the program
    }
// If the force field is not supported
else {string text = "The rotamer_rotamer.cpp program does not support the "
    + initial.forceField + " force field.";
    cout << text << endl; exit(EXIT_FAILURE);}
// End the program
return 0;}
