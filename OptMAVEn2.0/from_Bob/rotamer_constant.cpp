/* Written in 2013 by Robert Pantazes of the Costas Maranas Lab in the Chemical
 * Engineering Department of the Pennsylvania State University
 
 * This file calculate the energies between rotamers and the constant portions
 * of a DesignGroup. */

// Include needed C++ files
# include <iostream>
# include <fstream>
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
// Include the DesignGroup and Residue classes
using MOLECULES::DesignGroup;
using MOLECULES::Residue;
using MOLECULES::Atom;
// Include the exclusions and IE_calculation functions
using ENERGY_FUNCTIONS::exclusions;
using ENERGY_FUNCTIONS::IE_calculation;

// Open the file that contains the constant portions of the system
ifstream in1("CONSTANT.txt");
// Create the variables needed to read in the contents of the file
vector<Atom> atoms; string line;
// Read in the file
getline(in1, line);
while(!in1.eof()){
    // Store the line as an Atom
    Atom atom; atom.load(line); atoms.push_back(atom);
    // Read in the next line
    getline(in1, line);}
// Close the file
in1.close();
// Create the DesignGroup
DesignGroup group; group.load(atoms);

// Read in the rotamers and store them in this vector
ifstream in2("ROTAMERS.txt"); vector<Residue> rotamers;
// Clear the atoms vector
atoms.clear();
// Create the first Atom in the file
getline(in2, line); Atom temp; temp.load(line);
// Get its labelling information. Do the rest of this differently depending on
// the file format of the Atoms
if (temp.fileFormat == "PDB"){
    int resNum = temp.residueNumber; int rotNum = temp.rotamerNumber;
    // Read in the rotamers
    while(!in2.eof()){
        // Create an Atom out of the line
        Atom atom; atom.load(line);
        // If this Atom is the start of a new rotamer
        if ((atom.residueNumber != resNum) || (atom.rotamerNumber != rotNum)){
            // Make a Residue for the rotamer
            Residue res; res.load(atoms); 
            // Store the rotamer and clear the atoms vector
            rotamers.push_back(res); atoms.clear();
            // Update the numbering information
            resNum = atom.residueNumber; rotNum = atom.rotamerNumber;}
        // Store the Atom and get the next line in the file
        atoms.push_back(atom); getline(in2, line);}
    // Close the input file
    in2.close();
    // Store the last rotamer
    Residue res; res.load(atoms); rotamers.push_back(res);}
// If the file format isn't supported
else {string text = "The rotamer_constant.cpp program does not support the "
    + temp.fileFormat + " file format.";
    cout << text << endl; exit(EXIT_FAILURE);}

// Handle the remainder of this file differently depending on the force field
// used by the Atoms
if (group.forceField == "CHARMM"){
    // Create an array to store an energy for each rotamer
    float energies [rotamers.size()];
    // Loop through the rotamers
    for(int I=0; I<rotamers.size(); I++){
        // CHARMM calculates three energies: VDW, electrostatics, and solvation.
        // Create an array to store those energies for this rotamer
        float energy [3] = {0, 0, 0};
        // Loop through the rotamer's Atoms
        for(int J=0; J<rotamers[I].atoms.size(); J++){
            // Extract the Atom
            Atom atom1 = rotamers[I].atoms[J];
            // Now loop through the DesignGroup
            for(int i=0; i<group.molecules.size(); i++){
            for(int j=0; j<group.molecules[i].residues.size(); j++){
            for(int k=0; k<group.molecules[i].residues[j].atoms.size(); k++){
                // Extract the Atom
                Atom atom2 = group.molecules[i].residues[j].atoms[k];
                // Determine which energy terms should be calculated
                int skip = exclusions(atom1, atom2);
                // Calculate the energy between these two Atoms
                IE_calculation(atom1, atom2, energy, skip);
                // End all 4 for loops
                }}}}
        // Store the total energy of the rotamer
        energies[I] = energy[0] + energy[1] + energy[2];
        // End the rotamers for loop
        }
    // Now create an output of the calculated energies
    ofstream out1("RC_ENERGIES.txt");
    for(int i=0; i<rotamers.size(); i++){
        out1 << rotamers[i].number << " " << rotamers[i].rotamerNumber << " "
             << energies[i] << endl;}
    out1.close();
    // End the CHARMM force field calculation section
    }
// If the force field is not supported
else {string text = "The rotamer_constant.cpp program does not support the "
    + group.forceField + " force field.";
    cout << text << endl; exit(EXIT_FAILURE);}

// End the program
return 0;}
