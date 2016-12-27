/* Written in 2013 by Robert Pantazes of the Costas Maranas Lab in the Chemical
 * Engineering Department of the Pennsylvania State University
 
 * This file contains functions for calculating the energies between MAPs
 * database structures and an antigen */

// Include needed C++ files
# include <iostream>
# include <sstream>
# include <vector>
# include <string>

// Include IPRO Suite C++ header information
# include "MOLECULES.h"
# include "ENERGY_FUNCTIONS.h"

// Use the standard name space
using namespace std;

// Start the program
int main(){
// Include the Molecule and Atom classes
using MOLECULES::Molecule;
using MOLECULES::Atom;
// Include the ability to do energy calculations
using ENERGY_FUNCTIONS::IE_calculation;

// Load the antigen
ifstream in1("Antigen.txt");
// Be prepared to store its information
vector<Atom> atoms; string line;
getline(in1, line);
// Go through every Atom in the antigen
while(!in1.eof()){
    // Store the Atom
    Atom atom; atom.load(line); atoms.push_back(atom);
    // Get the next line
    getline(in1, line);}
in1.close();
// Make a Molecule
Molecule antigen; antigen.load(atoms);

// Store the calculated energies in this string stream
stringstream results;
// Read in the types of domains being considered
vector<string> kinds;
ifstream in2("domains.txt");
in2 >> line;
while(!in2.eof()){kinds.push_back(line); in2 >> line;}
in2.close();
// The different maps parts
string regions [3] = {"V", "CDR3", "J"};
// Go through the different region / part combinations
for(int i=0; i<kinds.size(); i++){
    for(int j=0; j<3; j++){
        // Make the name of the expected part and file
        string partName = kinds[i] + regions[j];
        string fileName = partName + ".txt";
        // Open the expected file
        ifstream in3(fileName.c_str());
        // Store the parts in this vector
        vector<Molecule> parts;
        // Clear the Atoms vector
        atoms.clear();
        // Get the information about the first Atom in the file
        getline(in3, line);
        Atom a; a.load(line);
        int oldN = a.rotamerNumber;
        // Read in the Atoms in the file
        while(!in3.eof()){
            Atom atom; atom.load(line);
            // Determine if this Atom is the first in a new part or not
            if(atom.rotamerNumber != oldN){
                // Store the part
                Molecule part; part.load(atoms); parts.push_back(part);
                // Clear the Atom list
                atoms.clear();
                // Store the new part number
                oldN = atom.rotamerNumber;}
            // Store the atom
            atoms.push_back(atom);
            // Get the next line
            getline(in3, line);}
        // Close the file
        in3.close();
        // Store the last part
        Molecule last; last.load(atoms); parts.push_back(last);
        // Calculate the energy between each part and the antigen
        for(int A=0; A<parts.size(); A++){
            // The energies between this structure and the antigen
            float energies [3] = {0, 0, 0};
            // Go through the Residues and Atoms of the part
            for(int a=0; a<parts[A].residues.size(); a++){
                for(int b = 0; b<parts[A].residues[a].atoms.size(); b++){
                    Atom atom1 = parts[A].residues[a].atoms[b];
                    // Go through the Residues and Atoms of the antigen
                    for(int Z=0; Z<antigen.residues.size(); Z++){
                        for(int z=0; z<antigen.residues[Z].atoms.size(); z++){
                            Atom atom2 = antigen.residues[Z].atoms[z];
                            IE_calculation(atom1, atom2, energies, 0);}}}}
            // Calculate the total energy for this MAPs structure
            float energy = energies[0] + energies[1] + energies[2];
            // Store this information in the results stream
            results << partName << " " << A+1 << " " << energy << "\n";
            }}}

// Write the results to a file
ofstream out1("MAPs_Energies.txt");
int count = 0;
while(!results.eof()){
    // Make sure end line characters are included
    if(count == 3){out1 << "\n"; count = 0;}
    // output the information for each result
    string text; results >> text; out1 << text; count++;}
out1.close();

// End the program
return 0;}
