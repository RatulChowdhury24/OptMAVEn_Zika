/* Written in 2013 by Robert Pantazes of the Costas Maranas Lab in the Chemical
 * Engineering Department of the Pennsylvania State University
 
 * This file contains functions that are used during docking. There's no reason
 * they couldn't go into the docking.cpp code, other than that it would make
 * that file much longer. So I put them here instead. */

// Have some preprocessor commands so the file is only loaded once
#ifndef DOCKING_FUNCTIONS_CHECK
#define DOCKING_FUNCTIONS_CHECK 1
// Include C++ modules
# include <cmath>
# include <string>
# include <vector>
# include <cstdlib>
# include <iostream>
# include <sstream>
# include <fstream>
// Include the Molecules and Energy functions
# include "MOLECULES.h"
# include "ENERGY_FUNCTIONS.h"

// Use the standard name space
using namespace std;

// Define a name space for docking functions
namespace DOCKING_FUNCTIONS
{

void calculate_rmatrix(float angle, char direction, float rmatrix [][3]){
    /* This function is analogous to the function in the MOLECULES.py module
     * that calculates a rotation matrix using Rodriguez's rotation formula.
     * However, here we account for the fact that we're only going to be
     * rotating around an axis */
    // Calculate some frequently used values
    float c = cos(angle);
    float s = sin(angle);
    float v = 1 - c;
    // Get the proper unit vector
    float r [3] = {0, 0, 0};
    if (direction == 'x'){r [0] = 1.0;}
    else if(direction == 'y'){r[1] = 1.0;}
    else if(direction == 'z'){r[2] = 1.0;}
    // Calculate the entries of the matrix
    rmatrix[0][0] = c + v*r[0]*r[0];
    rmatrix[0][1] = -s*r[2] + v*r[0]*r[1];
    rmatrix[0][2] = s*r[1] + v*r[0]*r[2];
    rmatrix[1][0] = s*r[2] + v*r[1]*r[0];
    rmatrix[1][1] = c + v*r[1]*r[1];
    rmatrix[1][2] = -s*r[0] + v*r[1]*r[2];
    rmatrix[2][0] = -s*r[1] + v*r[2]*r[0];
    rmatrix[2][1] = s*r[0] + v*r[2]*r[1];
    rmatrix[2][2] = c + v*r[2]*r[2];
    }

int load_docking_data(vector<vector<MOLECULES::Molecule> > &Docking, 
                       vector<vector<float> > &parameters){
    /* This function loads the molecules that are eligible to MOVE during
     * docking and stores them in the Docking vector */
    // Open the input file
    ifstream in1("docking_information.txt");
    // Read in the number of docking groups
    int GN; in1 >> GN;
    // Load each docking group
    for(int i=0; i<GN; i++){
        // Read in the number of Molecules in this Docking Group
        int groupSize; in1 >> groupSize;
        // A vector to store the Molecules of this Docking Group
        vector<MOLECULES::Molecule> molecules;
        // Load each Molecule
        for(int j=0; j<groupSize; j++){
            // The Molecule's file's name and Design Group number (stored as a
            // string)
            string fileName, groupNumber;
            in1 >> fileName >> groupNumber;
            // Create a Molecule class object
            MOLECULES::Molecule molecule;
            // Store it's name and group number
            molecule.fileName = fileName;
            molecule.groupNumber = groupNumber;
            // Make the name useable for opening a file
            char FN [molecule.fileName.size()];
            for(int k=0; k<molecule.fileName.size(); k++){
                FN[k] = molecule.fileName[k];}
            // open that file
            ifstream in2(FN);
            // Create the variables to read this in
            string line; vector<MOLECULES::Atom> atoms;
            getline(in2, line);
            while(!in2.eof()){
                // Make an Atom of the current line
                MOLECULES::Atom atom; atom.load(line);
                // Store the Atom
                atoms.push_back(atom);
                // Get the next line
                getline(in2, line);}
            // Close the file
            in2.close();
            // Load the Molecule
            molecule.load(atoms);
            // Store the Molecule
            molecules.push_back(molecule);}
        // Store the Docking Group
        Docking.push_back(molecules);}
    // Read in the number of docking iterations
    int iterations; in1 >> iterations;
    // Read in the parameters for docking
    for(int i=0; i<iterations; i++){
        // There's one set of random movements for each Docking Group for each
        // iteration
        for(int j=0; j<Docking.size(); j++){
            // Read in the X rotation, Y rotation, Z rotation, X perturbation, Y
            // perturbation, Z perturbation, simulated annealing 'temperature'
            // (actually T * GasConstant), and the random number for docking
            vector<float> data;
            for(int k=0; k<8; k++){
                float parameter; in1 >> parameter; data.push_back(parameter);}
            // Store the parameters
            parameters.push_back(data);}}
    // Close the input file
    in1.close();
    return iterations;}

void load_constant_molecules(vector<vector<MOLECULES::Molecule> > &Docking, 
                            vector<MOLECULES::DesignGroup> &constants){
    /* Load all of the Molecules that don't move during docking but are needed
     * for energy calculations to determine whether or not to keep the results
     * of docking */
    // Figure out what those Design Groups are
    vector<string> groupNumbers;
    // Loop through the Molecules
    for(int i=0; i<Docking.size(); i++){for(int j=0; j<Docking[i].size(); j++){
        // Use a boolean value to determine if this Molecule's group number has
        // been stored yet
        bool need = true;
        for(int k=0; k<groupNumbers.size(); k++){
            if(groupNumbers[k] == Docking[i][j].groupNumber){
                need = false; break;}}
        // If this is a new group number, store it
        if(need){groupNumbers.push_back(Docking[i][j].groupNumber);}}}
    // Loop through the Design Groups
    for(int i=0; i<groupNumbers.size(); i++){
        // Make the name of the file
        string fileName = "constant" + groupNumbers[i] + ".txt";
        // Make it into something useable for opening a file
        char FN[fileName.size()];
        for(int j=0; j<=fileName.size(); j++){FN[j] = fileName[j];}
        // Open that file
        ifstream in1(FN);
        // Read in the Atoms
        string line; vector<MOLECULES::Atom> atoms;
        getline(in1, line);
        while(!in1.eof()){
            // Make an Atom out of that line
            MOLECULES::Atom atom; atom.load(line);
            // Store the Atom
            atoms.push_back(atom);
            // Get the next line
            getline(in1, line);}
        // Close the file
        in1.close();
        // Create a Design Group
        MOLECULES::DesignGroup group; group.load(atoms, groupNumbers[i]);
        // Store that Design Group
        constants.push_back(group);}
    // And that ends the function
    }

void move_docking_group(vector<MOLECULES::Molecule> &molecules, 
                        vector<float> &data){
    /* Carry out a set of random motions to a group of Molecules that all need
     * to move together. */
    // First, we have to center the group of Molecules. Begin by calculating the
    // average coordinates of the atoms
    float coors [3] = {0, 0, 0};
    // Loop through every Atom, keeping track of how many there are
    int count = 0;
    // Molecules
    for(int i=0; i<molecules.size(); i++){
    // Residues
    for(int j=0; j<molecules[i].residues.size(); j++){
    // Atoms
    for(int k=0; k<molecules[i].residues[j].atoms.size(); k++){
        // Increment the count
        count += 1;
        coors[0] += molecules[i].residues[j].atoms[k].x;
        coors[1] += molecules[i].residues[j].atoms[k].y;
        coors[2] += molecules[i].residues[j].atoms[k].z;}}}
    // If there are no Atoms, end the program
    if(count == 0){string temp = "Docking has failed because one of the groups "
    "of Molecules does not contain any Atoms."; cout << temp << endl;
    exit(EXIT_FAILURE);}
    // Calculate the average coordinates
    for(int i=0; i<3; i++){coors[i] /= float(count);}
    // Move each Molecule
    for(int i=0; i<molecules.size(); i++){molecules[i].move(coors, '-');}
    // Calculate the rotation matrix for X and do that rotation
    float rmatrix [3][3];
    calculate_rmatrix(data[0], 'x', rmatrix);
    for(int i=0; i<molecules.size(); i++){molecules[i].rotate(rmatrix);}
    // Now do Y and Z
    calculate_rmatrix(data[1], 'y', rmatrix);
    for(int i=0; i<molecules.size(); i++){molecules[i].rotate(rmatrix);}
    calculate_rmatrix(data[2], 'z', rmatrix);
    for(int i=0; i<molecules.size(); i++){molecules[i].rotate(rmatrix);}
    // Now we need to do the cartesian perturbation. Those coordinates are being
    // ADDED to the Molecules' current coordinates. So just include them in
    // coors and then move the Molecules back where they were
    coors[0] += data[3]; coors[1] += data[4]; coors[2] += data[5];
    for(int i=0; i<molecules.size(); i++){molecules[i].move(coors, '+');}}

float calculate_energy(vector<MOLECULES::Molecule> &molecules, int I, 
                       vector<vector<MOLECULES::Molecule> > &best, 
                       vector<MOLECULES::DesignGroup> &constants){
    /* Calculate the energy of this docking group with the other Molecules in
     * their Design Groups, including the best copies of the other Molecules
     * that are moving during docking */
    // Store the calculated energy here
    float energies [3] = {0, 0, 0};
    // Loop through the Atoms of the Molecules in the Docking Group
    for(int i=0; i<molecules.size(); i++){
    for(int j=0; j<molecules[i].residues.size(); j++){
    for(int k=0; k<molecules[i].residues[j].atoms.size(); k++){
        // Loop through the constant Design Groups
        for(int w=0; w<constants.size(); w++){
            // Only continue if this Design Group matches the group number of
            // the molecule
            if (constants[w].number == molecules[i].groupNumber){
            // Go through the Molecules, Residues, and Atoms
            for(int x=0; x<constants[w].molecules.size(); x++){
            for(int y=0; y<constants[w].molecules[x].residues.size(); y++){
            for(int z=0; z<constants[w].molecules[x].residues[y].atoms.size(); z++){
                // Calculate the energy between these Atoms
                ENERGY_FUNCTIONS::IE_calculation(molecules[i].residues[j].atoms[k], 
                constants[w].molecules[x].residues[y].atoms[z], energies, 0,
                true);
            }}}}}
        // Also calculate energies with other molecules that move during docking
        // that are NOT in this docking group
        for(int w = 0; w<best.size(); w++){if(I != w){
        for(int x = 0; x<best[w].size(); x++){
        // Only do the calculation if the two Molecules are in the same Design
        // Group
        if(best[w][x].groupNumber == molecules[i].groupNumber){
        for(int y=0; y<best[w][x].residues.size(); y++){
        for(int z=0; z<best[w][x].residues[y].atoms.size(); z++){
        // Calculate the energy between this pair of Atoms
        ENERGY_FUNCTIONS::IE_calculation(molecules[i].residues[j].atoms[k], 
        best[w][x].residues[y].atoms[z], energies, 0, true);
        }}}}}}}}}
    // Calculate the total interaction energy
    float energy = energies[0] + energies[1] + energies[2];
    return energy;}

void copy_molecules(vector<MOLECULES::Molecule> &molecules, int I,
                    vector<vector<MOLECULES::Molecule> > &Docking){
    /* Copy the Molecule structures in molecules into the appropriate copy in
     * the Docking vector. This is used to start an iteration of docking or
     * update the Best copies */
    // Loop through the structures in molecules
    for(int i=0; i<molecules.size(); i++){
        // Copy the Molecule
        Docking[I][i].copy(molecules[i]);}}

int keep_choice(float currentEnergy, float bestEnergy, float temperature, 
                float rnd){
    /* Use simulated annealing to determine whether or not to keep the results
     * of these random motions */
    // Store the choice in this variable
    int choice = -1;
    // If this set of movements is better
    if (currentEnergy <= bestEnergy){choice = 0;}
    // If it is kept by simulated annealing
    else if (rnd <= exp((bestEnergy - currentEnergy)/temperature)){choice = 1;}
    // Otherwise don't keep it
    else {choice = 2;}
    return choice;}

// End the DOCKING FUNCTIONS name space
}
// End the preproccesor commands that make sure this file is only loaded once
#endif
