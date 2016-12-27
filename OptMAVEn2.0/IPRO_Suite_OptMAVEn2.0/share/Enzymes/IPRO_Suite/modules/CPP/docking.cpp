/* Written in 2013 by Robert Pantazes of the Costas Maranas Lab in the Chemical
 * Engineering Department of the Pennsylvania State University
 
 * This file carries out the movements and energy calculations of docking */

// Include general C++ files
# include <iostream>
# include <fstream>
# include <sstream>
# include <cstdlib>
# include <vector>
# include <string>
# include <cmath>
// Include the MOLECULES and ENERGY_FUNCTIONS header files
# include "MOLECULES.h"
# include "ENERGY_FUNCTIONS.h"
// And most importantly, include the DOCKING header file
# include "DOCKING_FUNCTIONS.h"

// use the standard name space
using namespace std;

// Start the program
int main(){
// Include the relevant Docking Functions
using DOCKING_FUNCTIONS::load_docking_data;
using DOCKING_FUNCTIONS::load_constant_molecules;
using DOCKING_FUNCTIONS::calculate_energy;
using DOCKING_FUNCTIONS::copy_molecules;
using DOCKING_FUNCTIONS::move_docking_group;
using DOCKING_FUNCTIONS::keep_choice;
// Create the Best and parameters vectors that will store structures and
// information for docking
vector<vector<MOLECULES::Molecule> > Best;
vector<vector<float> > parameters;
// Load that information
int iterations;
iterations = load_docking_data(Best, parameters);
// Load the constant structures, as well
vector<MOLECULES::DesignGroup> constants;
load_constant_molecules(Best, constants);
// Calculate initial energies for each Docking Group
float bestEnergies [Best.size()];
for(int I=0; I<Best.size(); I++){
    bestEnergies[I] = calculate_energy(Best[I], I, Best, constants);}
// Store an initial copy of these energies
float initialEnergies [Best.size()];
for(int I=0; I<Best.size(); I++){initialEnergies[I] = bestEnergies[I];}
// Initialize an array for the best movements so far as well as whatever the
// kept movements have been. 
float bestMovements [Best.size()][6];
float currentMovements [Best.size()][6];
for(int i=0; i<Best.size(); i++){for(int j=0; j<6; j++){
    bestMovements[i][j] = 0; currentMovements[i][j] = 0;}}
// Create the "Current" copies of the Molecules that is used as a reference
// point at the start of an iteration
vector<vector<MOLECULES::Molecule> > Current;
for(int i=0; i<Best.size(); i++){
    vector<MOLECULES::Molecule> molecules;
    for(int j=0; j<Best[i].size(); j++){
        MOLECULES::Molecule molecule;
        molecule.copy(Best[i][j]);
        molecules.push_back(molecule);}
    Current.push_back(molecules);}

// Keep track of which set of parameters is being used
int P = -1;
// Do the iterations of docking
for(int i=0; i<iterations; i++){
    // Loop through the Docking Groups
    for(int I=0; I<Best.size(); I++){
        // Increment the parameter index
        P += 1;
        // Make a copy of the Molecules that should be used in this iteration
        vector<MOLECULES::Molecule> molecules;
        for(int j=0; j<Best[I].size(); j++){
            MOLECULES::Molecule molecule; 
            molecule.copy(Best[I][j]);
            molecules.push_back(molecule);}
        // Do their random movements
        move_docking_group(molecules, parameters[P]);
        // Calculate their energies with other Molecules
        float energy = calculate_energy(molecules, I, Best, constants);
        // Determine what to do
        int keep = keep_choice(energy, bestEnergies[I], parameters[P][6],
                               parameters[P][7]);
        // If the result is the best so far or kept by simulated annealing,
        // update the Current structures and movements
        if (keep < 2){
            // Copy the Molecules
            copy_molecules(molecules, I, Current);
            // Increment the movements information
            for(int j=0; j<6; j++){currentMovements[I][j] += parameters[P][j];}}
        // If these are the best results so far, update the best information
        if (keep == 0){
            copy_molecules(Current[I], I, Best);
            // Current is updated as it goes, but best should just be assigned
            // the values in current
            for(int j=0; j<6; j++){bestMovements[I][j]=currentMovements[I][j];}
            // And update the bestEnergies
            bestEnergies[I] = energy;}
        // That's all that needs to happen during an iteration, so end those
        }}

// Now we just need to output the relevant information. Open the file
ofstream out1("docking_results.txt");
// Loop through the Design Groups
for(int i=0; i<Best.size(); i++){
    // Output the energy and perturbation information
    for(int j=0; j<6; j++){out1 << bestMovements[i][j] << " ";}
    out1 << bestEnergies[i] << " " << initialEnergies[i] << endl;
    // And output each Molecule's structure
    for(int j=0; j<Best[i].size(); j++){Best[i][j].write();}}
out1.close();
// End the program
return 0;}
