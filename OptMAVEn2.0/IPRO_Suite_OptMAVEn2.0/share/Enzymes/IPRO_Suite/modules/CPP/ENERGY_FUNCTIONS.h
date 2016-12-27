/* Written in 2013 by Robert Pantazes of the Costas Maranas Lab in the Chemical
 * Engineering Department of the Pennsylvania State University
 
 * This file contains functions used to calculate the energy between Atoms in
 * non-bonded energy calculations in the IPRO suite of programs (i.e. during
 * rotamer optimizations and docking). */

// Have some preprocessor commands so that this header file is only included
// once in any C++ program
#ifndef ENERGY_FUNCTIONS_CHECK
#define ENERGY_FUNCTIONS_CHECK 1
// Include C++ modules
# include <cmath>
# include <string>
# include <vector>
# include <cstdlib>
// Include the MOLECULES.h header file, which contains the Atom class
# include "MOLECULES.h"

// Declare the use of the standard namespace
using namespace std;

// Store the contents of this file in a namespace so they can be specifically
// referenced in other C++ files
namespace ENERGY_FUNCTIONS
{

// Declare various functions for calculating the energies between Atoms
float VDW_CHARMM(MOLECULES::Atom& atom1, MOLECULES::Atom& atom2, float dis, 
                 int skip){
    // Calculate the VDW energy between 2 Atoms when the CHARMM force field is
    // being used. There is only a single modification relative to the energy
    // function from CHARMM: the use of a "modify" parameter to "soften" the VDW
    // radii, as this function deals with rotamers.
    // Make sure that both Atoms have the same vdw attribute
    if (atom1.vdw != atom2.vdw){return 0.0;}
    // If VDW energy should not be calculated.
    if (atom1.vdw == false){return 0.0;}
    // Otherwise calculate the energy
    float modify = 1.00;
    float eps0, rmin;
    if (skip != 4){
        eps0 = sqrt(atom1.vdwEpsilon * atom2.vdwEpsilon);
        rmin = modify*(atom1.vdwRadius + atom2.vdwRadius);}
    else {
        eps0 = sqrt(atom1.vdw14E * atom2.vdw14E);
        rmin = modify*(atom1.vdw14R + atom2.vdw14R);}
    float phi = atom1.ljPhi;
    // If the interatomic distance is greater than the energetic minimum
    // distance, then do not soften the radius
    if (dis >= rmin){phi = 1.0;}
    // Make sure the two atoms have the same softening term
    if (atom1.ljPhi != atom2.ljPhi){string temp = "The Lennard-Jones softening "
        "value must be equal between atoms.";
        cout << temp << endl; exit(EXIT_FAILURE);}
    float energy = (-eps0+phi*eps0)+phi*eps0*(pow(rmin/dis, 12) - 2.0*pow(rmin/dis, 6));
    return energy;}

float ELEC_CHARMM(MOLECULES::Atom& atom1, MOLECULES::Atom& atom2, float dis, 
                  float cutoff){
    // Calculate the VDW energy between 2 Atoms when the CHARMM force field is
    // being used. This is the same function as is used in CHARMM.
    // If both Atoms don't have the same elec attribute
    if (atom1.elec != atom2.elec){return 0.0;}
    // If electrostatics energy should not be calculated
    if (atom1.elec == false){return 0.0;}
    // Otherwise calculate the energy
    float CCELEC = 331.843; float eps = 1.0;
    // Calculate the parts of this function
    float numerator = atom1.charge*atom2.charge*CCELEC;
    float denominator = dis*eps;
    float energy = pow(1.0 - pow(dis/cutoff, 2), 2)*numerator/denominator;
    return energy;}

float GB_CHARMM(MOLECULES::Atom& atom1, MOLECULES::Atom& atom2, float dis){
    // Calculate the Generalized-Born implicit solvation between 2 Atoms
    // This function comes from:
    // B. Dominy and C.L. Brooks, III. Development of a Generalized Born Model
    // Parameterization for Proteins and Nucleic Acids. J. Phys. Chem. 103,
    // 3765-3773 (1999).
    // Make sure that both Atoms have identical GB attributes
    if (atom1.gb != atom2.gb){string temp = "A Generalized-Born implicit "
        "solvation energy will not be calculated between two Atoms that do not "
        "have the same gb attribute.";
        cout << temp << endl; exit(EXIT_FAILURE);}
    // If the solvation energy should not be calculated
    if (atom1.gb == false){return 0.0;}
    // otherwise do the calculation
    float CCELEC = 331.843; float eps_prot = 80.0;
    float Dij = pow(dis, 2)/(4.0 * atom1.gbRadius * atom2.gbRadius);
    float numerator = atom1.charge * atom2.charge;
    float denominator = sqrt(pow(dis, 2) + atom1.gbRadius * atom2.gbRadius *
                        exp(-Dij));
    float energy = -CCELEC * (1.0 - (1.0/eps_prot)) * numerator/denominator;
    return energy;}

float LK_CHARMM(MOLECULES::Atom& atom1, MOLECULES::Atom& atom2, float dis){
    // Calculate the Lazaridis-Karplus implicit solvation between 2 Atoms
    // Make sure the Atoms have the same lk attribute
    if (atom1.lk != atom2.lk){return 0.0;}
    // If LK solvation should not be calculated, either because it is not being
    // used or either Atom is a hydrogen
    if ((atom1.lk == false) || ((atom1.lkRadius <= 1.34) || (atom2.lkRadius <=
        1.34))){return 0.0;}
    // Otherwise calculate the energy
    float x01 = (dis - atom1.lkRadius)/atom1.lkLambda;
    float x02 = (dis - atom2.lkRadius)/atom2.lkLambda;
    float t1 = -0.08979/pow(dis, 2);
    float t2 = exp(-pow(x01, 2))*atom1.lkGibbs*atom2.lkVolume/atom1.lkLambda;
    float t3 = exp(-pow(x02, 2))*atom2.lkGibbs*atom1.lkVolume/atom2.lkLambda;
    float energy = t1*(t2 + t3);
    return energy;}

int exclusions(MOLECULES::Atom& atom1, MOLECULES::Atom& atom2){
    // Determine whether non-bonded energy values should be calculated between 2
    // Atoms. This is important during rotamer energy calculations
    // Have different methods for each supported file format
    if (atom1.fileFormat != atom2.fileFormat){string temp = "The exclusions C++"
        " function cannot work for two C++ Atoms that do not have the same file"
        " format."; cout << temp << endl; exit(EXIT_FAILURE);}
    // The PDB file format
    if (atom1.fileFormat == "PDB"){
        // If the two Atom's aren't in the same Molecule
        if (atom1.moleculeName != atom2.moleculeName){return 0;}
        // The order of the Atoms matters - Atom1 must be from a rotamer
        // If Atom2 is in the Residue preceeding Atom1
        if (atom2.residueNumber + 1 == atom1.residueNumber){
            // Rotamers are only side chain Atoms, so only Atoms in the beta
            // list can be in range of the C in the previous residue
            if ((atom2.name == "C") && ((atom1.name == "CB") ||
                (atom1.name == "HA1"))){return 4;}
            // Proline has special handling because its CD is bound to its N
            if (atom1.residueKind == "PRO"){
                if (atom2.name == "C"){
                    if (atom1.name == "CD"){return 3;}
                    else if ((atom1.name == "CG") || ((atom1.name == "HD1") ||
                             (atom1.name == "HD2"))){return 4;}
                    }
                else if ((atom1.name == "CD") && ((atom2.name == "CA") ||
                         (atom2.name == "O"))){return 4;}
                }}
        // If Atom2 is in the Residue after Atom1
        else if (atom2.residueNumber - 1 == atom1.residueNumber){
            // If Atom1 is in the beta position, see if atom2 is N
            if ((atom2.name == "N") && ((atom1.name == "CB") ||
                (atom1.name == "HA1"))){return 4;}
            }
        // If they are in the same Residue
        else if (atom2.residueNumber == atom1.residueNumber){
            // Instead of LONG nested if statements, search through arrays of
            // the permitted Atoms, defined based on their distance from the
            // alpha carbon
            string alpha [1] = {"CA"};
            string next [4] = {"C", "N", "HA", "HA2"};
            string third [9] = {"O", "OT1", "OT2", "HN", "HN1", "HN2", "HT1",
                                "HT2", "HT3"};
            string beta [2] = {"CB", "HA1"};
            string gamma [10] = {"CG", "CG1", "CG2", "HB", "HB1", "HB2", "HB3", 
                                 "OG", "OG1", "SG"};
            string delta [17] = {"CD", "CD1", "CD2", "HG", "HG1", "HG2", "HG11",
                                 "HG12", "HG13", "HG21", "HG22", "HG23", "OD1", 
                                 "OD2", "ND1", "ND2", "SD"};
            // Handle special cases for proline first
            if (atom1.residueKind == "PRO"){
                // N terminal Residues that are not proline but are having
                // proline considered as a rotamer will have an extra Atom, HT3.
                // Skip this without question, using the smallest return value
                // possible that indicates skipping an Atom.
                if (atom2.name == "HT3"){return 2;}
                // The normal search starts from the alpha carbon. Do this one
                // from the N
                if (atom2.name == "N"){
                    if (atom1.name == "CD"){return 2;}
                    else if ((atom1.name == "CG") || ((atom1.name == "HD1") ||
                             (atom1.name == "HD2"))){return 3;}
                    // The CB is a 3, but that will get picked up in the normal
                    // check
                    else if ((atom1.name == "HG1") || (atom1.name == "HG2")){
                        return 4;}
                    }
                // Atoms bonded to the N Atom
                else if (((atom2.name == "CA") || ((atom2.name == "HN") ||
                    (atom2.name == "HN1"))) || ((atom2.name == "HN2") ||
                    ((atom2.name == "HT1") || (atom2.name == "HT2")))){
                    if (atom1.name == "CD"){return 3;}
                    else if ((atom1.name == "CG") && (atom2.name == "CA")){
                        return 3;}
                    else if (((atom1.name == "HD1") || (atom1.name == "HD2"))
                    || (atom1.name == "CG")){
                        return 4;}
                    }
                // Backbone Atoms bonded to those Atoms in THIS Residue
                else if ((atom2.name == "HA") || (atom2.name == "HA2")){
                    if (atom1.name == "CD"){return 4;}
                    }}
            // Now that special Proline exceptions have been identified, do the
            // normal comparisons for amino acids
            // If Atom 1 is at the beta position
            for(int i=0; i<2; i++){
                if (beta[i] == atom1.name){
                    // If Atom2 is at the alpha position
                    if (atom2.name == "CA"){return 2;}
                    // Or bonded to the CA
                    for (int j=0; j<4; j++){if (next[j]==atom2.name){return 3;}}
                    // Or bonded to those Atoms and in this Residue
                    for (int j=0; j<9; j++){if(third[j]==atom2.name){return 4;}}
                    return 0;}}
            // If Atom1 is at the gamma position
            for(int i=0; i<10; i++){
                if (gamma[i] == atom1.name){
                    if (atom2.name == "CA"){return 3;}
                    for(int j=0;j<4;j++){if (next[j]==atom2.name){return 4;}}
                    return 0;}}
            // Delta position
            for(int i=0; i<17; i++){
                if (delta[i] == atom1.name){
                    if (atom2.name == "CA"){return 4;}
                    return 0;}}
            return 0;
            // Close the if statement for equal Residues and the one for the PDB
            // file format
            }}
    // If the file format is not supported
    else{string temp = "The C++ exclusions function does not support the "
         + atom1.fileFormat + " file format.";
         cout << temp << endl; exit(EXIT_FAILURE);}
    // End the function
    return 0;}

void IE_calculation(MOLECULES::Atom& atom1, MOLECULES::Atom& atom2, 
                    float energies [], int skip = 0, bool solvation = true){
    // Calculate the interaction energy between two Atoms
    // First, confirm that both Atoms have the same force field attribute
    if (atom1.forceField != atom2.forceField){string temp = "The "
        "IE_calculation function requires that both C++ Atoms given to it use "
        "the same force field attribute.";
        cout << temp << endl; exit(EXIT_FAILURE);}
    // Have a different method for each supported force field
    if (atom1.forceField == "CHARMM"){
        // Calculate the distance between the Atoms
        float dis = sqrt(pow(atom1.x - atom2.x, 2) + pow(atom1.y - atom2.y, 2) +
                         pow(atom1.z - atom2.z, 2));
        // If appropriate, calculate VDW energies
        if (((atom1.vdw) && (dis <= 12.0)) && ((skip == 0) || (skip == 4))){
            energies[0] += VDW_CHARMM(atom1, atom2, dis, skip);}
        // Electrostatics
        if (((atom1.elec) && (dis <= 12.0)) && ((skip == 0) || (skip == 4))){
            energies[1] += ELEC_CHARMM(atom1, atom2, dis, 12.0);}
        // Implicit solvation
        if (((atom1.gb) && (dis <= 12.0)) && ((skip == 0) || (skip == 4)) && solvation){
            energies[2] += GB_CHARMM(atom1, atom2, dis);}
        else if (((atom1.lk) && (dis <= 9.0)) && ((skip == 0) || (skip == 4)) && solvation){
            energies[2] += LK_CHARMM(atom1, atom2, dis);}
        // End the CHARMM portion of this function
        }
    // If the force field is not supported
    else {string temp = "The IE_calculation function does not support the "
        + atom1.forceField + " force field.";
        cout << temp << endl; exit(EXIT_FAILURE);}
    // End the function
    }

// close the ENERGY_FUNCTIONS namespace
}
// End the preprocessor commands that prevent this header file from being
// included more than once in any C++ program
#endif
