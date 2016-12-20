/* Written in 2013 by Robert Pantazes of the Costas Maranas Lab in the Chemical
 * Engineering Department of the Pennsylvania State University
 
 * This file contains classes that are similar to the ones in the MOLECULES.py
 * module, but these are for use in non-bonded energy calculations that are not
 * easy to do in CHARMM. */

// Have preprocessor commands so that this header file is only included once in
// any C++ program
#ifndef MOLECULE_CHECK
#define MOLECULE_CHECK 1
// Include C++ modules
# include <iostream>
# include <sstream>
# include <cstdlib>
# include <fstream>
# include <string>
# include <vector>
# include <cmath>

// Declare the use of the standard namespace
using namespace std;

// Store the contents of this file in a unique namespace
namespace MOLECULES
{
// Create a template function for converting integers and floats to strings,
// because I can't figure out how C++ regularly does that
template <class T>
string make_string(T number){
    // Convert a number to a string
    // Declare a string stream
    stringstream s;
    // Put the number into the string stream
    s << number;
    // Put that into a string
    string text; s >> text;
    return text;}

// The Atom class contains information about a single Atom.
class Atom{
    // Everything in the Atom is public and can be accessed by any function
    public:
        // Every Atom will have a file format and a force field, as well as
        // coordinates
        string fileFormat, forceField;
        float x, y, z;
        // Attributes associated with the PDB file format
        string moleculeName, residueKind, name;
        int rotamerNumber, residueNumber;
        // Attributes associated with the CHARMM force field.
        bool vdw, elec, gb, lk;
        float vdwRadius, vdwEpsilon, charge, gbRadius;
        float lkRadius, lkLambda, lkVolume, lkGibbs, ljPhi;
        // And the 1-4 VDW terms we didn't used to have
        float vdw14R, vdw14E;
        // Now declare functions of the Atom class

        // The function that loads the contents of the Atom from a string of
        // text
        void load(string text){
            // Set the boolean values to false
            // Put the string of text into a stringstream
            stringstream line(text);
            // Read in the attributes of the Atom, starting with its file format
            line >> fileFormat;
            // Have a different method for each supported file format
            if (fileFormat == "PDB"){
                line >> moleculeName >> residueNumber >> residueKind;
                line >> rotamerNumber >> name >> x >> y >> z;}
            // If the file format is not supported
            else {string temp = "The load function of the C++ Atom class does "
                "not support the " + fileFormat + " file format.";
                cout << temp << endl; exit(EXIT_FAILURE);}
            // Now get the force field
            line >> forceField;
            // Have a unique method for each forceField
            if (forceField == "CHARMM"){
                // Set the boolean values to a default of false
                vdw = false; elec = false; gb = false; lk = false; 
                // Read until the end of the stream. Store the parameter type in
                // this string
                string type;
                line >> type;
                while (!line.eof()){
                    // Respond appropriately based on the content of type
                    if (type == "VDW"){
                        vdw = true; 
                        line >> vdwRadius >> vdwEpsilon >> vdw14R >> vdw14E >>
                        ljPhi;}
                    else if(type == "ELEC"){elec = true; line >> charge;}
                    else if(type == "GB"){gb = true; line >> gbRadius;}
                    else if(type == "LK"){
                        lk = true; 
                        line >> lkRadius >> lkLambda >> lkVolume >> lkGibbs;}
                    // If the type is something that isn't recognized
                    else {string temp = "The load function of the C++ Atom "
                        "class does not recognize " + type + " as an energy "
                        "term.";
                        cout << temp << endl; exit(EXIT_FAILURE);}
                    // Read in the next energy function type
                    line >> type;}}
            // If the force field is not supported
            else {string temp = "The load function of the C++ Atom class does "
                "not support the " + forceField + " force field.";
                cout << temp << endl; exit(EXIT_FAILURE);}
            // End the load function
            }

        // Have an explicit function to copy an Atom's contents to another Atom
        void copy(Atom& atom){
            // Copy the coordinates
            x = atom.x; y = atom.y; z = atom.z;
            // Copy the file format
            fileFormat = atom.fileFormat;
            // Copy attributes based on the file format
            if (fileFormat == "PDB"){
                moleculeName = atom.moleculeName; 
                residueKind = atom.residueKind;
                residueNumber = atom.residueNumber;
                rotamerNumber = atom.rotamerNumber;
                name = atom.name;}
            else {string temp = "The copy function of the C++ Atom class does "
                "not support the " + fileFormat + " file format.";
                cout << temp << endl; exit(EXIT_FAILURE);}
            // The force field
            forceField = atom.forceField;
            // Copy attributes based on the force field
            if (forceField == "CHARMM"){
                // Copy the boolean values
                vdw = atom.vdw; elec = atom.elec;
                gb = atom.gb; lk = atom.lk;
                // Copy the appropriate attributes
                if (vdw){
                    vdwRadius = atom.vdwRadius; vdwEpsilon = atom.vdwEpsilon;
                    vdw14R = atom.vdw14R; vdw14E = atom.vdw14E; 
                    ljPhi = atom.ljPhi;}
                if (elec){charge = atom.charge;}
                if (gb){gbRadius = atom.gbRadius;}
                if (lk){
                    lkRadius = atom.lkRadius; lkLambda = atom.lkLambda;
                    lkVolume = atom.lkVolume; lkGibbs = atom.lkGibbs;}}
            else {string temp = "The copy function of the C++ Atom class does "
                "not support the " + forceField + " force field.";
                cout << temp << endl; exit(EXIT_FAILURE);}
            // End the copy function
            }
        
        // Create a summary output of the Atom
        string output(){
            // Store the contents of the Atom in this line
            string line = "";
            // Do this differently for each supported file format
            if (fileFormat == "PDB"){
                // Include the labelling information of the Atom
                line += moleculeName + " " + make_string(residueNumber) + " " +
                        residueKind + " " + make_string(rotamerNumber) + " " +
                        name + " " + make_string(x) + " " + make_string(y) +
                        " " + make_string(z) + "\n";
                return line;}
            // If the file format is not supported, throw an exception
            else {string temp = "The output function of the C++ Atom class does"
                " not support the " + fileFormat + " file format.";
                cout << temp << endl; exit(EXIT_FAILURE);}
            // End the output function
            }

        // Have a function for moving the Atom
        void move(float coors [], char sign = '-'){
            // Do this differently for different signs
            if (sign == '-'){x -= coors[0]; y -= coors[1]; z -= coors[2];}
            else if(sign == '+'){x += coors[0]; y += coors[1]; z += coors[2];}
            else {string temp = "The move function of the C++ Atom class does "
                "not recognize this sign: ";
                cout << temp << sign << endl; exit(EXIT_FAILURE);}}

        // And for rotating the Atom
        void rotate(float rmatrix [][3]){
            // Calculate new coordinates for the Atom
            float a = rmatrix[0][0]*x + rmatrix[0][1]*y + rmatrix[0][2]*z;
            float b = rmatrix[1][0]*x + rmatrix[1][1]*y + rmatrix[1][2]*z;
            float c = rmatrix[2][0]*x + rmatrix[2][1]*y + rmatrix[2][2]*z;
            // Store the new coordinates
            x = a; y = b; z = c;}

        // End the declaration of the Atom class
        };

class Residue{
// Everything in the Residue class is public
public:
    // Every Residue will have a file format and a force field
    string fileFormat, forceField;
    // Every Residue will have a vector of Atoms
    vector<Atom> atoms;
    // Attributes associated with the PDB file format
    string moleculeName, kind;
    int rotamerNumber, number;
    // Attributes associated with the CHARMM force field
    bool vdw, elec, gb, lk;

    // Have functions for the Residue class. Start with one that validates that
    // Atoms are acceptable
    void validate_Atoms(vector<Atom>& ATOMS){
        // Check attributes associated with the file format
        if (fileFormat == "PDB"){
            if (ATOMS.size() == 0){string temp = "The PDB file format requires "
                "that a C++ Residue contain at least one C++ Atom.";
                cout << temp << endl; exit(EXIT_FAILURE);}
            // Each Atom must use this file format, have a unique name, and
            // match the Residue's moleculeName, kind, number, and rotamerNumber
            // Loop through the Atoms
            for(int i=0; i<ATOMS.size(); i++){
                // check to see if it uses the PDB file format
                if (ATOMS[i].fileFormat != "PDB"){string temp = "The PDB file "
                    "format requires that each C++ Atom in a C++ Residue use "
                    "the Residue's file format attribute.";
                    cout << temp << endl; exit(EXIT_FAILURE);}
                // Check the Atom's name
                bool flag = false;
                for(int j=0; j<i; j++){
                    if (ATOMS[j].name == ATOMS[i].name){flag = true; break;}}
                if (flag){string temp = "The PDB file format requires that each"
                    " C++ Atom in a C++ Residue have a unique name attribute.";
                    cout << temp << endl; exit(EXIT_FAILURE);}
                // moleculeName
                if (ATOMS[i].moleculeName != moleculeName){string temp = "The "
                    "PDB file format requires that each C++ Atom in a C++ "
                    "Residue match the Residue's moleculeName attribute.";
                    cout << temp << endl; exit(EXIT_FAILURE);}
                // kind
                if (ATOMS[i].residueKind != kind){string temp = "The PDB file "
                    "format requires that each C++ Atom in a C++ Residue use "
                    "the Residue's kind attribute.";
                    cout << temp << endl; exit(EXIT_FAILURE);}
                // number
                if (ATOMS[i].residueNumber != number){string temp = "The PDB file "
                    "format requires that each C++ Atom in a C++ Residue use "
                    "the Residue's number attribute.";
                    cout << temp << endl; exit(EXIT_FAILURE);}
                // rotamer number
                if (ATOMS[i].rotamerNumber != rotamerNumber){string temp = "The"
                    " PDB file format requires that each C++ Atom in a C++ "
                    "Residue use the Residue's rotamerNumber attribute.";
                    cout << temp << endl; exit(EXIT_FAILURE);}
                // End the for loop and the if statement
                }}
        // If the file format is not supported, exit the program
        else {string temp = "The validate_Atoms function in the C++ Residue "
            "class does not support the " + fileFormat + " file format.";
            cout << temp << endl; exit(EXIT_FAILURE);}
        // Check force field attributes
        if (forceField == "CHARMM"){
            // Check the boolean values of the Atoms
            for(int i=0; i<ATOMS.size(); i++){
                if (ATOMS[i].vdw != vdw){string temp = "The CHARMM force field "
                    "requires that each C++ Atom in a C++ Residue match the "
                    "Residue's vdw attribute.";
                    cout << temp << endl; exit(EXIT_FAILURE);}
                if (ATOMS[i].elec != elec){string temp = "The CHARMM force "
                    "field requires that each C++ Atom in a C++ Residue use the"
                    " Residue's elec attribute.";
                    cout << temp << endl; exit(EXIT_FAILURE);}
                if (ATOMS[i].gb != gb){string temp = "The CHARMM force field "
                    "requires that each C++ Atom in a C++ Residue match the "
                    "Residue's gb attribute.";
                    cout << temp << endl; exit(EXIT_FAILURE);}
                if (ATOMS[i].lk != lk){string temp = "The CHARMM force field "
                    "requires that each C++ Atom in a C++ Residue match the "
                    "Residue's lk attribute.";
                    cout << temp << endl; exit(EXIT_FAILURE);}
                // End the for loop and if statement
                }}
        // If the force field is not supported
        else {string temp = "The validate_Atoms function in the C++ Residue "
            "class does not support the " + forceField + " force field.";
            cout << temp << endl; exit(EXIT_FAILURE);}
        // End the function
        }

    // Have a function to load a Residue from a group of Atoms
    void load(vector<Atom>& ATOMS){
        // Confirm that there ARE atoms
        if (ATOMS.size() == 0){string temp = "The load function of the C++ "
            "Residue class requires the Residue to contain at least one C++ "
            "Atom"; cout << temp << endl; exit(EXIT_FAILURE);}
        // Since there are Atoms, store the file format and force field
        fileFormat = ATOMS[0].fileFormat;
        forceField = ATOMS[0].forceField;
        // Store additional attributes based on the file format
        if (fileFormat == "PDB"){
            kind = ATOMS[0].residueKind; 
            moleculeName = ATOMS[0].moleculeName;
            number = ATOMS[0].residueNumber;
            rotamerNumber = ATOMS[0].rotamerNumber;}
        // If the file format is not supported
        else {string temp = "The load function of the C++ Residue class does "
            "not support the " + fileFormat + " file format.";
            cout << temp << endl; exit(EXIT_FAILURE);}
        // Store other attributes based on the force field
        if (forceField == "CHARMM"){
            vdw = ATOMS[0].vdw; elec = ATOMS[0].elec;
            gb = ATOMS[0].gb; lk = ATOMS[0].lk;}
        else {string temp = "The load function of the C++ Residue class does "
            "not support the " + forceField + " force field.";
            cout << temp << endl; exit(EXIT_FAILURE);}
        // Validate the Atoms
        validate_Atoms(ATOMS);
        // With the Atoms validated, copies of them can be stored
        atoms.clear();
        for (int i=0; i<ATOMS.size(); i++){
            Atom atom; atom.copy(ATOMS[i]); atoms.push_back(atom);}
        // End the load function
        }

    // Explicitly control copying of a Residue
    void copy(Residue& residue){
        // Copy the Residue's file format and force field
        fileFormat = residue.fileFormat;
        forceField = residue.forceField;
        // Have different methods for each supported file format
        if (fileFormat == "PDB"){
            number = residue.number; rotamerNumber = residue.rotamerNumber;
            moleculeName = residue.moleculeName; kind = residue.kind;}
        // If the file format is not supported
        else {string temp = "The copy function of the C++ Residue class does "
            "not support the " + forceField + " force field.";
            cout << temp << endl; exit(EXIT_FAILURE);}
        // And by force field
        if (forceField == "CHARMM"){
            vdw = residue.vdw; elec = residue.elec;
            gb = residue.gb; lk = residue.lk;}
        else {string temp = "The copy function of the C++ Residue class does "
            "not support the " + forceField + " force field.";
            cout << temp << endl; exit(EXIT_FAILURE);}
        // Validate the Residue's AToms
        validate_Atoms(residue.atoms);
        // Clear the current atoms vector, than copy the contents of the residue
        atoms.clear();
        for(int i=0; i<residue.atoms.size(); i++){
            Atom atom; atom.copy(residue.atoms[i]); atoms.push_back(atom);}
        // End the function
        }

    // Have a function that creates a string of output for the Residue's Atoms
    string output(){
        // Store the text here
        string text = "";
        // Loop through the Atoms
        for(int i=0; i<atoms.size(); i++){text += atoms[i].output();}
        return text;}

    // Have functions for moving and rotating the Residue
    void move(float coordinates [], char sign = '-'){
        for(int i=0; i<atoms.size(); i++){atoms[i].move(coordinates, sign);}}

    void rotate(float rmatrix [][3]){
        for(int i=0; i<atoms.size(); i++){atoms[i].rotate(rmatrix);}}
    // That's the end of the Residue class
    };

class Molecule{
// Everything in the Molecule class is public
public:
    // Every Molecule has a file format and force field
    string fileFormat, forceField;
    // The Residues that make up the Molecule
    vector<Residue> residues;
    // Attributes associated with the PDB file format
    string name;
    // Attributes associated with the CHARMM force field
    bool vdw, elec, gb, lk;
    // The Molecule's file name is only used during docking
    string fileName;
    // Same thing for its group number, which is only a pretend string
    string groupNumber;

    void validate_Residues(vector<Residue>& RESIDUES){
        // Make sure that the RESIDUES vector is not empty
        if (RESIDUES.size() == 0){string temp = "The validate_Residues function"
            " of the C++ Molecule class requires that the Molecule contain at "
            "least one C++ Residue.";
            cout << temp << endl; exit(EXIT_FAILURE);}
        // Check the Molecule's Residues for compatibility with the Molecule.
        // Have different checks for each supported
        if (fileFormat == "PDB"){
            // Loop through the Residues
            for(int i=0; i<RESIDUES.size(); i++){
                // Confirm that the Residue uses the right file format
                if (RESIDUES[i].fileFormat != fileFormat){string temp = "The "
                    "PDB file format requires that all C++ Residues in a C++ "
                    "Molecule use the Molecule's file format attribute.";
                    cout << temp << endl; exit(EXIT_FAILURE);}
                // Check that the Residue's moleculeName is right
                if (RESIDUES[i].moleculeName != name){string temp = "The PDB "
                    "file format requires that each C++ Residue in a C++ "
                    "Molecule use the Molecule's name attribute.";
                    cout << temp << endl; exit(EXIT_FAILURE);}
                // Check that the Residue's number and rotamer number
                // combination are unique
                for(int j=0; j<i; j++){
                    if ((RESIDUES[i].number == RESIDUES[j].number) && 
                       (RESIDUES[i].rotamerNumber ==RESIDUES[j].rotamerNumber)){
                        string temp = "The PDB file format requires that each "
                        "C++ Residue in a C++ Molecule have a unique "
                        "combination of name and rotamerNumber attributes.";
                        cout << temp << endl; exit(EXIT_FAILURE);}}
                // End the for loop and the if statement
                }}
        // If the file format is not supported
        else {string temp = "The validate_Residue function of the C++ Molecule "
            "class does not support the " + fileFormat + " file format.";
            cout << temp << endl; exit(EXIT_FAILURE);}
        // Also check by force field
        if (forceField == "CHARMM"){
            // Loop through the Residues
            for(int i=0; i<RESIDUES.size(); i++){
                // Confirm that the vdw, elec, gb, and lk attributes match
                if (RESIDUES[i].vdw != vdw){string temp = "The CHARMM force "
                    "field requires that all C++ Residues in a C++ Molecule "
                    "match the Molecule's vdw attribute.";
                    cout << temp << endl; exit(EXIT_FAILURE);}
                if (RESIDUES[i].elec != elec){string temp = "The CHARMM force "
                    "field requires that all C++ Residues in a C++ Molecule "
                    "match the Molecule's elec attribute.";
                    cout << temp << endl; exit(EXIT_FAILURE);}
                if (RESIDUES[i].gb != gb){string temp = "The CHARMM force field"
                    " requires that all C++ Residues in a C++ Molecule match "
                    "the Molecule's gb attribute.";
                    cout << temp << endl; exit(EXIT_FAILURE);}
                if (RESIDUES[i].lk != lk){string temp = "The CHARMM force field"
                    " requires that all C++ Residues in a C++ Molecule match "
                    "the Molecule's lk attribute.";
                    cout << temp << endl; exit(EXIT_FAILURE);}
                // End the for loop and if statement
                }}
        // If the force field is not supported
        else {string temp = "The validate_Residues function of the C++ Molecule"
            " class does not support the " + forceField + " force field.";
            cout << temp << endl; exit(EXIT_FAILURE);}
        // Now go through the Residues and validate each one, to make sure all
        // the Atoms are OK
        for(int i=0; i<RESIDUES.size(); i++){
            RESIDUES[i].validate_Atoms(RESIDUES[i].atoms);}
        // That ends the function
        }

    // Have a function to load a Molecule from a set of Atoms
    void load(vector<Atom>& ATOMS){
        // Confirm that there are Atoms
        if (ATOMS.size() == 0){string temp = "The load function of the C++ "
            "Molecule class requires that a Molecule contain at least one C++ "
            "Atom."; cout << temp << endl; exit(EXIT_FAILURE);}
        // Extract the first Atom's file format, force field, and other
        // information
        fileFormat = ATOMS[0].fileFormat; forceField = ATOMS[0].forceField;
        // Attributes based on the file format
        if (fileFormat == "PDB"){name = ATOMS[0].moleculeName;}
        else {string temp = "The load function of the C++ Molecule class does "
            "not support the " + fileFormat + " file format.";
            cout << temp << endl; exit(EXIT_FAILURE);}
        // Attributes based on the force field
        if (forceField == "CHARMM"){
            vdw = ATOMS[0].vdw; elec = ATOMS[0].elec;
            gb = ATOMS[0].gb; lk = ATOMS[0].lk;}
        else {string temp = "The load function of the C++ Molecule class does "
            "not support the " + forceField + " force field.";
            cout << temp << endl; exit(EXIT_FAILURE);}
        // Clear the residues vector
        residues.clear();
        // Create Residues
        vector<Atom> residue;
        // Have a different method for each supported file format
        if (fileFormat == "PDB"){
            int n1 = ATOMS[0].residueNumber;
            int n2 = ATOMS[0].rotamerNumber;
            // Loop through the Atoms
            for(int i=0; i<ATOMS.size(); i++){
                // If this Atom belongs in a new residue
                if ((ATOMS[i].residueNumber != n1)||(ATOMS[i].rotamerNumber!=n2)){
                    // Update the numbering information
                    n1 = ATOMS[i].residueNumber;
                    n2 = ATOMS[i].rotamerNumber;
                    // Make a Residue, load the Atoms, and store it
                    Residue res; res.load(residue); residues.push_back(res);
                    // Clear the residue vector
                    residue.clear();}
                // Store the Atom in the residue vector
                residue.push_back(ATOMS[i]);}
            // Store the last set of Atoms
            Residue res; res.load(residue); residues.push_back(res);
            // End the PDB method of the function
            }
        else {string temp = "The load function of the C++ Molecule class does "
            "not support the " + fileFormat + " file format.";
            cout << temp << endl; exit(EXIT_FAILURE);}
        // Validate the Molecule's Residues
        validate_Residues(residues);}

    void copy(Molecule& molecule){
        // Copy the contents of one Molecule to another. Start with the file
        // format
        fileFormat = molecule.fileFormat;
        // Include the Molecule's group number and file name
        fileName = molecule.fileName;
        groupNumber = molecule.groupNumber;
        // Copy attributes based on the value of the file format
        if (fileFormat == "PDB"){name = molecule.name;}
        // If the file format is not supported
        else {string temp = "The copy function of the C++ Molecule class does "
            "not support the " + fileFormat + " file format.";
            cout << temp << endl; exit(EXIT_FAILURE);}
        // Copy attributes based on the force field
        forceField = molecule.forceField;
        if (forceField == "CHARMM"){
            vdw = molecule.vdw; elec = molecule.elec;
            gb = molecule.gb; lk = molecule.lk;}
        // If the force field is not supported
        else {string temp = "The copy function of the C++ Molecule class does "
            "not support the " + forceField + " force field.";
            cout << temp << endl; exit(EXIT_FAILURE);}
        // Validate the residues of the molecule
        validate_Residues(molecule.residues);
        // clear the residues vector
        residues.clear();
        // Loop through the Molecule's Residues
        for (int i=0; i<molecule.residues.size(); i++){
            // copy the Residue
            Residue res; res.copy(molecule.residues[i]);
            // Store this new Residue
            residues.push_back(res);}
        // End the function
        }

    string output(){
        // Create a formatted string for the Molecule's Atoms.
        string text = "";
        for (int i=0; i<residues.size(); i++){text += residues[i].output();}
        return text;}

    void write(){
        // Write the Molecule's data to its file
        string text = output();
        // Create the name of the file
        char outputName [fileName.size()];
        for(int i=0; i<fileName.size(); i++){outputName[i] = fileName[i];}
        ofstream out1(outputName);
        out1 << text; out1.close();}

    void move(float coordinates [], char sign = '-'){
        // Move the Molecule
        for(int i=0; i<residues.size(); i++){
            residues[i].move(coordinates, sign);}}

    void rotate(float rmatrix [][3]){
        for(int i=0; i<residues.size(); i++){residues[i].rotate(rmatrix);}}
    // End the Molecule class declaration
    };

class DesignGroup{
// Everything in the Design Group is public
public:
    // Every DesignGroup has a number, but they're stored as strings
    string number;
    // They also have a file format and force field
    string fileFormat, forceField;
    // There are no attributes associated with the PDB file format
    // These are the attributes associated with CHARMM
    bool vdw, elec, gb, lk;
    // Store the Molecules in this vector
    vector<Molecule> molecules;

    // Have a function for validating that every Molecule in a DesignGroup is
    // acceptable
    void validate_Molecules(vector<Molecule>& MOLECULES){
        // Make sure the MOLECULES vector is not empty
        if (MOLECULES.size() == 0){string temp = "The validate_Molecules "
            "function of the C++ DesignGroup class requires that the "
            "DesignGroup contain at least one C++ Molecule.";
            cout << temp << endl; exit(EXIT_FAILURE);}
        // Have a different method for each supported file format
        if (fileFormat == "PDB"){
            // Loop through the Molecules and make sure they all use the PDB
            // file format
            for(int i=0; i<MOLECULES.size(); i++){
                if (MOLECULES[i].fileFormat != fileFormat){string temp = "The "
                    "PDB file format requires that each C++ Molecule in a C++ "
                    "DesignGroup use the DesignGroup's file format attribute.";
                    cout << temp << endl; exit(EXIT_FAILURE);}
                // Also make sure that they each have a unique name
                for(int j=0; j<i; j++){
                    if (MOLECULES[j].name == MOLECULES[i].name){
                        string temp = "The PDB file format requires that each "
                        "C++ Molecule in a C++ DesignGroup have a unique name.";
                        cout << temp << endl; exit(EXIT_FAILURE);}}
            // End the for loop and the if statement
            }}
        else {string temp = "The validate_Molecules function of the C++ "
            "DesignGroup class does not support the " + fileFormat + " file "
            "format."; cout << temp << endl; exit(EXIT_FAILURE);}
        // Check force field related attributes
        if (forceField == "CHARMM"){
            for(int i=0; i<MOLECULES.size(); i++){
                if (MOLECULES[i].vdw != vdw){string temp = "The CHARMM force "
                    "field requires that each C++ Molecule in a C++ DesignGroup"
                    " use the DesignGroup's vdw attribute.";
                    cout << temp << endl; exit(EXIT_FAILURE);}
                if (MOLECULES[i].elec != elec){string temp = "The CHARMM force "
                    "field requires that each C++ Molecule in a C++ DesignGroup"
                    " use the DesignGroup's elec attribute.";
                    cout << temp << endl; exit(EXIT_FAILURE);}
                if (MOLECULES[i].gb != gb){string temp = "The CHARMM force "
                    "field requires that each C++ Molecule in a C++ DesignGroup"
                    " use the DesignGroup's gb attribute.";
                    cout << temp << endl; exit(EXIT_FAILURE);}
                if (MOLECULES[i].lk != lk){string temp = "The CHARMM force "
                    "field requires that each C++ Molecule in a C++ DesignGroup"
                    " use the DesignGroup's lk attribute.";
                    cout << temp << endl; exit(EXIT_FAILURE);}
            // End the for loop and if statement
            }}
        else {string temp = "The validate_Molecules function of the C++ "
            "DesignGroup class does not support the " + forceField + " force "
            "field."; cout << temp << endl; exit(EXIT_FAILURE);}
        // Validate each individual Molecule
        for(int i=0; i<MOLECULES.size(); i++){
            MOLECULES[i].validate_Residues(MOLECULES[i].residues);}
        // End the function
        }

    // Have a function to load a DesignGroup from a set of Atoms
    void load(vector<Atom>& ATOMS, string NUM = "1"){
        // Confirm that the ATOMS vector is not empty
        if (ATOMS.size() == 0){string temp = "The load function of the C++ "
            "DesignGroup class requires that the DesignGroup contain at least "
            "one C++ Atom."; cout << temp << endl; exit(EXIT_FAILURE);}
        // Store the number
        number = NUM;
        // Store appropriate attributes from the first Atom
        fileFormat = ATOMS[0].fileFormat;
        // There are currently no DesignGroup attributes that are dependent on
        // the file format
        forceField = ATOMS[0].forceField;
        if (forceField == "CHARMM"){
            vdw = ATOMS[0].vdw; elec = ATOMS[0].elec;
            gb = ATOMS[0].gb; lk = ATOMS[0].lk;}
        else {string temp = "The load function of the C++ DesignGroup class "
            "does not support the " + forceField + " force field.";
            cout << temp << endl; exit(EXIT_FAILURE);}
        // Separate the Atoms into Molecules, doing so differenlty depending on
        // the file format
        vector<Atom> molecule;
        if (fileFormat == "PDB"){
            // Get the name of the first Atom's molecule
            string name = ATOMS[0].moleculeName;
            // Loop through the Atoms
            for (int i=0; i<ATOMS.size(); i++){
                // If this Atom belongs in a new Molecule
                if (ATOMS[i].moleculeName != name){
                    // Store the new molecule name
                    name = ATOMS[i].moleculeName;
                    // Create a Molecule out of the molecule vector
                    Molecule mol; mol.load(molecule); molecules.push_back(mol);
                    // Clear the molecule vector
                    molecule.clear();}
                // Store the Atom in the molecule vector
                molecule.push_back(ATOMS[i]);}
            // Store the last Molecule
            Molecule mol; mol.load(molecule); molecules.push_back(mol);
            // End the PDB method of this function
            }
        else {string temp = "The load function of the C++ DesignGroup class "
            "does not support the " + fileFormat + " file format.";
            cout << temp << endl; exit(EXIT_FAILURE);}
        // Validate the Molecules of the DesignGroup
        validate_Molecules(molecules);
        // End the function
        }

    // Have a function for copying one DesignGroup into another.
    void copy(DesignGroup& group){
        // Copy the number
        number = group.number;
        // Copy the file format
        fileFormat = group.fileFormat;
        // There are currently no file format related attributes
        // copy the force field
        forceField = group.forceField;
        if (forceField == "CHARMM"){
            vdw = group.vdw; elec = group.elec;
            gb = group.gb; lk = group.lk;}
        else {string temp = "The copy function of the C++ DesignGroup class "
            "does not support the " + forceField + " force field.";
            cout << temp << endl; exit(EXIT_FAILURE);}
        // Validate the other group's molecules
        validate_Molecules(group.molecules);
        // Clear the molecules vector
        molecules.clear();
        // Copy and store each individual Molecule
        for(int i=0; i<group.molecules.size(); i++){
            Molecule mol; mol.copy(group.molecules[i]);
            molecules.push_back(mol);}
        // End the function
        }

    // Create a string of text containing the information of all of the
    // DesignGroup's Atoms
    string output(){string text = "";
        for(int i=0; i<molecules.size(); i++){
            for(int j=0; j<molecules[i].residues.size(); j++){
                text += molecules[i].residues[j].output();}}
        return text;}
     
    // Have a function for moving the Design Group
    void move(float coors [], char sign = '-'){
        for(int i=0; i<molecules.size(); i++){molecules[i].move(coors, sign);}}

    // And another for rotating the group
    void rotate(float rmatrix [][3]){
        for(int i=0; i<molecules.size(); i++){molecules[i].rotate(rmatrix);}}

    // And that ends the DesignGroup class
    };

// End the MOLECULES namespace declaration
}
// End the preprocessor commands for including this file in other C++ files
#endif
