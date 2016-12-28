
#ifndef INITIALANTIGEN_H
#define INITIALANTIGEN_H


#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <utility>
#include "MOLECULES.h"
//#include "ENERGY_FUNCTIONS.h"
//#include "DOCKING_FUNCTIONS.h"

using namespace std;
using namespace MOLECULES;
//using namespace ENERGY_FUNCTIONS;
//using namespace DOCKING_FUNCTIONS;
class InitialAntigen {
    
    public:
        
        //Constructor
        InitialAntigen();
        
        //Deconstructor
        virtual ~InitialAntigen();

        // Int to string
        string intToString(int x);
        // Load parameters
        map<string, float> load_parameters(string parameterFile);
        
        // Load antigen
        Molecule load_antigen(string antigenFile);
        
        // Load parts
        map<string, Molecule> load_parts(string domainFile);
        
        // Load translation_xyz.txt
        vector< vector<float> > load_translations(string translationFile);
        
        // Load rotation_angle.txt
        vector<float> load_rotations(string rotationFile);
        
        // Generate random float number with mean and std
        float generate_random_number(float a, float b);
        float generate_gauss_random_number(float mean, float std);
        // Generate rotation matrix
        void calculate_rmatrix(float angle, char direction, float rmatrix [3][3]);
        
        // Calculate the coordinates center of the selected molecule
        void calculate_molecule_center(Molecule molecule, float average [3]);
        void calculate_residue_center(vector<Residue> residues, float average [3]);
        
        // Move the center of the molecules to the orgin
        Molecule move_molecule_orgin(Molecule molecule);
        Molecule move_molecule_orgin(Molecule molecule, float average [3]);
        
        // Calculate the distance between two atoms
        float calculate_atom_distance(Atom atom1, Atom atom2);
        
        // Count the total molecule atoms
        int count_molecule_atom(Molecule molecule);
        
        // Measure the overlap between two molecules
        float measure_molecule_overlap(Molecule receptor, Molecule ligand);
        
        // Filter the antigens and remove those with high overlapped with antibody parts
        vector<Molecule> molecule_overlap_filter(vector<Molecule> antibodyParts, Molecule antigen, float maxOverlap);
        
        // Move the antigen so that the epitope residues have the most negative Z coordinates in the antigen
        Molecule move_antigen_epitope_znegative(Molecule antigen, vector<int> epitopeResiduesNames);
        
        // General the antigen initial position relative to the antibody
        vector<Molecule> initial_antigen_position(Molecule antigen, vector<int> epitopeResiduesNames);
        
        // General the antigen initial position relative to the antibody
        vector<Molecule> initial_antigen_position(Molecule antigen, vector< vector<float> > translations, vector<float> rotations);
                 
        // Format the atom line to pdb format atom line
        string format_pdb_line(string atomLine, int atomNum);
        // Ouput the molecule to pdb file
        void output_molecule_pdb(Molecule molecule, string pdbName);
        // Calculate the interaction energy between two molecules
        float calculate_interaction_energy(Molecule part, Molecule antigen, bool solvation);
        //float calculate_interaction_energy(vector<Molecule> parts, vector<Molecule> antigens);
        // Output the energies to MAPs_Energies.txt
        void output_maps_energies(map<string, Molecule> partsMap, Molecule antigen, string outFile, bool solvation);
        // Perturb the antigen 
        void antigens_perturbation(vector<Molecule> antigens, float cartesian_step, float angle_step);

        //Random antigen movement
        Molecule antigen_random_move(Molecule antigen, float cartesianStep, float angleStep);
        
        // Output the optimized antibody part conbination  after selecting the best antibody parts based on the interaction energy between antibody and antigen using python cplex script
        float output_optimize_part(string mapEnergyFile, string resultsFile);  
        //Refinment the antigen
        vector<string> antigen_position_refinement(Molecule & antigen, string & domainsName, string & summary);
        Molecule antigen_position_refinement(Molecule & antigen, Molecule & antibodyH, Molecule & antibodyL);

       // vector<Molecule> antigens_refinement(map<string, float> parameters, vector<Molecule> parts, vector<Molecule> antigens);


    private:


};

#endif
