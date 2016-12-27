
#include <fstream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include <initializer_list>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <random>
#include "InitialAntigen.h"
//#include "MOLECULES.h"
#include "DOCKING_FUNCTIONS.h"
#include "ENERGY_FUNCTIONS.h"
//#include <boost/accumulators/framework/accumulator_set.hpp>
//#include <boost/accumulators/statistics/mean.hpp>
//#include <boost/accumulators/statistics/variance.hpp>
//using namespace boost::accumulators;
//using namespace MOLECULES;
//using namespace DOCKING_FUNCTIONS;
//Constructor
using namespace std;

InitialAntigen::InitialAntigen() 
{
}

//Deconstructor
InitialAntigen::~InitialAntigen() {}

using namespace std;

string InitialAntigen::intToString(int x)
{
    /**************************************/
    /* This function is similar to itoa() */
    /* "integer to alpha", a non-standard */
    /* C language function. It takes an   */
    /* integer as input and as output,    */
    /* returns a C++ string.              */
    /* itoa()  returned a C-string (null- */
    /* terminated)                        */
    /* This function is not needed because*/
    /* the following template function    */
    /* does it all                        */
    /**************************************/   
    string r;
    stringstream s;
    s << x;
    r = s.str();
    return r;

}

map<string, float> InitialAntigen::load_parameters(string parameterFile)
{
    map<string, float> parameters;
    ifstream in(parameterFile.c_str());
    string line;
    getline(in, line);
    while(!in.eof())
    {
        vector<string> strs;
        boost::split(strs, line, boost::is_any_of(":"));
        parameters[strs[0]] = atof(strs[1].c_str());

    }

    return parameters;

}

Molecule InitialAntigen::load_antigen(string moleculeFile)
{
    ifstream in1(moleculeFile.c_str());
    vector<Atom> atoms;
    string line;
    getline(in1, line);
    while(!in1.eof())
    {
        Atom atom;
        atom.load(line);
        atoms.push_back(atom);
        getline(in1, line);
    }
    in1.close();
    Molecule molecule;
    molecule.load(atoms);
    return molecule;
}

map<string, Molecule> InitialAntigen::load_parts(string domainFile)
{
    map<string, Molecule> parts;
    vector<string> kinds;
    vector<Atom> atoms;
    ifstream in1(domainFile.c_str());
    string line;
    in1 >> line;
    while(!in1.eof())
    {
        kinds.push_back(line); 
        in1 >> line;
    }

    in1.close();
    // The different maps parts
    string regions [3] = {"V", "CDR3", "J"};
    // Go through the different region / part combinations
    for(int i=0; i<kinds.size(); i++)
    {
        for(int j=0; j<3; j++)
        {
            // Make the name of the expected part and file
            string partName = kinds[i] + regions[j];
            string fileName = partName + ".txt";
            // Open the expected file
            ifstream in2(fileName.c_str());
            // Clear the Atoms vector
            atoms.clear();
            // Get the information about the first Atom in the file
            getline(in2, line);
            Atom a; 
            a.load(line);
            int oldN = a.rotamerNumber;
            // Read in the Atoms in the file
            while(!in2.eof())
            {
                Atom atom; 
                atom.load(line);
                // Determine if this Atom is the first in a new part or not
                if(atom.rotamerNumber != oldN)
                {
                    // Store the part
                    Molecule part; 
                    part.load(atoms); 
                    string name = partName + "_" + to_string(oldN); 
                    parts[name] = part;
                    //parts.push_back(part);
                    // Clear the Atom list
                    atoms.clear();
                    // Store the new part number
                    oldN = atom.rotamerNumber;
                }
                // Store the atom
                atoms.push_back(atom);
                // Get the next line
                getline(in2, line);
             }
            // Close the file
            in2.close();
            // Store the last part
            Molecule last; 
            last.load(atoms); 
            string name = partName + "_" + to_string(oldN); 
            parts[name] = last;
            //parts.push_back(last);
        }
    }
    return parts;
}

vector< vector<float> > InitialAntigen::load_translations(string translationFile)
{
    vector< vector<float> > translations;
    ifstream in1(translationFile.c_str());
    string line;
    in1 >> line;
    vector<string> strs;
    vector<float> trans;
    while(!in1.eof())
    {
        strs.clear();
        trans.clear();
        boost::trim(line);
        boost::split(strs, line, boost::is_any_of(":"));
        //boost::split(strs, line, boost::is_space());
        //cout << strs.size() << endl;
        for(vector<string>::iterator it = strs.begin(); it != strs.end(); it++)
        {
            //cout << atof((*it).c_str()) << endl;
            trans.push_back(atof((*it).c_str()));
        }
        //cout << trans.size() << endl;
        translations.push_back(trans);
        in1 >> line;
    }
    in1.close();

    return translations;

}

vector<float> InitialAntigen::load_rotations(string rotationFile)
{
    vector<float> rotations;
    ifstream in1(rotationFile.c_str());
    string line;
    in1 >> line;
    while(!in1.eof())
    {
        //float rot = atof(boost::trim(line).c_str());
        float rot = atof(line.c_str());
        rotations.push_back(rot);
        in1 >> line;
    }
    in1.close();

    return rotations;
}

float InitialAntigen::generate_random_number(float a, float b)
{
    //srand(time(NULL));
    //return ((b-a)*((float)rand()/RAND_MAX))+a;
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(a, b);
    return dis(gen);

}

float InitialAntigen::generate_gauss_random_number(float mean, float std)
{
    // Generate random float number with mean and std
    // The source code needs to be compiled with the option -std=c++0x
    std::normal_distribution<float> distribution(mean, std);
    std::default_random_engine generator(time(0));
    float number = distribution(generator);
    cout << number << endl; 
    return number;
}

void InitialAntigen::calculate_rmatrix(float angle, char direction, float rmatrix [3][3])
{
    float c = cos(angle);
    float s = sin(angle);
    float v = 1 - c;
    // Get the proper unit vector
    float r [3] = {0, 0, 0};
    if (direction == 'x')
    {
        r [0] = 1.0;
    }
    else if(direction == 'y')
    {
        r[1] = 1.0;
    }
    else if(direction == 'z')
    {
        r[2] = 1.0;
    }
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

void InitialAntigen::calculate_molecule_center(Molecule molecule, float average [3])
{
    average[0] = 0.0;
    average[1] = 0.0;
    average[2] = 0.0;
    int count = 0;
    for(int residueIndex=0; residueIndex < molecule.residues.size(); residueIndex++)
    {
        for(int atomIndex=0; atomIndex < molecule.residues[residueIndex].atoms.size(); atomIndex++)
        {
            Atom atom = molecule.residues[residueIndex].atoms[atomIndex];
            count++;
            average[0] += atom.x;
            average[1] += atom.y;
            average[2] += atom.z;
            
        }
    }

    for(int j=0; j < 3; j++)
        {
            if (count == 0)
            {
                exit(EXIT_FAILURE);
            }
            else 
            {
                average[j] /= count;
            }
        }
    //cout << average[0] << " " << average[1] << " " << average[2] << endl;
}

void InitialAntigen::calculate_residue_center(vector<Residue> residues, float average [3])
{
    average[0] = 0.0;
    average[1] = 0.0;
    average[2] = 0.0;
    int count = 0;
    for(int residueIndex=0; residueIndex < residues.size(); residueIndex++)
    {
        for(int atomIndex=0; atomIndex < residues[residueIndex].atoms.size(); atomIndex++)
        {
            Atom atom = residues[residueIndex].atoms[atomIndex];
            count++;
            average[0] += atom.x;
            average[1] += atom.y;
            average[2] += atom.z;
            
        }
    }

    for(int j=0; j < 3; j++)
        {
            if (count == 0)
            {
                exit(EXIT_FAILURE);
            }
            else 
            {
                average[j] /= count;
            }
        }

    //cout << average[0] << " " << average[1] << " " << average[2] << endl;
}
Molecule InitialAntigen::move_molecule_orgin(Molecule molecule)
{
    float average [3] ;
    calculate_molecule_center(molecule, average);
    //cout << average[0] << " " << average[1] << " " << average[2] << endl;
    for(int residueIndex=0; residueIndex < molecule.residues.size(); residueIndex++)
    {
        for(int atomIndex=0; atomIndex < molecule.residues[residueIndex].atoms.size(); atomIndex++)
        {
            
            //cout << molecule.residues[residueIndex].atoms[atomIndex].x << endl;
            molecule.residues[residueIndex].atoms[atomIndex].move(average, '-');
            //cout << molecule.residues[residueIndex].atoms[atomIndex].x << endl;
        }
    }
    return molecule;
}


Molecule InitialAntigen::move_molecule_orgin(Molecule molecule, float average [3])
{
    for(int residueIndex=0; residueIndex < molecule.residues.size(); residueIndex++)
    {
        for(int atomIndex=0; atomIndex < molecule.residues[residueIndex].atoms.size(); atomIndex++)
        {
            molecule.residues[residueIndex].atoms[atomIndex].move(average, '-');
        }
    }
    return molecule;
}

int InitialAntigen::count_molecule_atom(Molecule molecule)
{
    int count = 0;
    for(int residueIndex=0; residueIndex < molecule.residues.size(); residueIndex++)
    {
        for(int atomIndex=0; atomIndex < molecule.residues[residueIndex].atoms.size(); atomIndex++)
        {
            count++;
        }
    }
    return count;

}

float InitialAntigen::measure_molecule_overlap(Molecule receptor, Molecule ligand)
{
    int length = count_molecule_atom(ligand);
    float overlop = 0.0;
}

Molecule InitialAntigen::move_antigen_epitope_znegative(Molecule antigen, vector<int> epitopeResiduesNames)
{
    // Move the antigen to the orgin
    antigen = move_molecule_orgin(antigen);
    float pi = 3.14;
    // Get the epitope residue
    vector<Residue> epitope;
    for(vector<int>::iterator it = epitopeResiduesNames.begin(); it != epitopeResiduesNames.end(); it++)
    {
        epitope.push_back(antigen.residues[*it]);
    }
    // Get the z coordinates for the epitope residues
    vector<double> caZ_best;
    for(int residueIndex=0; residueIndex < epitope.size(); residueIndex++)
    {
        for(int atomIndex=0; atomIndex < epitope[residueIndex].atoms.size(); atomIndex++)
        {
            //cout << "epitope residue number:  " << epitope[residueIndex].number << endl;
            if (epitope[residueIndex].atoms[atomIndex].name == "CA")
            {
                caZ_best.push_back(epitope[residueIndex].atoms[atomIndex].z);
            }
        }
    }
    // Get the summary, mean and std of the z coordinates from the inital antigen conformation as the best 
    double bestZ = accumulate(caZ_best.begin(), caZ_best.end(), 0.0);
    double zmean_best = bestZ / caZ_best.size();
    double accum_best = 0.0;
    accum_best = std::inner_product(caZ_best.begin(), caZ_best.end(), caZ_best.begin(), 0.0);
    double zstd_best = sqrt(accum_best / (caZ_best.size()-1));
    
    cout << "First " << "bestZ: " << bestZ << "zstd_best: " << zstd_best << endl;
    Molecule antigen_best;
    antigen_best.copy(antigen);
    int totalDegree = 360;
    int stepDegree = 3;
    int xstep = totalDegree / stepDegree;
    int ystep = totalDegree / stepDegree;
    int zstep = totalDegree / stepDegree;
    
    for (int x = 1; x < xstep; x++)
    {
        for (int y = 1; y < ystep; y++)
        {
            float anglex = (x * stepDegree) * pi / 180.0;
            float angley = (y * stepDegree) * pi / 180.0;
            cout << "x: " << x << " anglex: " << anglex << endl;
            float xrmatrix [3][3];
            calculate_rmatrix(anglex, 'x', xrmatrix);
            float yrmatrix [3][3];
            calculate_rmatrix(angley, 'y', yrmatrix);
            Molecule antigen_new;
            antigen_new.copy(antigen);
            antigen_new.rotate(xrmatrix);
            antigen_new.rotate(yrmatrix);
            vector<Residue> epitope_new;
            for(vector<int>::iterator it = epitopeResiduesNames.begin(); it != epitopeResiduesNames.end(); it++)
            {
                epitope_new.push_back(antigen_new.residues[*it]);
            }

            vector<double> caZ;
            for(int residueIndex=0; residueIndex < epitope_new.size(); residueIndex++)
            {
                for(int atomIndex=0; atomIndex < epitope_new[residueIndex].atoms.size(); atomIndex++)
                {
                    if (epitope_new[residueIndex].atoms[atomIndex].name == "CA")
                    {
                        caZ.push_back(epitope_new[residueIndex].atoms[atomIndex].z);
                    }
                }
            }
            double sumZ = accumulate(caZ.begin(), caZ.end(), 0.0);
            double zmean = sumZ / caZ.size();
            double accum = 0.0;
            accum= std::inner_product(caZ.begin(), caZ.end(), caZ.begin(), 0.0);
            double zstd = sqrt(accum / (caZ.size()-1));
            cout << "sumZ: " << sumZ << " zstd: " << zstd << endl;
            if (zstd < zstd_best)
            {
        //cout << "bestZ: " << bestZ << "sumZ: " << sumZ << endl;
                bestZ = sumZ;
                zmean_best = zmean;
                zstd_best = zstd;
                antigen_best.copy(antigen_new);
            }
        }
    }
    
    antigen.copy(antigen_best);
    
    for (int x = 1; x < xstep; x++)
    {
        float anglex = (x * stepDegree) * pi / 180.0;
        cout << "x: " << x << " anglex: " << anglex << endl;
        float xrmatrix [3][3];
        calculate_rmatrix(anglex, 'x', xrmatrix);
        Molecule antigen_new;
        antigen_new.copy(antigen);
        antigen_new.rotate(xrmatrix);
        vector<Residue> epitope_new;
        for(vector<int>::iterator it = epitopeResiduesNames.begin(); it != epitopeResiduesNames.end(); it++)
        {
            epitope_new.push_back(antigen_new.residues[*it]);
        }

        vector<double> caZ;
        for(int residueIndex=0; residueIndex < epitope_new.size(); residueIndex++)
        {
            for(int atomIndex=0; atomIndex < epitope_new[residueIndex].atoms.size(); atomIndex++)
            {
                if (epitope_new[residueIndex].atoms[atomIndex].name == "CA")
                {
                    caZ.push_back(epitope_new[residueIndex].atoms[atomIndex].z);
                }
            }
        }
        double sumZ = accumulate(caZ.begin(), caZ.end(), 0.0);
        double zmean = sumZ / caZ.size();
        double accum = 0.0;
        accum= std::inner_product(caZ.begin(), caZ.end(), caZ.begin(), 0.0);
        double zstd = sqrt(accum / (caZ.size()-1));
        cout << "sumZ: " << sumZ << " zstd: " << zstd << endl;
        if (sumZ < bestZ)
        {
            //cout << "bestZ: " << bestZ << "sumZ: " << sumZ << endl;
            bestZ = sumZ;
            zmean_best = zmean;
            zstd_best = zstd;
            antigen_best.copy(antigen_new);
        }
   } 
    antigen.copy(antigen_best);
    for (int y = 1; y < ystep; y++)
    {
        float angley = (y * stepDegree) * pi / 180.0;
        cout << "y: " << y << " angley: " << angley << endl;
        float yrmatrix [3][3];
        calculate_rmatrix(angley, 'y', yrmatrix);
        Molecule antigen_new;
        antigen_new.copy(antigen);
        antigen_new.rotate(yrmatrix);
        vector<Residue> epitope_new;
        for(vector<int>::iterator it = epitopeResiduesNames.begin(); it != epitopeResiduesNames.end(); it++)
        {
            epitope_new.push_back(antigen_new.residues[*it]);
        }

        vector<double> caZ;
        for(int residueIndex=0; residueIndex < epitope_new.size(); residueIndex++)
        {
            for(int atomIndex=0; atomIndex < epitope_new[residueIndex].atoms.size(); atomIndex++)
            {
                if (epitope_new[residueIndex].atoms[atomIndex].name == "CA")
                {
                    caZ.push_back(epitope_new[residueIndex].atoms[atomIndex].z);
                }
            }
        }
        double sumZ = accumulate(caZ.begin(), caZ.end(), 0.0);
        double zmean = sumZ / caZ.size();
        double accum = 0.0;
        accum= std::inner_product(caZ.begin(), caZ.end(), caZ.begin(), 0.0);
        double zstd = sqrt(accum / (caZ.size()-1));
        cout << "sumZ: " << sumZ << " zstd: " << zstd << endl;
        if (sumZ < bestZ)
        {
            //cout << "bestZ: " << bestZ << "sumZ: " << sumZ << endl;
            bestZ = sumZ;
            zmean_best = zmean;
            zstd_best = zstd;
            antigen_best.copy(antigen_new);
        }
    }

    cout << "bestZ: " << bestZ << " zstd_best: " << zstd_best << endl; 
            
    return antigen_best;
}    
    /*int iterations = 200;
    float pi = 3.142;
    //float startT = 1;
    //float endT = 0.01;
    //float delta = (endT - startT) / ( iterations -1);
    vector<Residue> epitope;
    antigen = move_molecule_orgin(antigen);
    for(vector<int>::iterator it = epitopeResiduesNames.begin(); it != epitopeResiduesNames.end(); it++)
    {
        epitope.push_back(antigen.residues[*it]);
    }

    double bestZ = 0.0;
    vector<double> caZ_best;
    //for(vector<Residue>::iteractor it = epitope.begin(); it != epitope.end(); it++)
    
    for(int residueIndex=0; residueIndex < epitope.size(); residueIndex++)
    {
        for(int atomIndex=0; atomIndex < epitope[residueIndex].atoms.size(); atomIndex++)
        {
            
            cout << "epitope residue number:  " << epitope[residueIndex].number << endl;
            //bestZ += epitope[residueIndex].atoms[atomIndex].z;
            if (epitope[residueIndex].atoms[atomIndex].name == "CA")
            {
                caZ_best.push_back(epitope[residueIndex].atoms[atomIndex].z);
            }
        }
    }
       
    bestZ = accumulate(caZ_best.begin(), caZ_best.end(), 0.0);
    double zmean_best = bestZ / caZ_best.size();
    double accum_best = 0.0;
    //std::for_each (std::begin(caZ_best), std::end(caZ_best), [&](const double d) 
    //{
    //        accum_best += (d - z_mean_best) * (d - zmean_best);
    //});
    accum_best = std::inner_product(caZ_best.begin(), caZ_best.end(), caZ_best.begin(), 0.0);
    double zstd_best = sqrt(accum_best / (caZ_best.size()-1));

    cout << "First " << bestZ << endl;
    Molecule antigen_best;
    antigen_best.copy(antigen);
    
    for(int i=0; i < iterations ; i++)
    {
       float anglex = generate_random_number(0.0, 10.0) * pi / 180.0;
       sleep(1);
       float angley = generate_random_number(0.0, 10.0) * pi / 180.0;
       sleep(1);
       float anglez = generate_random_number(0.0, 10.0) * pi / 180.0;
       //sleep(1);
       cout << anglex << " " << angley << " " << anglez << endl;
       float xrmatrix [3][3];
       calculate_rmatrix(anglex, 'x', xrmatrix);
       float yrmatrix [3][3];
       calculate_rmatrix(angley, 'y', yrmatrix);
       float zrmatrix [3][3];
       calculate_rmatrix(anglez, 'z', zrmatrix);
       Molecule antigen_new;
       antigen_new.copy(antigen_best);
       antigen_new.rotate(xrmatrix);
       antigen_new.rotate(yrmatrix);
       antigen_new.rotate(zrmatrix);
       vector<Residue> epitope_new;
    
       for(vector<int>::iterator it = epitopeResiduesNames.begin(); it != epitopeResiduesNames.end(); it++)
       {
           epitope_new.push_back(antigen_new.residues[*it]);
       }

       double sumZ = 0.0;
       vector<double> caZ;
       double zmean;
       double zstd; 
      //for(vector<Residue>::iteractor it = epitope.begin(); it != epitope.end(); it++)
       for(int residueIndex=0; residueIndex < epitope_new.size(); residueIndex++)
       {
           for(int atomIndex=0; atomIndex < epitope_new[residueIndex].atoms.size(); atomIndex++)
           {
               if (epitope_new[residueIndex].atoms[atomIndex].name == "CA")
               {
                   caZ.push_back(epitope_new[residueIndex].atoms[atomIndex].z);
                   //sumZ += epitope_new[residueIndex].atoms[atomIndex].z;
               }
           }
       }
    sumZ = accumulate(caZ.begin(), caZ.end(), 0.0);
    zmean = sumZ / caZ.size();
    double accum = 0.0;
    //std::for_each (std::begin(caZ), std::end(caZ), [&](const double d) 
    //{
    //        accum += (d - z_mean_best) * (d - zmean_best);
    //});
    accum= std::inner_product(caZ.begin(), caZ.end(), caZ.begin(), 0.0);
    zstd = sqrt(accum / (caZ.size()-1));

       cout << "bestZ: " << bestZ << " zmean_best: " << zmean_best << " zstd_best: " << zstd_best << endl;
       cout << "sumZ: " << sumZ << " zmean: " << zmean << " zstd: " << zstd << endl;

        //if (( sumZ < bestZ) && ( zmean < zmean_best) && (zstd < zstd_best) )
        if (zstd < zstd_best)
        {
            //cout << "bestZ: " << bestZ << "sumZ: " << sumZ << endl;
            bestZ = sumZ;
            zmean_best = zmean;
            zstd_best = zstd;
            antigen_best.copy(antigen_new);
        }
        else
        {
            continue;
            //float t = startT - i * delta;
            //float t = 300.0;
            //float diff = sumZ - bestZ;
            //float rnd = generate_random_number(0, 1.0);
            //cout << "rnd: " << rnd << endl;
            cout << "-diff/t: " << -diff/t << " exp: " << exp(-diff/t) << endl; 
            double boltz = round(-diff/t);
            float probability = exp(min(40.0,  max(-40.0, boltz)));  
            //if (exp(-diff/t) > rnd)
            if (probability > rnd)
            {
                bestZ = sumZ;
                antigen_best.copy(antigen_new);
            }
        }
    }*/


vector<Molecule> InitialAntigen::initial_antigen_position(Molecule antigen, vector<int> epitopeResiduesNames)
{
    //antigen = move_molecule_orgin(antigen);
    //string text = antigen.output();
    //outfile << text << endl;
    Molecule antigen_new = move_antigen_epitope_znegative(antigen, epitopeResiduesNames);
    output_molecule_pdb(antigen_new, "test.pdb");
    vector<Molecule> antigens;
    // Rotate 60 degree around z axis 
    float pi = 3.142;
    float angle = 0.0;
    float angleBegin = 0.0;
    float rotationAngle = 60.0;
    float radian = pi * rotationAngle / 180;

    for(int step = 1; step <= 6; step++)
    {
        angle = angleBegin + step * radian; 
        float zrmatrix [3][3];
        calculate_rmatrix(angle, 'z', zrmatrix);
        Molecule antigen_rotation;
        antigen_rotation.copy(antigen_new);
        antigen_rotation.rotate(zrmatrix);
    
        vector<Residue> epitope;
        for(vector<int>::iterator it = epitopeResiduesNames.begin(); it != epitopeResiduesNames.end(); it++)
        {
           epitope.push_back(antigen_rotation.residues[*it]);
        }
        
        float average [3] ;
        calculate_residue_center(epitope, average); 
        Molecule antigen_center = move_molecule_orgin(antigen_rotation, average);
     
        // Translate the antigen along z, y, x axis
        // z from 3.75 to 16.25  1.25 interval
        // y from -5 to 10   2.50 interval
        // x from -10 to 5   2.50 interval 
        float zmin = -16.25;
        float zmax = -3.75;
        float zinterval = 1.25;
        int zstep = (zmax - zmin) / zinterval;
        float ymin = -5;
        float ymax = 10;
        float yinterval = 2.50;
        int ystep = (ymax - ymin) / yinterval;
        float xmin = -10;
        float xmax = 5;
        float xinterval = 2.50;
        int xstep = (xmax - xmin) / xinterval;
        float zmove = 0.0;
        float ymove = 0.0;
        float xmove = 0.0;
        for(int z = 0; z <= zstep ; z++)
        {
            zmove = zmin + z * zinterval; 
            for(int y = 0; y <= ystep ; y++)
            {
                ymove = ymin + y * yinterval;
                for(int x = 0; x <= xstep ; x++)
                {
                    xmove = xmin + x * xinterval;
                    float coor [3] = {float(xmove), float(ymove), float(zmove)};
                    Molecule antigen_current;
                    antigen_current.copy(antigen_center);
                    antigen_current.move(coor, '-');
                    antigens.push_back(antigen_current);
                }
            }
        }
    }
    return antigens;
}

vector<Molecule> InitialAntigen::initial_antigen_position(Molecule antigen, vector< vector<float> > translations, vector<float> rotations)
{
    vector<Molecule> antigens;
    for(vector< vector<float> >::iterator it = translations.begin(); it != translations.end(); it++)
    {
        // Do the translations
        float coor [3] = {(*it)[0], (*it)[1], (*it)[2]};
        Molecule antigen_current;
        antigen_current.copy(antigen);
        antigen_current.move(coor, '-');
        // Rotate 60 degree around z axis 
        for(vector<float>::iterator rt = rotations.begin(); rt != rotations.end(); rt++)
        {           
            float zrmatrix [3][3];
            calculate_rmatrix(*rt, 'z', zrmatrix);
            antigen_current.rotate(zrmatrix);
            antigens.push_back(antigen_current);
        }
    }

    return antigens;
}
    
string InitialAntigen::format_pdb_line(string atomLine, int atomNum)
{
    vector<string> strs;
    boost::split(strs, atomLine, boost::is_any_of("\t "));
    //const char * atom_format = "ATOM  %5d %4s%1c%3s %1c%4d%1c   %8.3f%8.3f%8.3f%6.2f%6.2f            %2s";
    
    char line[81];
    char chainName = strs[0].c_str()[0];
    char alt=' ';
    int residueNum = atoi(strs[1].c_str());
    const char * residueName = strs[2].c_str();
    const char * atomName = strs[4].c_str();
    float x = atof(strs[5].c_str());
    float y = atof(strs[6].c_str());
    float z = atof(strs[7].c_str());
    float m = 1.01;
    float b = 0.01;
    //cout << atomNum << " " << atomName << " " << residueName << " " << chainName << " " << residueNum << " " << x << " " << y << " " << z << " " << m << " " << b << endl;
    sprintf(line, "ATOM  %5d %4s%1c%3s %1c%4d%1c   %8.3f%8.3f%8.3f%6.2f%6.2f", atomNum, atomName, alt, residueName, chainName, residueNum, alt, x, y, z, m, b);
    string pdbLine(line);
    return pdbLine;

}
void InitialAntigen::output_molecule_pdb(Molecule molecule, string pdbName)
{
    int count = 0;
    string pdb;
    for(int residueIndex=0; residueIndex < molecule.residues.size(); residueIndex++)
    {
        for(int atomIndex=0; atomIndex < molecule.residues[residueIndex].atoms.size(); atomIndex++)
        {
            count++;
            string text = molecule.residues[residueIndex].atoms[atomIndex].output();
            string pdbLine = format_pdb_line(text, count);
            pdb += pdbLine + "\n";
        }
    }
    ofstream out1(pdbName.c_str());
    out1 << pdb;
}

float InitialAntigen::calculate_interaction_energy(Molecule part, Molecule antigen, bool solvation)
{
    float energy;
    // The energies between this structure and the antigen
    float energies [3] = {0, 0, 0};
    // Go through the Residues and Atoms of the part
    for(int a=0; a<part.residues.size(); a++)
    {
        for(int b = 0; b<part.residues[a].atoms.size(); b++)
        {
            Atom atom1 = part.residues[a].atoms[b];
            // Go through the Residues and Atoms of the antigen
                for(int Z=0; Z<antigen.residues.size(); Z++)
                {
                    for(int z=0; z<antigen.residues[Z].atoms.size(); z++)
                    {
                        Atom atom2 = antigen.residues[Z].atoms[z];
                        //int skip = ENERGY_FUNCTIONS::exclusions(atom1, atom2);
                        ENERGY_FUNCTIONS::IE_calculation(atom1, atom2, energies, 0, solvation);
                    }
                }
        }
    }
    // Calculate the total energy for this MAPs structure
    energy = energies[0] + energies[1] + energies[2];
    return energy;
}

void InitialAntigen::output_maps_energies(map<string, Molecule> partsMap, Molecule antigen, string outFile, bool solvation)
{
    //ofstream of1("MAPs_Energies.txt");
    ofstream of1;
    of1.open(outFile.c_str(), ios::app);
    
    for(map<string, Molecule>::iterator it = partsMap.begin(); it != partsMap.end(); it++)
    {
        string partName = it->first;
        Molecule part = it->second;
        float energy = calculate_interaction_energy(part, antigen, solvation);
        
        //cout << energy << endl;
        of1 << partName << "  " << energy << "\n";

    }

    of1.close();

}

void InitialAntigen::antigens_perturbation(vector<Molecule> antigens, float cartesianStep, float angleStep)
{
    float pi = 3.14;
    float coordinates [3];
    // Generate random number of cartesian and angle movement
    coordinates[0] = generate_random_number(0.0, cartesianStep);
    coordinates[1] = generate_random_number(0.0, cartesianStep);
    coordinates[2] = generate_random_number(0.0, cartesianStep);
    float anglex = generate_random_number(0.0, angleStep) * pi / 180.0;
    float angley = generate_random_number(0.0, angleStep) * pi / 180.0;
    float anglez = generate_random_number(0.0, angleStep) * pi / 180.0;

    // Cartesian perturbation
    for (int molIndex = 0; molIndex < antigens.size(); molIndex++ )
    {
        antigens[molIndex].move(coordinates);
    }

    // Rotation

    float xrmatrix [3][3];
    calculate_rmatrix(anglex, 'x', xrmatrix);
    float yrmatrix [3][3];
    calculate_rmatrix(angley, 'y', yrmatrix);
    float zrmatrix [3][3];
    calculate_rmatrix(anglez, 'z', zrmatrix);

    for (int molIndex = 0; molIndex < antigens.size(); molIndex++ )
    {
        antigens[molIndex].rotate(xrmatrix);
        antigens[molIndex].rotate(yrmatrix);
        antigens[molIndex].rotate(zrmatrix);
    }

}

Molecule InitialAntigen::antigen_random_move(Molecule antigen, float cartesianStep, float angleStep)
{
    float pi = 3.14;
    float coordinates [3];
    // Generate random number of cartesian and angle movement
    //coordinates[0] = generate_gauss_random_number(0.0, cartesianStep);
    //coordinates[1] = generate_gauss_random_number(0.0, cartesianStep);
    //coordinates[2] = generate_gauss_random_number(0.0, cartesianStep);
    //float anglex = generate_gauss_random_number(0.0, angleStep) * pi / 180.0;
    //float angley = generate_gauss_random_number(0.0, angleStep) * pi / 180.0;
    //float anglez = generate_gauss_random_number(0.0, angleStep) * pi / 180.0;
    coordinates[0] = generate_random_number(-cartesianStep, cartesianStep);
    coordinates[1] = generate_random_number(-cartesianStep, cartesianStep);
    coordinates[2] = generate_random_number(-cartesianStep, cartesianStep);
    float anglex = generate_random_number(-angleStep, angleStep) * pi / 180.0;
    float angley = generate_random_number(-angleStep, angleStep) * pi / 180.0;
    float anglez = generate_random_number(-angleStep, angleStep) * pi / 180.0;
    cout << coordinates[0] << " " << coordinates[1] << " " << coordinates[2] << endl;
    cout << anglex << " " << angley << " " << anglez << endl;
    // Cartesian movement 
    antigen.move(coordinates);
    // Rotation movement
    float xrmatrix [3][3];
    calculate_rmatrix(anglex, 'x', xrmatrix);
    float yrmatrix [3][3];
    calculate_rmatrix(angley, 'y', yrmatrix);
    float zrmatrix [3][3];
    calculate_rmatrix(anglez, 'z', zrmatrix);
    antigen.rotate(xrmatrix);
    antigen.rotate(yrmatrix);
    antigen.rotate(zrmatrix);

    return antigen;
}

float InitialAntigen::output_optimize_part(string mapEnergyFile, string resultsFile)
{
   cout << mapEnergyFile << endl;
   cout << resultsFile << endl;
   char command[100];
   sprintf(command, "python output_optimize_parts.py %s %s", mapEnergyFile.c_str(), resultsFile.c_str());
   /*for(int i= 0; i < strlen(command); i++)
   {
       cout << command[i] << endl; 
   }*/
   system(command);

}


Molecule InitialAntigen::antigen_position_refinement(Molecule & antigen, Molecule & antibodyH, Molecule & antibodyL)
{
    
    float startT = 3640; 
    float endT = 2190; 
    int iterations = 500;
    float delta = (startT - endT)/(iterations - 1);
    float cartesianStep = 1.0; 
    float angleStep = 5.0; 
    float energyBest = 0.0;
    float gasConstant = 0.001986;
    Molecule antigenBest;
    antigenBest.copy(antigen);
    for(int i=0; i < int(iterations); i++)
    {
        // Do the random movement of the antigen 
        Molecule antigenCurrent = antigen_random_move(antigenBest, cartesianStep, angleStep);
        bool solvation = false;
        float energyCurrent = calculate_interaction_energy(antigenCurrent, antibodyH, solvation) + calculate_interaction_energy(antigenCurrent, antibodyL, solvation);
        cout << "energyBest: " << energyBest << " energyCurrent: " << energyCurrent << endl;
        if (energyCurrent < energyBest)
        {
            // Save the currentEnergy as bestEnergy so far
            energyBest = energyCurrent;
            antigenBest = antigenCurrent;
            cout << "Iteration " << i + 1 << endl;
            cout << "Best energy: " << energyBest << endl;

        } else
        {
            // Simulated annealing parameters
            float t = startT - i * delta;
            float diff = energyCurrent - energyBest;
            float rnd = generate_random_number(0.0, 1.0);
            cout << "t: " << t << " diff: " << diff << " rnd: " << rnd << endl; 
            // Simulated annealing to determinze whether to keep this perturbation
            if (exp(-diff / gasConstant * t) > rnd)
            {
            // Save the currentEnergy as bestEnergy so far
                energyBest = energyCurrent;
                antigenBest = antigenCurrent;
                cout << "Iteration " << i + 1 << endl;
                cout << "Best energy: " << energyBest << endl;
            }
        }
     }
     // Return the refined antigens
     return antigenBest;
}

vector<string> InitialAntigen::antigen_position_refinement(Molecule &antigen, string &domainsName, string &summary)
{
    map<string, Molecule> parts = load_parts(domainsName);
    float startT = 3640; 
    float endT = 2190; 
    int iterations = 500;
    float delta = (startT - endT)/(iterations - 1);
    float cartesianStep = 1.0; 
    float angleStep = 5; 
    float energyBest = 0.0;
    float gasConstant = 0.001986;
    Molecule antigenBest;
    antigenBest.copy(antigen);
    vector<string> partsBest;
    for(int i=0; i < int(iterations); i++)
    {
        // Do the random movement of the antigen 
        Molecule antigenCurrent = antigen_random_move(antigenBest, cartesianStep, angleStep);
        // Calculate and output the interaction between the antigne and antibody parts
        string energyFile = "antigen_" + to_string(i) + "_maps_energies.txt";
        bool solvation = false;
        output_maps_energies(parts, antigenCurrent, energyFile, solvation);
        // Select the optimized part using python script
        // Output the optimized energy and optimized parts to results.txt
        string results = "results.txt";
        output_optimize_part(energyFile, results);
        ifstream inf(results.c_str());
        string line;
        inf >> line;
        vector<string> temp;
        
        while(!inf.eof())
        {
            boost::trim(line);
            temp.push_back(line);
            inf >> line;
        }
        float energyCurrent = atof(temp[0].c_str());
        cout << "energyCurrent: " << energyCurrent << " energyBest: " << energyBest << endl;
        //string summary = "summary.txt";
        ofstream outfile(summary.c_str(), std::ios_base::app);
        //If the new energy is better than the best energy so far, keep it
        if (energyCurrent < energyBest)
        {
            // Save the currentEnergy as bestEnergy so far
            energyBest = energyCurrent;
            antigenBest = antigenCurrent;
            partsBest.clear();
            for (vector<string>::iterator it = temp.begin() + 1; it != temp.end(); it++)
            {
                partsBest.push_back(*it);
            }
            outfile << "Iteration " << i + 1 << endl;
            outfile << "Best energy: " << energyBest << endl;
            outfile << "Best parts: ";
            for (vector<string>::iterator it = partsBest.begin(); it != partsBest.end(); it++)
            {
                outfile << *it << " ";
            }
            outfile << endl << endl;

        } else
        {
            // Simulated annealing parameters
            float t = startT - i * delta;
            float diff = energyCurrent - energyBest;
            float rnd = generate_random_number(0.0, 1.0);
            cout << "t: " << t << " diff: " << diff << " rnd: " << rnd << endl; 
            // Simulated annealing to determinze whether to keep this perturbation
            if (exp(-diff / gasConstant * t) > rnd)
            {
            // Save the currentEnergy as bestEnergy so far
                energyBest = energyCurrent;
                antigenBest = antigenCurrent;
                partsBest.clear();
                for (vector<string>::iterator it = temp.begin() + 1; it != temp.end(); it++)
                {
                    partsBest.push_back(*it);
                }
                outfile << "Iteration " << i + 1 << endl;
                outfile << "Best energy: " << energyBest << endl;
                outfile << "Best parts: ";
                for (vector<string>::iterator it = partsBest.begin() + 1; it != partsBest.end(); it++)
                {
                    outfile << *it << " ";
                }   
                outfile << endl << endl;
            }
        }
     }
     // Return the refined antigens
     antigen = antigenBest;
     return partsBest;

}


