
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "Deimmunization.h"
#include <boost/algorithm/string.hpp>


using namespace std;

int main(int argc, char * argv[]) 
{

    string seqFile = argv[1];
    string spotsFile = argv[2];
    //string databaseName = argv[2];
    Deimmunization de = Deimmunization();
    string databaseName = "Human_9mer_Sequences.txt";
    ifstream in;
    in.open(seqFile.c_str());
    string sequence;
    in >> sequence;
    in.close();
    //cout << sequence << endl;
    in.open(spotsFile.c_str());
    vector<int> positions;
    string line;
    in >> line;
    while(!in.eof())
    {
        boost::trim(line);
        positions.push_back(atoi(line.c_str()));
        in >> line;
    }
    //cout << positions.size() << endl;
    //for (vector<int>::iterator it = positions.begin(); it != positions.end(); it++)
    //{
    //    cout << *it << endl;
    //}
    de.make_humanization_cuts(sequence, positions);

    return 0;
}

