
/* FFTDock is used for docking two proteins or two find the inital substrate positions
   for OptMAVEN 
*/

#ifndef FFTDOCK_H
#define FFTDOCK_H

#include <iostream>
#include <string>
#include <vector>


using namespace std;
using namespace MOLECULES;
using namespace ENERGY_FUNCTIONS;


class FFTDock {
    
    public:
    
        // Constructor
        FFTDock();

        // Decostructor
        virtual ~FFTDock();

        //Create FFT grid
        void createFFTGrid(Molecule & mol);


    private:
        
};

#endif
