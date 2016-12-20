#include "molvec.cpp"
#include "energy.cpp"
using namespace std;

int main(int argc, char* argv[])
{
	MolVec mv1, mv2;
	mv1.load(argv[1]);
	mv2.load(argv[2]);
	cout << interactionEnergy(mv1, mv2);
}
