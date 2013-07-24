#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>

using namespace std;

int main()
{

string input;

cout << "Please enter the inputfile: ";
cin >> input;

ifstream infile;
ofstream O_outfile, H_outfile;
infile.open(input.c_str());
O_outfile.open("output_O.dat");
H_outfile.open("output_H.dat");
int number_of_atoms;
cout << "Number of atoms: ";
cin >> number_of_atoms;
double x[number_of_atoms],y[number_of_atoms],z[number_of_atoms];
char atom[number_of_atoms];
int i = 0;
while(!infile.eof())
{
         infile >> atom[i] >> x[i] >> y[i] >> z[i];
	 i ++; 
}
for (int i = 0; i < number_of_atoms; i ++)
{
	if (atom[i] == 'O')
	{
		O_outfile << "  " << x[i]/100 << "  " << y[i]/100 << "  " << z[i]/100 << endl;
	}
	else
	{
		H_outfile << "  " <<  x[i]/100 << "  " << y[i]/100 << "  " << z[i]/100 << endl;
	}
}
infile.close();
O_outfile.close();
H_outfile.close();
return 0;
}  





