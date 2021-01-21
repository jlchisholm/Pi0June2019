// temporary to get photon energy column for tagger conversion

#include "iostream"
#include "fstream"
using namespace std;

void MakeFile()
{
    ifstream oldfile ("Tagg_Con2.txt");
    string line;

    ofstream newfile;
    newfile.open("new_Tagg_con.txt");

    double ch;
    double en;

    int i = 0;
    while (!oldfile.eof())
    {
        oldfile >> ch >> en;

        newfile << ch << "\t" << 450 - en << endl;

        ++i;
    }

    newfile.close();
    oldfile.close();
}
