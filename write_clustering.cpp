///////////////////////////////////////////////////////////////////////////////
//
//  Author  : Duygu Ucar
//  Date    : March 2010 
//  
//  Description : Given a co-clustering result
//  It generates a file that contains signals for the members of the cluster
//  It generates a file that can be read by matlap function for plotting 
//
//
///////////////////////////////////////////////////////////////////////////////

// My header files
#include "input_new.h"

// Standard header files
#include <stdio.h>
#include <unistd.h> 
#include <fstream>
#include <assert.h>
#include <iostream>
#include <math.h>
#include <cstdlib>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <algorithm>
#include <set>                                                                                 
#include <vector>
#include <map>
#include <ctime>
#include <string>
using namespace std;

//#define NAME_LEN	50
#define debug(...)   fprintf(stderr, __VA_ARGS__)
//example usage -> debug("Delay %d\n",delay);

#define DEBUG
#undef DEBUG

 void writeNames(Array3D& array3d, char* file_name)
    {
      char file[100], *p;
      TabInput in2;
      in2.open( file_name );

      for (int k=0; k<array3d.S; k++)
	{
	  
	  in2.getLine();
	  p = in2.get2tab();
	  strcpy(array3d.sName[k].str, p);
	}

      
      for (int k=0; k<array3d.G; k++)
	{
	  in2.getLine();
	  p = in2.get2tab();
	  strcpy(array3d.gName[k].str, p);
	}

    }

int main(int argc, char *argv[])
{
        char	filename[100];	
	char	cfilename[100];
	char	gfilename[100];
	int C;

	ofstream out;

        if(argc < 3) {
	  cout << endl << "Usage: " << argv[0] << " -fFile [other options]" << endl << endl;
	  cout << "==================== OPTIONS ========================"      << endl << endl;
	  cout << "File Name:           -fString	/*Original 3D Data File name"        << endl;
	  cout << "Clustering File:     -cString	/*Output biclustering result File name"        << endl;
	  cout << "Names File:          -gString        /*File that contains gene names and sample names"<<endl;
	  cout << "Cluster Number:      -nint	        /*Number of clusters"    << endl;
	  cout << "./qw -fcombined.tricluster -ccombined_1000_3.tab -n907" <<endl;
	  cout << "=====================================================" << endl << endl;
	  exit(1);
        }
	
        for(int i=1; i<argc; i++){
	  if( argv[i][0] == '-' )
	    switch( argv[i][1] ){
	    case 'f': strcpy( filename, argv[i]+2 );        break;
	    case 'c': strcpy( cfilename, argv[i]+2 );        break;
	    case 'g': strcpy( gfilename, argv[i]+2 );        break;
	    case 'n': C = (int) atoi(argv[i]+2); break;
	    default: cout << "Error Input Options: " << argv[i] << endl; exit(0); break;
	    }
	  else{	cout << "Error Input Options: " << argv[i] << endl;  exit(0); }
        }

#ifdef DEBUG 
	cout <<"number of clusters:"<< C <<endl;
#endif

	// start clock
	clock_t start;
	double diff;
	start = clock();

	// Define variables
	Array3D	array3d; // holds the 3d data n a single matrix of this order (G,S,T)
	Input	in(filename, array3d); 
	in.readTabFile(); // my data is loaded into the array3d matrix

	//writeNames(array3d,gfilename);


#ifdef DEBUG 
	cout <<"Debug Mode is ON."<<endl; 
	array3d.show_size();
#endif

	diff = (clock() - start ) / (double)CLOCKS_PER_SEC;
	//	cout<<"Elapsed time for reading in 3d matrices "<<diff <<" seconds\n";

	// Read the clustering file
	Clustering cl(cfilename, C); 
	cl.readTabFile();
	cl.show_all(array3d);

	std::string histones[39] = {"CD4-H2AK5ac","CD4-H3K14ac","CD4-H3K9ac","H3K36me1","H3K79me2","H3R2me2","CD4-H2AK9ac","CD4-H3K18ac","CD4-H4K12ac","H2AZ","H3K36me3","H3K79me3","H4K20me1","CD4-H2BK120ac","CD4-H3K23ac","CD4-H4K16ac","H2BK5me1","H3K4me1","H3K9me1","H4K20me3","CD4-H2BK12ac","CD4-H3K27ac","CD4-H4K5ac","H3K27me1","H3K4me2","H3K9me2","H4R3me2","CD4-H2BK20ac","CD4-H3K36ac","CD4-H4K8ac","H3K27me2","H3K4me3","H3K9me3","CD4-H2BK5ac","CD4-H3K4ac","CD4-H4K91ac","H3K27me3","H3K79me1","H3R2me1"};


	// Write clusters into a file
	
		start = clock();
		/*		char* str = "BiCluster"; 
	char result[80];	
	
	for(int i = 0; i < C ; i++){
	  sprintf( result, "%s%d", str, i );
	  cl.clusters[i].writeFile(array3d,result);
	  }
	
	  diff = (clock() - start ) / (double)CLOCKS_PER_SEC;*/
	//cout<<"Elapsed time for writing clusters to files "<<diff <<" seconds\n"; 

	return 0;

}


