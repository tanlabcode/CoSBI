///////////////////////////////////////////////////////////////////////////////
//
//  Author  : Duygu Ucar
//  Date    : March 2010 
//  
//  Description : Given a co-clustering result
//  It calculates the in-cluster similarity score
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

int main(int argc, char *argv[])
{
        char	filename[100];	
	char	cfilename[100];	
	int C;

	ofstream out;

        if(argc < 3) {
	  cout << endl << "Usage: " << argv[0] << " -fFile [other options]" << endl << endl;
	  cout << "==================== OPTIONS ========================"      << endl << endl;
	  cout << "File Name:           -fString	/*Original 3D Data File name"        << endl;
	  cout << "Clustering File:     -cString	/*Output biclustering result File name"        << endl;
	  cout << "Cluster Number:      -nint	        /*Number of clusters"    << endl;
	  cout << "./qw -fcombined.tricluster -cG250_S3_06_06.tab -n230" <<endl;
	  cout << "=====================================================" << endl << endl;
	  exit(1);
        }
	
        for(int i=1; i<argc; i++){
	  if( argv[i][0] == '-' )
	    switch( argv[i][1] ){
	    case 'f': strcpy( filename, argv[i]+2 );        break;
	    case 'c': strcpy( cfilename, argv[i]+2 );        break;
	    case 'n': //strcpy( C, argv[i]+2 );        break;
	      C = (int) atoi(argv[i]+2); break;
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

#ifdef DEBUG 
	cout <<"Debug Mode is ON."<<endl; 
	array3d.show_size();
#endif

	diff = (clock() - start ) / (double)CLOCKS_PER_SEC;
	//cout<<"Elapsed time for reading in 3d matrices "<<diff <<" seconds\n";

	// Read the clustering file
	Clustering cl(cfilename, C); 
	cl.readTabFile();
	//		cl.show();

	std::string histones[39] = {"CD4-H2AK5ac","CD4-H3K14ac","CD4-H3K9ac","H3K36me1","H3K79me2","H3R2me2","CD4-H2AK9ac","CD4-H3K18ac","CD4-H4K12ac","H2AZ","H3K36me3","H3K79me3","H4K20me1","CD4-H2BK120ac","CD4-H3K23ac","CD4-H4K16ac","H2BK5me1","H3K4me1","H3K9me1","H4K20me3","CD4-H2BK12ac","CD4-H3K27ac","CD4-H4K5ac","H3K27me1","H3K4me2","H3K9me2","H4R3me2","CD4-H2BK20ac","CD4-H3K36ac","CD4-H4K8ac","H3K27me2","H3K4me3","H3K9me3","CD4-H2BK5ac","CD4-H3K4ac","CD4-H4K91ac","H3K27me3","H3K79me1","H3R2me1"};


	// Calculate in-clusters similarity per cluster
	
	start = clock();
	
	double similarity[C];
	for(int i = 0; i < C ; i++){	  
	  similarity[i] = cl.clusters[i].calculate_similarity(array3d);
	  cout <<i<<"\t"<<similarity[i]<<endl;
	}
	
	diff = (clock() - start ) / (double)CLOCKS_PER_SEC;
	//cout<<"Elapsed time for writing clusters to files "<<diff <<" seconds\n"; 

	return 0;

	
}



