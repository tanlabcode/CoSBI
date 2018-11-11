///////////////////////////////////////////////////////////////////////////////
//
//  Author  : Duygu Ucar
//  Date    : Oct - Nov 2009 
//  
//  Description : Given a time series microarray data (M X N X T) - M genes, N conditions, and T time samples,   
//  the algorithm idenfities coherent gene groups over sample sets
//
//  1) Pre-processing: 
//  1-a) The program first calculates a M by M binary matrix for each gene 
//  which indicates if the gene is coherent for the corresponding sample pair, matrix C.
//  1-b) Next, cliques within each of these C matrics are identified. 
//
//  2) 
//
///////////////////////////////////////////////////////////////////////////////

// My header files
#include "input_new.h"

//#include "timeutil.h"
/*#include <stdio.h>
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
#include <ctime>*/
using namespace std;

Array3D	array3d; // holds the 3d data n a single matrix of this order (G,S,T)
int G, S, AT;
float T, C;

void print_C_array(int C[], int size){
  for(int i = 0; i< size ; i++){
    for(int j = 0; j< size ; j++){
      cout<<C[i*size+j]<<",";
    }
    cout<<endl;
  }
}



/*-------------------------------------------------------------*/
//                    -> QHu starts here <-                    //
void searchClique(int* cMatrix, vector<int> head, vector<int> tail, vector< vector<int> >* maxCliq, int size, int sizeThr);
void recurSearch(vector<int> head, vector<int> tail, vector< pair< vector<int>, vector<int> > >* gsSet, int minS, int minG, vector< vector<pair<int,int> > >* invList, vector<pair<int,int> > ILintsec);

vector< vector<int> > find_cliques (int* cMatrix, int size, int sizeThr)
{
    vector< vector<int> > maxCliq;
    for(int n=0; n<=size-sizeThr; n++)
    {
        vector<int> head (1, n);
        
        //get tail
        vector<int> tail;
        for(int m = n+1; m < size; m++)
        {
            tail.push_back(m);
        }
        
        searchClique(cMatrix, head, tail, &maxCliq, size, sizeThr);
    }
    return maxCliq;
}

void searchClique(int* cMatrix, vector<int> head, vector<int> tail, vector< vector<int> >* maxCliq, int size, int sizeThr)
{
            
    //remove elements from tail according to the c matrix
    int i = head[head.size()-1];
    for(int n = 0; n < tail.size(); n++)
    {
        int j = tail[n];
        if (cMatrix[i*size+j] == 0)
        {
            tail.erase(tail.begin()+n);
            n--;    //the index should move one step back after an element is deleted
        }
    }
         
    //check the total length
    if(head.size()+tail.size() < sizeThr)
        return;
    
    //check the situation of union
    vector<int> container;
    for(int n = 0; n < (*maxCliq).size(); n++)
    {
        for(int m = 0; m < ((*maxCliq)[n]).size(); m++)
        {
            container.push_back(((*maxCliq)[n]).at(m));
        }
    }
    sort (container.begin(), container.end());
    
    vector<int> continent;
    for(int n = 0; n < head.size(); n++)
    {
        continent.push_back(head[n]);
    }
    for(int n = 0; n < tail.size(); n++)
    {
        continent.push_back(tail[n]);
    }
    sort (continent.begin(), continent.end());
    
    if (includes(container.begin(), container.end(), continent.begin(), continent.end()) )
        return;

    //the situation of |tail|=0
    if(tail.size() == 0)
    {
        (*maxCliq).push_back(head);
    }
    else
    {
        while(tail.size() != 0)
        {
            vector<int> headTemp = head;
            headTemp.push_back(tail[0]);
            tail.erase(tail.begin());
            searchClique(cMatrix, headTemp, tail, maxCliq, size, sizeThr);
        }
    }
    return;
}

vector< pair< vector<int>, vector<int> > > SampleGeneSearch(int size, int minS, int minG,vector< vector<pair<int,int> > >* invList)
{
    vector< pair< vector<int>, vector<int> > > gsSet;
    for (int n=0; n<=size-minS; n++)
    {
        vector<int> head (1, n);

        //get tail
        vector<int> tail;
        for (int m = n+1; m < size; m++)
        {
            if ((*invList)[m].size() >= minG)

                tail.push_back(m);
        }
        recurSearch(head, tail, &gsSet, minS, minG, invList, (*invList)[n]);
    }
    return gsSet;
}
    
void recurSearch(vector<int> head, vector<int> tail, vector< pair< vector<int>, vector<int> > >* gsSet, int minS, int minG, vector< vector<pair<int,int> > >* invList, vector<pair<int,int> > ILintsec)
{
    vector<pair<int,int> > g ((*invList)[head[0]]);
    //check the g size
    //for(int n=1; n<head.size(); n++)
    {
        //get the intersection
        int n = head.size()-1;
        //vector<pair<int,int> > tempG (g);
        //g.resize(min(tempG.size(), (*invList)[head[n]].size()));
        g.resize(set_intersection(ILintsec.begin(), ILintsec.end(), (*invList)[head[n]].begin(), (*invList)[head[n]].end(), g.begin())-g.begin());
    }
    vector<int> gSet;
    for(int n=0; n<g.size(); n++)
    {
        gSet.push_back(g[n].first);
    }
    gSet.resize(unique(gSet.begin(), gSet.end())-gSet.begin());

    if(gSet.size()<minG)
        return;

    //remove samples as described in 4.1.2
    /*for(int n = 0; n < tail.size(); n++)
    {
        if((*invList)[tail[n]].size() < minG)
        {
            tail.erase(tail.begin()+n);
            n--;
        
    }*/
    //check the total length
    if(head.size()+tail.size() < minS)
        return;
        
    while(tail.size() != 0)
    {
        vector<int> headTemp = head;
        headTemp.push_back(tail[0]);
        tail.erase(tail.begin());
        recurSearch(headTemp, tail, gsSet, minS, minG, invList, g);
    }
    
    if(head.size()>=minS)
    {
        for(int n=0; n<(*gsSet).size(); n++)
        {
            if((*gsSet)[n].first.size()==gSet.size() && (*gsSet)[n].second.size()>=head.size())
            {
                //if(includes((*gsSet)[n].first.begin(), (*gsSet)[n].first.end(), gSet.begin(), gSet.end()))
                {
                    if(includes((*gsSet)[n].second.begin(), (*gsSet)[n].second.end(), head.begin(), head.end()))
                    {
                        return;
                    }
                }
            }
            else if((*gsSet)[n].first.size()==gSet.size() && (*gsSet)[n].second.size()<head.size())
            {
                //if(includes(gSet.begin(), gSet.end(), (*gsSet)[n].first.begin(), (*gsSet)[n].first.end()))
                {
                    if(includes(head.begin(), head.end(), (*gsSet)[n].second.begin(), (*gsSet)[n].second.end()))
                    {
                        (*gsSet).erase((*gsSet).begin()+n);
                        n--;
                    }
                }
            }
        }
        
	pair< vector<int>, vector<int> > gs (gSet, head);
	// if any of the sample pairs is similar enough
	if (AT ==1){
	  for(int i = 1; i < gSet.size(); i++)
	    {
	      for(int k = 0; k < i; k++)            
		{
		  for(int j = 0; j < head.size(); j++)
		    {
		      if(cross_cor(array3d.get3d(gSet[i],head[j]), array3d.get3d(gSet[k],head[j]), array3d.T, array3d.T) >= C)
		      //if(pearsoncorrelation(array3d.get3d(gSet[i],head[j]), array3d.get3d(gSet[k],head[j]), array3d.T) >= C)
                        goto GoOn;
		    }
		  return;
                GoOn:
		  ;
		}
	    }
	  //  pair< vector<int>, vector<int> > gs (gSet, head);
	  (*gsSet).push_back(gs);
	}

	// if all sample pairs are similar enough
	else if(AT ==2){
	  for(int i = 1; i < gSet.size(); i++)
	    {
	      for(int k = 0; k < i; k++)            
		{
		  for(int j = 0; j < head.size(); j++)
		    {
		      if(cross_cor(array3d.get3d(gSet[i],head[j]), array3d.get3d(gSet[k],head[j]), array3d.T, array3d.T) < C)
			//if(pearsoncorrelation(array3d.get3d(gSet[i],head[j]), array3d.get3d(gSet[k],head[j]), array3d.T) >= C)
                        goto GoOn2;
		    }
		}
	    }
	  //pair< vector<int>, vector<int> > gs (gSet, head);
	  (*gsSet).push_back(gs);
	GoOn2:
	  ;
	}

	else {
	  cout << "Invalid option for analysis type. Cannot set to "<<AT<<endl;
	}

        //pair< vector<int>, vector<int> > gs (gSet, head);
        
        /*for(int i=0; i<gs.first.size(); i++)
        {
            cout<<gs.first.at(i)<<' ';
        }
        cout<<endl;
        for(int i=0; i<gs.second.size(); i++)
        {
            cout<<gs.second.at(i)<<' ';
        }
        cout<<endl;*/
        
        //(*gsSet).push_back(gs);
    }
}

//                    -> QHu ends here <-                      //
/*-------------------------------------------------------------*/



int main(int argc, char *argv[])
{
        char	filename[100];		
	ofstream out;

        if(argc < 2) {
	  cout << endl << "Usage: " << argv[0] << " -fFile [other options]" << endl << endl;
	  cout << "==================== OPTIONS =====================" << endl << endl;
	  cout << "File Name:           -fString	/*File name"         << endl;
	  cout << "Gene Count:           -g1000	        /*Gene set size"         << endl;
	  cout << "Sample Count:           -s3	        /*Sample set size"         << endl;
	  cout << "Within sim thr:         -t0.7	/*Similarity threshold"         << endl;
	  cout << "Between sim thr:        -c0.5	/*Similarity threshold"         << endl;
	  cout << "Type:                   -k1	        /* 1-min, 2-max, 3-average"         << endl;
	  cout << "==================== OPTIONS =====================" << endl << endl;
	  cout << "Eg, ./cocluster+ -fcombined.tricluster -g1000 -s3 -t0.7 -c0.5 -k1" << endl;
	  cout << "==================================================" << endl << endl;
	  exit(1);
        }
	
        for(int i=1; i<argc; i++){
	  if( argv[i][0] == '-' )
	    switch( argv[i][1] ){
	    case 'f': strcpy( filename, argv[i]+2 );        break;
	    case 'g':
	      G = (int) atoi(argv[i]+2); break;
	    case 's':
	      S = (int) atoi(argv[i]+2); break;
	    case 't':
	      T = (float) atof(argv[i]+2); break;
	    case 'c':
	      C = (float) atof(argv[i]+2); break;
	    case 'k':
	      AT = (int) atoi(argv[i]+2); break;
	    default: cout << "Error Input Options: " << argv[i] << endl; exit(0); break;
	    }
	  else{	cout << "Error Input Options: " << argv[i] << endl;  exit(0); }
        }

	cout <<"Analyzing "<< filename<<"..."<<endl;
	// start clock
	time_t start,end;
	double dif;

	time (&start);
	// Define variables
	//Array3D	array3d; // holds the 3d data n a single matrix of this order (G,S,T)
	Input	in(filename, array3d); 
	in.readTabFile(); // my data is loaded into the array3d matrix

        float threshold = T; // to convert the correlation scores into a binary value

	int** matrix; // This will hold the Cross Correlation matrices for each gene
	int numCols = array3d.G;
	int numRows = array3d.S*array3d.S;

	cout<<"Number of Genes:"<<array3d.G<<endl;
	cout<<"Number of Samples:"<<array3d.S<<endl;
	cout<<"Length of Signal:"<<array3d.T<<endl;
	//array3d.show();

	/* Initialize  2d array */
	/* Allocate pointer memory for the first dimension of the matrix[][]; */
	matrix = (int **) malloc(numCols * sizeof(int *));
	if(NULL == matrix){free(matrix); printf("Memory allocation failed while allocating for matrix[].\n"); exit(-1);}
	
	/* Allocate integer memory for the second dimension of the matrix[][]; */
	for(int x = 0; x < numCols; x++){
	  matrix[x] = (int *) malloc(numRows * sizeof(int));
	  if(NULL == matrix[x]){free(matrix[x]); printf("Memory allocation failed while allocating for matrix[x][].\n"); exit(-1);}
	}
	
	for(int y = 0; y < numRows; y++)
	  for(int x = 0; x < numCols; x++)
	    matrix[x][y] = 0;
	
	
	// For each gene construct a S x S matrix
	for(int i = 0; i < (array3d.G); i++){
#ifdef DEBUG 
	  cout <<"Gene"<<i<<endl; 
#endif
	  for(int j = 0; j <array3d.S ; j++){
	    for(int k = 0 ; k <array3d.S; k++){

#ifdef DEBUG
	      cout<<j<<" "<<k<<endl;
#endif
	      

	      float temp = cross_cor( array3d.get3d(i,j), array3d.get3d(i,k), array3d.T, array3d.T);
	      //float temp = pearsoncorrelation(array3d.get3d(i,j), array3d.get3d(i,k), array3d.T);
	      if (temp > threshold){
		matrix[i][(j*(array3d.S)+k)] = 1;
	      }
	    }
	  }
	}
	
	time (&end);
	dif = (double)(difftime (end,start)/ (double) 60);
	printf ("Elapsed time for constructing CC matrices with cross correlation: %.2lf minutes.\n", dif );


    //--QHu--//

	time(&start);
    
    vector< vector<int> > maxCliq[(const int)array3d.G];
    for(int i = 0; i < (array3d.G); i++)
    {
        //print_C_array(matrix[i],array3d.S);
        //find_cliques(C);
        maxCliq[i] = find_cliques(matrix[i], array3d.S, S);
        //Display each set in a line
        //cout<<"g"<<i<<": ";
        for(int n = 0; n < maxCliq[i].size(); n++)
        {
            for(int m = 0; m < maxCliq[i].at(n).size(); m++)
            {
                int s = (maxCliq[i]).at(n).at(m);
                //cout<<'s'<<s<<',';
            }
            //cout<<" / ";
        }
        //cout<<endl;
    }
    
    //cout<<endl;
    vector< vector<pair<int,int> > > invList;
    for(int eachS = 0; eachS < 39; eachS++)
    {
        //cout<<'s'<<eachS<<": ";
        vector<pair<int,int> > temp;
        for(int i = 0; i < (array3d.G); i++)
        {
            for(int n = 0; n < maxCliq[i].size(); n++)
            {
                if (binary_search ((maxCliq[i]).at(n).begin(), (maxCliq[i]).at(n).end(), eachS))
                {
                    //cout<<'g'<<i<<",b"<<n<<", ";
                    temp.push_back(make_pair (i,n));
                }
            }
        }
        //cout<<endl;
        invList.push_back(temp);
    }    

    time (&end);
    dif = (double)(difftime (end,start)/ (double) 60);
    printf ("Elapsed time for finding cliques: %.2lf minutes.\n", dif );


    time (&start);
    //SampleGeneSearch(39,3,1000,&invList);
    vector< pair< vector<int>, vector<int> > > gsSet = SampleGeneSearch(39,S,G,&invList);

    //Print results to a file
    //Clustering(vector< pair< vector<int>, vector<int> > > gs, int c)
      
    Clustering cl(gsSet, gsSet.size());   
    
    char* str = "BiCluster_new"; 
    char result[80];	
	
    cout<<gsSet.size()<<" biclusters are identified"<<endl;

    for(int i = 0; i < gsSet.size(); i++){
      sprintf( result, "%s%d", str, i );
      cl.clusters[i].writeFile(array3d,result);
      }


     
       for (int n=0; n<gsSet.size(); n++)
	{

        for(int i=0; i<gsSet[n].first.size(); i++)
        {
	  cout<<gsSet[n].first.at(i);
	  cout<<' ';
        }
        cout<<endl;
        for(int i=0; i<gsSet[n].second.size(); i++)
        {
            cout<<gsSet[n].second.at(i)<<' ';
        }
        cout<<endl;
	
	}
        

    //--QHu--//

    time (&end);
    dif = (double)(difftime (end,start)/ (double) 60);
    printf ("Elapsed time for generating GS sets: %.2lf minutes.\n", dif );
    
    return 0;
}


