#ifndef INPUT_NEW_H
#define INPUT_NEW_H

#include <time.h>
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

using namespace std;

#define NAME_LEN	200
#define LENGTH          25
#define debug(...)   fprintf(stderr, __VA_ARGS__)
//example usage -> debug("Delay %d\n",delay);

#define DEBUG
#undef DEBUG

/* Given two sequences calculates cross correlation *
 * Returns the max of all pairs                     */
float cross_cor(float* x, float* y, int n, int maxdelay){

#ifdef DEBUG
  cout<<"first series: ";
  
  for (int k = 0; k < n; k++)
    {
      cout<<x[k]<<",";
    }
  cout<<endl;
  
  cout<<"second series: ";
  
  
  for (int k = 0; k < n; k++)
    {
      cout<<y[k]<<",";
    }
  cout<<endl;
#endif    

  int i,j;
  float mx,my,sx,sy,sxy,denom,r;
  float max = 0;
  int nzeros1 = 0;
  int nzeros2 = 0;
 
  /* Calculate the mean of the two series x[], y[] */
  mx = 0;
  my = 0;   
  for (i=0;i<n;i++) {
    mx += x[i];
    my += y[i];
    if(x[i] != 0){
        nzeros1++;
    }
    if(y[i] != 0){
        nzeros2++;
    }
  }
  mx /= n;
  my /= n;

  /* Return zero if there are more than 10 zeros in any sequence*/
  if (nzeros1 <= 5 || nzeros2 <= 5 ){
    return 0;
  }
 
  /* Calculate the denominator */
  sx = 0;
  sy = 0;
  for (i=0;i<n;i++) {    sx += (x[i] - mx) * (x[i] - mx);
    sy += (y[i] - my) * (y[i] - my);
  }
  denom = sqrt(sx*sy);

  /* Return zero if denominator is zero */
  if (denom == 0){
    return 0;
  }

  /* Calculate the correlation series */
  for (int delay=-maxdelay;delay<maxdelay;delay++) {
    sxy = 0;
    for (i=0;i<n;i++) {
      j = i + delay;

      //Ignore cases where there are no points
      if (j < 0 || j >= n)
	continue;
      else
	sxy += (x[i] - mx) * (y[j] - my);
      
      // Replace with zeros
      /* if (j < 0 || j >= n)
	 sxy += (x[i] - mx) * (-my);
	 else
	 sxy += (x[i] - mx) * (y[j] - my);
      */
    }

    r = sxy / denom;

#ifdef DEBUG
    cout<<"At delay:"<<delay<<" corr "<<r<<endl;
#endif
    
    if(r>max)
      max = r;
    /* r is the correlation coefficient at "delay" */
  }

#ifdef DEBUG
  cout <<"MAX : "<<max<<endl;  
#endif
  
     return max;
  
}

/* calculates Pearson Correlation Coefficient for a given pair */
/* returns the absolute value of the coefficient */
float pearsoncorrelation(float* x, float* y, int size)
{

    float xmean = 0;
    float ymean = 0;
    float result = 0;
    int nzeros1 = 0;
    int nzeros2 = 0;

    // calculate mean
    for (int i = 0; i < size; i++)
    {
        xmean = xmean+x[i];
        ymean = ymean+y[i];

	if(x[i] != 0){
	  nzeros1++;
	}
	if(y[i] != 0){
	  nzeros2++;
	}

    }
    xmean = xmean/size;
    ymean = ymean/size;

    if (nzeros1 <= 5 || nzeros2 <= 5 ){
      return 0;
    }

    // calculate denominator
    double s = 0;
    double xv = 0;
    double yv = 0;
    double t1;
    double t2;
    for (int i = 0; i < size; i++)
    {
        t1 = x[i]-xmean;
        t2 = y[i]-ymean;
        xv = xv+t1*t1;
        yv = yv+t2*t2;
        s = s+t1*t2;
    }
    if ( xv==0||yv==0 )
    {
        result = 0;
    }
    else
    {
        result = s/(sqrt(xv)*sqrt(yv));
    }
    //  return fabs(result);//absolute value of PC
    return result;
}

struct Name
{
    char str[NAME_LEN];
};

struct Profile
{
    float p[LENGTH];
};

struct Plot  //QHu
{
    vector<int> csample;
    vector< vector<int> > submatrix;
};

class Array3D
{
    float	*data;

public:

    Name	*tName, *sName, *gName;
    int	T, S, G;

    void init( int t, int s, int g )
    {
        T=t;
        G=g;
        S=s;
        data = new float[T*S*G];
        assert( data != NULL );
        tName = new Name[T];
        assert( tName != NULL );
        sName = new Name[S];
        assert( sName != NULL );
        gName = new Name[G];
        assert( gName != NULL );
    }

    /*        float& dat( int t, int s, int g ){
                  return data[t*(S*G)+s*G+g];
          }

      float* get2d( int t ){
                  return data + t*(S*G);
      }


          float& operator()( int t, int s, int g ){
                  return data[t*(S*G)+s*G+g];
          }
    */

    float& dat( int g, int s, int t)
    {
        return data[g*(T*S)+s*T+t];
    }

    float* get2d( int g )
    {
        return data + g*(T*S);
    }

    float* get3d( int g, int s)
    {
        return data + g*(T*S) + s*T;
    }


    float& operator()( int g, int s, int t )
    {
        return data[g*(T*S)+s*T+t];
    }


    void clear()
    {
        delete   data;
        delete[] tName;
        delete[] sName;
        delete[] gName;
    }

    void show_size()
    {
        cout << endl;
        cout << "Total Times:\t"   << T << endl;
        cout << "Total Samples:\t" << S << endl;
        cout << "Total Genes:\t"   << G << endl;
    }

    void show()
    {
        cout << endl;
        cout << "Total Times:\t"   << T << endl;
        cout << "Total Samples:\t" << S << endl;
        cout << "Total Genes:\t"   << G << endl;

        for (int i=0; i<T; i++)
        {
            cout << "Time\t" << i << endl;
            cout << "ID\tNAME";
            for (int k=0; k<S; k++)
                cout << "\tS-" << k;
            cout << endl;

            for (int j=0; j<G; j++)
            {
                cout << j << "\tG-" << j;
                for (int k=0; k<S; k++)
                    printf( "\t%2.2f", dat(j,k,i) );
                cout << endl;
            }
        }
        cout << endl;
    }
};



class TabInput
{
#define MAX_LINE	100000
#define MAX_WORD	1000
    char	buf[MAX_LINE], sub[MAX_WORD];
    int	len, k;
    ifstream infile;

public:
    void open(char *file)
    {
        infile.open(file);
    }

    void close()
    {
        infile.close();
    }

    void getLine()
    {
        infile.getline(buf, MAX_LINE, '\n');
        len = strlen( buf );
        assert(len <= MAX_LINE);
        k=0;
    }

    char* get2tab()
    {
        int p=0;

        if (buf[k] == '\t' || buf[k]=='\n')
        {
            sub[0] = '\0';
            if (buf[k] == '\t')
                k++;
        }
        else
        {
            do
            {
                sub[p++] = buf[k++];
            }
            while ( buf[k]!='\t' && buf[k]!=0xd && k<len ); //oxd for ^M
            k++;
            sub[p] = '\0';
            //assert(p <= MAX_WORD);
            assert(p <= NAME_LEN);
        }
        return (char *)sub;
    }
};


class GeneCluster
{
    //  #define LENGTH	25
    float* centroid;
    vector<Profile> genes;
    vector<int> ids;

public:

    GeneCluster()
    {
        centroid = new float[LENGTH];
        for (int i = 0 ; i < LENGTH ; i++)
        {
            centroid[i] = 0;
        }
    }

    GeneCluster(float* g)
    {
        centroid = new float[LENGTH];
        for (int i = 0 ; i < LENGTH ; i++)
        {
            centroid[i] = g[i];
        }
    }

    float* getCentroid()
    {
        return centroid;
    }

    float getWithinDissimilarity()
    {
        //cout<<"Within dis is ";
        float wd = 0;
        float pair_count = 0;
        for (int i = 0; i < (genes.size()-1); i++)
        {
            for (int k = i; k < genes.size(); k++)
            {

                Profile pr = genes.at(i);
                Profile pr2 = genes.at(k);

                float dis = pearsoncorrelation(pr.p, pr2.p, LENGTH);
                wd = wd + dis;
                pair_count++;
            }
        }

        wd = 1-(wd/pair_count);
        //cout <<wd<<endl;
        return wd;
    }

    float getWithinDissimilarity(int g)
    {
        float dis = 0;
        for (int i = 0; i < (genes.size()); i++)
        {
            if (i != g)
            {
                Profile pr = genes.at(i);
                Profile pr2 = genes.at(g);
                dis = dis + pearsoncorrelation(pr.p, pr2.p, LENGTH);
            }
        }
        return (dis/(genes.size()-1));
    }

    float getBetweenDissimilarity(GeneCluster* cluster2)
    {
        //cout<<"Between dis is ";
        float bd = 0;
        float pair_count = 0;
        for (int i = 0; i < genes.size(); i++)
        {
            for (int k = 0; k < cluster2->genes.size(); k++)
            {

                Profile pr = genes.at(i);
                Profile pr2 = cluster2->genes.at(k);

                float dis = pearsoncorrelation(pr.p, pr2.p, LENGTH);
                bd = bd + dis;
                pair_count++;
            }
        }

        bd = 1-(bd/pair_count);
        //cout <<bd<<endl;
        return bd;
    }


    float calculateSSE()
    {
        float dis = 0;
        for (int i = 0; i < genes.size(); i++)
        {
            Profile pr = genes.at(i);
            float temp = 1 - pearsoncorrelation(pr.p, centroid, LENGTH);
            dis = dis + temp*temp;
        }
        return dis;
    }

    void setCentroid(float* temp)
    {
        for (int i = 0 ; i < LENGTH ; i++)
        {
            centroid[i] = temp[i];
        }
    }

    void insertGene(Profile g, int gid)
    {
        genes.push_back(g);
        ids.push_back(gid);
    }

    void removeElements()
    {
        while (!genes.empty())
        {
            genes.pop_back();
        }
        while (!ids.empty())
        {
            ids.pop_back();
        }
    }

    void recalculateCentroid()
    {
        float center[LENGTH];
        for (int i = 0 ; i < LENGTH ; i++)
        {
            center[i] = 0;
        }

        for (int i = 0; i < genes.size(); i++)
        {
            Profile pr = genes.at(i);
            for (int j = 0 ; j < LENGTH ; j++)
            {
                center[j] = center[j] + pr.p[j];
            }
        }

        for (int i = 0 ; i < LENGTH ; i++)
        {
            centroid[i] = center[i]/genes.size();
        }

    }

    vector<int> getIds()
    {
        return ids;
    }

    void displayGenes()
    {
        cout <<genes.size()<< " genes available."<<endl;
        for (int i = 0; i < genes.size(); i++)
        {
            Profile pr = genes.at(i);
            for (int i = 0 ; i < LENGTH ; i++)
            {
                cout << pr.p[i] << " ";
            }
            cout<<endl;
        }
    }

    void displaySize()
    {
        cout <<genes.size()<< " genes available."<<endl;
        for (int i = 0; i < ids.size(); i++)
        {
            cout<<ids.at(i)<<" ";
        }
        cout<<endl;
    }


    void displayCentroid()
    {
        cout<<"Cluster Size:"<<genes.size()<<endl;
        cout<<"Centroid elements: ";
        for (int i = 0 ; i < LENGTH ; i++)
        {
            cout<<centroid[i]<<" ";
        }
        cout<<endl;
    }

};


class Cluster
{
    vector <int> cgene;
    vector <int> csample;
    friend class Array3D;

public:

    int	G, S;

    Cluster()
    {
        G = 0;
        S = 0;
    }

    //QHu
    Cluster(vector <int> gene, vector <int> sample)
    {
        cgene = gene;
        csample = sample;
        G = gene.size();
        S = sample.size();
    }

    vector<int> get_genes()
    {
        return cgene;
    }

    vector<int> get_samples()
    {
        return csample;
    }

    void init( int s, int g )
    {
        G=g;
        S=s;
    }

    void add_gene(int a)
    {
        cgene.push_back(a);
        G++;
    }

    void add_sample (int a)
    {
        csample.push_back(a);
        S++;
    }

    float* get_mean_profile(Array3D& array3d, int id, float *result)
    {
        for (int j=0; j<array3d.T; j++)
        {
            result[j]  = 0;
        }

        //cout<<S<<" samples"<<endl;
        for (int i =0; i< S; i++)
        {
            float * submatrix = array3d.get3d(cgene[id], csample[i]);
            for (int j=0; j<array3d.T; j++)
            {
                result[j]  = result[j] + submatrix[j];
                //cout<<submatrix[j]<<" ";
            }
            //cout<<endl;
        }

        for (int j=0; j<array3d.T; j++)
        {
            result[j]  = result[j]/S;
            //      cout<<result[j]<<" ";
        }
        return result;
    }

    float* get_mean_profile(Array3D& array3d, int id)
    {
      float* result;
      for (int j=0; j<array3d.T; j++)
        {
	  result[j]  = 0;
        }
      
      //cout<<S<<" samples"<<endl;
      for (int i =0; i< S; i++)
        {
	  float * submatrix = array3d.get3d(cgene[id], csample[i]);
	  for (int j=0; j<array3d.T; j++)
            {
	      result[j]  = result[j] + submatrix[j];
	      //cout<<submatrix[j]<<" ";
            }
	  //cout<<endl;
        }
      
      for (int j=0; j<array3d.T; j++)
        {
	  result[j]  = result[j]/S;
            //      cout<<result[j]<<" ";
        }
      return result;
    }
    
    

    void show_summary()
    {
        cout << endl;
        cout << "Total Genes:\t"   << G << endl;
        cout << "Total Samples:\t" << S << endl;

    }

    void show(Array3D& array3d)
    {
        //cout << endl;
        //cout << "Total Genes:\t"   << G << endl;
        //cout << "Total Samples:\t" << S << endl;
      for (int i=0; i<(S-1); i++)
        {
	  //cout<<csample[i]<<"\t";
	  cout<<array3d.sName[csample[i]].str<<"\t";
        }
      cout<<array3d.sName[csample[(S-1)]].str<<"\n";

      for (int i=0; i<(G-1); i++)
        {
	  cout<<array3d.gName[cgene[i]].str<<"\t";
	  //cgene[i]<<"\t";
        }
      cout<<array3d.gName[cgene[(G-1)]].str<<"\n";


    }

    float kmeans(Array3D& array3d, int num_clusters, vector<GeneCluster *> gclusters, char* file_name, int write_to_file)
    {

        int num_repeats = 20;
        float mean[array3d.T];
        for (int j = 0 ; j< array3d.T ; j++)
        {
            mean[j] = 0;
        }

        //initialization -assign first k genes as centroids
        for (int i = 0; i < num_clusters; i++)
        {
            get_mean_profile(array3d,i,mean);
            gclusters[i]->setCentroid(mean);
            //gclusters[i]->displayCentroid();
        }

        for (int a = 0; a < num_repeats; a++)
        {

            //cout<<"Repeat number:" << a <<endl;
            // remove all elements in the cluster
            for (int m = 0; m < num_clusters; m++)
            {
                gclusters[m]-> removeElements();
            }
            //cout<<"Cluster is now empty - centorids are updated."<<endl;

            //put each document in the appropriate cluster
            for (int i = 0; i < G; i++)
            {
                //cout<<"Processing gene "<<i<<endl;
                int max_index = 0;
                float max_similarity = 0;

                for (int j = 0; j < num_clusters; j++)
                {

                    float* temp = gclusters[j]->getCentroid();
                    float t = pearsoncorrelation(temp, get_mean_profile(array3d, i, mean), array3d.T);

                    if ( t > max_similarity)
                    {
                        max_index = j;
                        max_similarity = t;
                    }
                    //cout << j << " " <<t << endl;
                }

                // Add gene to the correct cluster
                Profile temp_array;
                float* p =  get_mean_profile(array3d, i, mean);
                for (int b=0;b < LENGTH ; b++)
                    temp_array.p[b] = p[b];
                gclusters[max_index]->insertGene(temp_array,cgene[i]);

            }

#ifdef DEBUG
            for (int m = 0; m < num_clusters; m++)
            {
                cout<<"Cluster "<<m<<endl;
                gclusters[m]->displaySize();
            }
#endif

            // Recalculate centroids
            for (int m = 0; m < num_clusters; m++)
            {

#ifdef DEBUG
                cout<<"Before step" <<a<<endl;
                gclusters[m]->displayCentroid();
#endif

                gclusters[m]->recalculateCentroid();

#ifdef DEBUG
                cout<<"After step" <<a<<endl;
                gclusters[m]->displayCentroid();
#endif

            }
        }


        // calculate SSE
        float sse = 0;
        for (int n = 0; n < num_clusters; n++)
        {
            sse = sse + gclusters[n]->calculateSSE();
        }
       
        if (write_to_file != 0)
        {
            for (int m = 0; m < num_clusters; m++)
            {

#ifdef DEBUG
                gclusters[m]->displaySize();
#endif
                char result[80];
                sprintf( result, "%d", m );
                writeKmeansClustersFile(gclusters[m],array3d,result);
            }
        } //if



        // calculate Dunn's index and return to it
        // calculate maximum of intracluster distance
        //float max_within = 0;
        //float min_between = 10;
        /*for(int m = 0; m < num_clusters; m++){
        float temp = gclusters[m]->getWithinDissimilarity();
        if(temp > max_within){
        max_within = temp;
        }
        }

        float min_between = 10;
        for(int m = 0; m < (num_clusters -1); m++){
        for(int n = m; n < num_clusters; n++){
        float temp = gclusters[m]->getBetweenDissimilarity(gclusters[n]);
        temp = temp / max_within;
        if (temp < min_between){
        min_between = temp;
        }
        }
        }*/
        return sse;

    }

    void writeKmeansClustersFile(GeneCluster* gcluster, Array3D& array3d, char* file_name)
    {
        ofstream myfile;
        myfile.open (file_name);

        myfile << S<<"\n";
        for (int k = 0 ; k < S ; k++)
        {
            myfile<<csample[k]<<"\t";
        }
        myfile <<"\n";

        vector<int> mygenes = gcluster->getIds();
        int mysize = mygenes.size();

        myfile << mysize<<"\n";

        for (int j = 0; j < mysize ; j++)
        {
            //      float * mean = get_mean_profile(array3d, j);

            //for(int i=0; i<array3d.T; i++){
            //myfile << mean[i]<<" ";
            //}
            //myfile << endl;

            for (int k = 0 ; k < S ; k++)
            {
                float * submatrix = array3d.get3d(mygenes.at(j), csample[k]);

                for (int i=0; i<array3d.T; i++)
                {
                    myfile << submatrix[i]<<" ";
                }
                myfile << endl;
            }
        }
        myfile.close();
    }


    void writeFile(Array3D& array3d, char* file_name)
    {
        ofstream myfile;
        myfile.open (file_name);
        myfile << S<<"\n";
        for (int k = 0 ; k < S ; k++)
        {
            myfile<<csample[k]<<"\t";
        }
        myfile <<"\n";
        myfile << G<<"\n";
	for (int k = 0 ; k < G ; k++)
        {
            myfile<<cgene[k]<<"\t";
        }
	myfile <<"\n";

        for (int j = 0; j < G ; j++)
        {
            //      float * mean = get_mean_profile(array3d, j);

            //for(int i=0; i<array3d.T; i++){
            //myfile << mean[i]<<" ";
            //}
            //myfile << endl;

            for (int k = 0 ; k < S ; k++)
            {
                float * submatrix = array3d.get3d(cgene[j], csample[k]);

                for (int i=0; i<array3d.T; i++)
                {
                    myfile << submatrix[i]<<" ";
                }
                myfile << endl;
            }
        }
        myfile.close();
    }


        
    float calculate_similarity(Array3D& array3d)
    {

      float meanCor = 0; 
      int counter = 0;
      // For every pair of gene 
      // Get the mean profile
      // Calculate min(cross(m1,m2))
      for (int j = 0; j < (G-1) ; j++){
	  for(int k = j ; k < G ; k++){


	    float mean[array3d.T];
	    float mean2[array3d.T];
	    for (int i = 0 ; i< array3d.T ; i++)
	      {
		mean[i] = 0;
		mean2[i] = 0;
	      }
	    

            get_mean_profile(array3d,j,mean);
	    get_mean_profile(array3d,k,mean2);


	    /*	    
	    for (int i=0; i<array3d.T; i++)
	      {
		cout << mean[i]<<" ";
	      }
	    cout<<endl;

	    for (int i=0; i<array3d.T; i++)
	      {
		cout << mean2[i]<<" ";
	      }
	      cout<<endl;*/

	    meanCor  = meanCor + cross_cor(mean, mean2, LENGTH, LENGTH);
	    counter++;
	  }
      }
      
      return double(meanCor/counter);
    }
      
};


class Clustering
{
    char file[100], *p;
    int C;
    TabInput in;

public:

    Cluster* clusters;

    Clustering(char * name, int c)
    {
        clusters = new Cluster[c];

        for (int i=0; i < c ; i++)
        {
            clusters[i] = Cluster();
        }

        strcpy (file,name);
        C = c;
    }


    //QHu
    Clustering(vector< pair< vector<int>, vector<int> > > gs, int c)
    {
        clusters = new Cluster[c];

        for (int i=0; i < c ; i++)
        {
            clusters[i] = Cluster(gs[i].first, gs[i].second);
        }

        C = c;
    }

    void readTabFile()
    {

        in.open(file);

        for (int i=0; i<C; i++)
        {

            // get first line
            in.getLine();
            p = in.get2tab();
            while (p[0] !='\0' && p[0] != 0xd)
            {
                clusters[i].add_gene(atoi(p));
                //cout<<"Add "<<p<<" into Cluster"<<i<<endl;
                p = in.get2tab();
            }

            // get Second Line
            in.getLine();
            p = in.get2tab();
            while (p[0] !='\0' && p[0] != 0xd)
            {
                //cout<<"Add "<<p<<" into Cluster"<<i<<endl;
                clusters[i].add_sample(atoi(p));
                p = in.get2tab();
            }

        }
    }


    void show()
    {
        for (int i=0; i < C ; i++)
        {
            clusters[i].show_summary();
        }
    }

    void show_all(Array3D& array3d)
    {
        for (int i=0; i < C ; i++)
        {
            clusters[i].show(array3d);
        }
    }


};



class Input
{
    char file[100], *p;
    Array3D* pArray3d;
    TabInput in;
public:
    Input( char* name, Array3D& array )
    {
        strcpy( file, name );
        pArray3d = (Array3D*) &array;
    }

    void readTabFile()
    {
        int T, S, G;

        in.open( file );

        in.getLine();
        in.get2tab();
        p=in.get2tab(); // T
        pArray3d->T = T = atoi( p );

        in.getLine();
        in.get2tab();
        p=in.get2tab(); // S
        pArray3d->S = S = atoi( p );

        in.getLine();
        in.get2tab();
        p=in.get2tab(); // G
        pArray3d->G = G = atoi( p );

        pArray3d->init( T, S, G );

        for (int i=0; i<T; i++)
        {
            in.getLine();
            in.get2tab();
            p=in.get2tab();
            strcpy(pArray3d->tName[i].str, p);

            in.getLine();

            if ( i==0 )
            {
                in.get2tab();
                in.get2tab(); // ID and NAME
                for (int k=0; k<S; k++)
                {
                    p = in.get2tab();
                    strcpy(pArray3d->sName[k].str, p);
                }
            }

            for (int j=0; j<G; j++)
            {
                in.getLine();
                in.get2tab();
                p=in.get2tab(); // ID and NAME
                if ( i==0 )
                    strcpy(pArray3d->gName[j].str, p);

                for (int k=0; k<S; k++)
                {
                    p = in.get2tab();
                    // assert p is a float;
                    pArray3d->dat(j,k,i) = atof(p);
                    // TMP times 10
                    //if( k!= 0 );
                    //	cout << "\t";
                    //printf( "%2d", int(atof(p)*100) );
                }
                // TMP times 100
                //cout << endl;
            }
        }
        in.close();
    }
};


#endif
