#ifndef results_H
#define results_H
#include "Graphs.h"
#include "Trees.h"
#include "parameter.h"
#include <ctime>
using namespace std;
void testKSP(Graph& G)
{
	int s,d;
	{
		for (int t=1;t<=8;t++)
		{
			cout<<"\n"<<"Input: source dest"<<"\n";
			if(cin>>s&&cin>>d&&s!=d)
			{	s--;d--;
				G.YenKSP(s,d);
				for (int k=0;k<G.n_path[s][d];k++){
					for (unsigned short i=1;i<=G.Kpaths[s][d][k][0]+1;i++)
						cout<<G.Kpaths[s][d][k][i]+1<<"-";
					cout<<"\n";
				}
			}
		}
	}

}


void output2txt(float* x,int size,char *filename)
{
	// open input file
	ofstream outFile(filename);

	if (!outFile) {
	cerr << "Cannot open" <<filename<<" for output\n";
	}
        // copy file contents to cout

	int i=0;
    while (i<size) {
         	outFile<<x[i]<<" ";
		 i++;
    }
	outFile<<endl;
	outFile.close();
}
// results from dynamic simulation
void output2txt(int size,float t,float load[],
				float max_utilization[],float success_rate[],float bw[],float RC[],char *filename)
{
	// open input file
	ofstream outFile(filename);

	if (!outFile) {
	cerr << "Cannot open" <<filename<<" for output\n";
	}
#ifdef _FatTree
	outFile<<"\n"<<"FatTree"<<"\n";
#endif
#ifdef _L2VL2
	outFile<<"\n"<<"L2VL2"<<"\n";
#endif
#ifdef _VL2
	outFile<<"\n"<<"VL2"<<"\n";
#endif
#ifdef _Bcube
	outFile<<"\n"<<"Bcube"<<"\n";
#endif

	outFile<<"\n"<<"N=( "<<minN<<" , "<<maxN<<" )";
	outFile<<"\n"<<"B=( "<<minB<<" , "<<maxB<<" )";
	outFile<<"\n"<<"max_Slot= "<<maxSlot;
#ifndef _Tree
	outFile<<"\n"<<"process time is "<<t<<"ms"; 
#else
	outFile<<"\n"<<"process time is "<<t<<"us"; 
#endif
	int i=0;
	outFile<<"\n"<<"load = ";
    while (i<size) {
         	outFile<<load[i]<<" ";
		 i++;
    }
	i=0;
	outFile<<"\n"<<"MLU = ";
    while (i<size) {
         	outFile<<max_utilization[i]<<" ";
		 i++;
    }
	i=0;
	outFile<<"\n"<<"success_rate = ";
    while (i<size) {
         	outFile<<success_rate[i]<<" ";
		 i++;
    }
	i=0;
	outFile<<"\n"<<"R/C = ";
    while (i<size) {
         	outFile<<RC[i]<<" ";
		 i++;
    }
	i=0;
	outFile<<"\n"<<"bandwidth_cost = ";
    while (i<size) {
         outFile<<bw[i]<<" ";
		 i++;
    }
	outFile<<endl;
	outFile.close();
}
//K paths results from results from static simulation:N=10,varying B
void output2txt(Graph& G,int size,float t[],int numOfVM,float avgBw[],
				float max_utilization[],float success_rate[],float bw[],char *filename)
{
	// open input file
	//ofstream outFile(filename);
	ofstream outFile(filename, ios_base::out|ios_base::app);
	if (!outFile) {
	cerr << "Cannot open" <<filename<<" for output\n";
	}
#ifdef _FatTree
	outFile<<"\n"<<"FatTree"<<"\n";
#endif
#ifdef _L2VL2
	outFile<<"\n"<<"L2VL2"<<"\n";
#endif
#ifdef _VL2
	outFile<<"\n"<<"VL2"<<"\n";
#endif
#ifdef _Bcube
	outFile<<"\n"<<"Bcube"<<"\n";
#endif

	outFile<<"\n"<<"Kmax= "<<G.Kmax;
	outFile<<"\n"<<"N= "<<numOfVM;
	//outFile<<"\n"<<"B=( "<<minB<<" , "<<maxB<<" )";
	outFile<<"\n"<<"max_Slot= "<<maxSlot;
	
	int i=0;
	
	outFile<<"\n"<<"process time= ";
	
    while (i<size) {
         outFile<<t[i]<<"ms ";
		 i++;
    }
	i=0;
	outFile<<"\n"<<"avgBw= ";
    while (i<size) {
         outFile<<avgBw[i]<<" ";
		 i++;
    }
	i=0;
	outFile<<"\n"<<"MLU = ";
    while (i<size) {
         outFile<<max_utilization[i]<<" ";
		 i++;
    }
	i=0;
	outFile<<"\n"<<"success_rate = ";
    while (i<size) {
         outFile<<success_rate[i]<<" ";
		 i++;
    }
	i=0;
	outFile<<"\n"<<"bandwidth_cost = ";
    while (i<size) {
         outFile<<bw[i]<<" ";
		 i++;
    }
	
	outFile<<endl;
	outFile.close();
}

// results from static simulation,N={2,4,6,8,10},B=100-700
void output2txt(Graph& G,int size,float t[],int numOfVM[],float minBw,float maxBw,
				float max_utilization[],float success_rate[],float bw[],char *filename)
{
	// open input file
	ofstream outFile(filename, ios_base::out|ios_base::app);//ofstream outFile(filename);

	if (!outFile) {
	cerr << "Cannot open" <<filename<<" for output\n";
	}
#ifdef _FatTree
	outFile<<"\n"<<"FatTree"<<"\n";
#endif
#ifdef _L2VL2
	outFile<<"\n"<<"L2VL2"<<"\n";
#endif
#ifdef _VL2
	outFile<<"\n"<<"VL2"<<"\n";
#endif
#ifdef _Bcube
	outFile<<"\n"<<"Bcube"<<"\n";
#endif
	outFile<<"\n"<<"Kmax= "<<G.Kmax;
	//outFile<<"\n"<<"N=( "<<minN<<" , "<<maxN<<" )";
	outFile<<"\n"<<"B=( "<<minBw<<" , "<<maxBw<<" )";
	outFile<<"\n"<<"max_Slot= "<<maxSlot;
	
	int i=0;
	
	outFile<<"\n"<<"process time= ";
	
    while (i<size) {
         outFile<<t[i]<<"ms ";
		 i++;
    }
	i=0;
	outFile<<"\n"<<"VM= ";
    while (i<size) {
         outFile<<numOfVM[i]<<" ";
		 i++;
    }
	i=0;
	outFile<<"\n"<<"MLU = ";
    while (i<size) {
         outFile<<max_utilization[i]<<" ";
		 i++;
    }
	i=0;
	outFile<<"\n"<<"success_rate = ";
    while (i<size) {
         outFile<<success_rate[i]<<" ";
		 i++;
    }
	i=0;
	outFile<<"\n"<<"bandwidth_cost = ";
    while (i<size) {
         outFile<<bw[i]<<" ";
		 i++;
    }
	
	outFile<<endl;
	outFile.close();
}
void output2txt(Graph& G,int size,float t[],int numberOfGroup,int numOfVM[],float minBw,float maxBw,
				float max_utilization[],float success_rate[],float bw[],char *filename)
{
	// open input file
	ofstream outFile(filename, ios_base::out|ios_base::app);//ofstream outFile(filename);

	if (!outFile) {
	cerr << "Cannot open" <<filename<<" for output\n";
	}
#ifdef _FatTree
	outFile<<"\n"<<"FatTree"<<"\n";
#endif
#ifdef _L2VL2
	outFile<<"\n"<<"L2VL2"<<"\n";
#endif
#ifdef _VL2
	outFile<<"\n"<<"VL2"<<"\n";
#endif
#ifdef _Bcube
	outFile<<"\n"<<"Bcube"<<"\n";
#endif
	outFile<<"\n"<<"Kmax= "<<G.Kmax;
	//outFile<<"\n"<<"N=( "<<minN<<" , "<<maxN<<" )";
	outFile<<"\n"<<"B=( "<<minBw<<" , "<<maxBw<<" )";
	outFile<<"\n"<<"max_Slot= "<<maxSlot;

	outFile << "\n"<<"number of group ="<<numberOfGroup;
	
	int i=0;
	
	outFile<<"\n"<<"process time= ";
	
    while (i<size) {
         outFile<<t[i]<<"ms ";
		 i++;
    }
	i=0;
	outFile<<"\n"<<"VM= ";
    while (i<size) {
         outFile<<numOfVM[i]<<" ";
		 i++;
    }
	i=0;
	outFile<<"\n"<<"MLU = ";
    while (i<size) {
         outFile<<max_utilization[i]<<" ";
		 i++;
    }
	i=0;
	outFile<<"\n"<<"success_rate = ";
    while (i<size) {
         outFile<<success_rate[i]<<" ";
		 i++;
    }
	i=0;
	outFile<<"\n"<<"bandwidth_cost = ";
    while (i<size) {
         outFile<<bw[i]<<" ";
		 i++;
    }
	
	outFile<<endl;
	outFile.close();
}


void output2txt(Tree& G,int size,float t[],int numOfVM[],float minBw,float maxBw,
				float max_utilization[],float success_rate[],float bw[],char *filename)
{
	// open input file
	ofstream outFile(filename, ios_base::out|ios_base::app);//ofstream outFile(filename);

	if (!outFile) {
	cerr << "Cannot open" <<filename<<" for output\n";
	}
#ifdef _Tree
	outFile<<"\n"<<"Tree"<<"\n";
#endif

	//outFile<<"\n"<<"N=( "<<minN<<" , "<<maxN<<" )";
	outFile<<"\n"<<"B=( "<<minBw<<" , "<<maxBw<<" )";
	outFile<<"\n"<<"max_Slot= "<<maxSlot;
	
	int i=0;
	
	outFile<<"\n"<<"process time= ";
	
    while (i<size) {
         outFile<<t[i]<<"ms ";
		 i++;
    }
	i=0;
	outFile<<"\n"<<"VM= ";
    while (i<size) {
         outFile<<numOfVM[i]<<" ";
		 i++;
    }
	i=0;
	outFile<<"\n"<<"MLU = ";
    while (i<size) {
         outFile<<max_utilization[i]<<" ";
		 i++;
    }
	i=0;
	outFile<<"\n"<<"success_rate = ";
    while (i<size) {
         outFile<<success_rate[i]<<" ";
		 i++;
    }
	i=0;
	outFile<<"\n"<<"bandwidth_cost = ";
    while (i<size) {
         outFile<<bw[i]<<" ";
		 i++;
    }
	
	outFile<<endl;
	outFile.close();
}


void output2txt(Tree& G,int size,float t[],int numOfVM,float avgBw[],
				float max_utilization[],float success_rate[],float bw[],char *filename)
{
	// open input file
	//ofstream outFile(filename);
	ofstream outFile(filename, ios_base::out|ios_base::app);
	if (!outFile) {
	cerr << "Cannot open" <<filename<<" for output\n";
	}
#ifdef _Tree
	outFile<<"\n"<<"Tree"<<"\n";
#endif

	outFile<<"\n"<<"N= "<<numOfVM;
	//outFile<<"\n"<<"B=( "<<minB<<" , "<<maxB<<" )";
	outFile<<"\n"<<"max_Slot= "<<maxSlot;
	
	int i=0;
	
	outFile<<"\n"<<"process time= ";
	
    while (i<size) {
         outFile<<t[i]<<"ms ";
		 i++;
    }
	i=0;
	outFile<<"\n"<<"avgBw= ";
    while (i<size) {
         outFile<<avgBw[i]<<" ";
		 i++;
    }
	i=0;
	outFile<<"\n"<<"MLU = ";
    while (i<size) {
         outFile<<max_utilization[i]<<" ";
		 i++;
    }
	i=0;
	outFile<<"\n"<<"success_rate = ";
    while (i<size) {
         outFile<<success_rate[i]<<" ";
		 i++;
    }
	i=0;
	outFile<<"\n"<<"bandwidth_cost = ";
    while (i<size) {
         outFile<<bw[i]<<" ";
		 i++;
    }
	
	outFile<<endl;
	outFile.close();
}

void output2txt(Graph& G,int size,float t_RA[],int numOfVM[],float minBw,float maxBw,char *filename)
{
	// open input file
	//ofstream outFile(filename);
	ofstream outFile(filename, ios_base::out|ios_base::app);
	if (!outFile) {
	cerr << "Cannot open" <<filename<<" for output\n";
	}
#ifdef _FatTree
	outFile<<"\n"<<"FatTree"<<"\n";
#endif
#ifdef _L2VL2
	outFile<<"\n"<<"L2VL2"<<"\n";
#endif
#ifdef _VL2
	outFile<<"\n"<<"VL2"<<"\n";
#endif
#ifdef _Bcube
	outFile<<"\n"<<"Bcube"<<"\n";
#endif

	outFile<<"\n"<<"Kmax= "<<G.Kmax;
	outFile<<"\n"<<"B=( "<<minBw<<" , "<<maxBw<<" )";
	outFile<<"\n"<<"max_Slot= "<<maxSlot;
	
	int i=0;

	outFile<<"\n"<<"VM= ";
    while (i<size) {
         outFile<<numOfVM[i]<<" ";
		 i++;
    }
	i=0;
	outFile<<"\n"<<"process time= ";
	
    while (i<size) {
         outFile<<t_RA[i]<<"ms ";
		 i++;
    }
	

    
	
	outFile<<endl;
	outFile.close();
}


//oversub test
void output2txt(int size,float t[],float sub_tor[],float sub_agg[],int numOfVM1,int numOfVM2,float minBw,float maxBw,
				float max_utilization[],float success_rate[],float bw[],char *filename)
{
	// open input file
	ofstream outFile(filename, ios_base::out|ios_base::app);//ofstream outFile(filename);

	if (!outFile) {
	cerr << "Cannot open" <<filename<<" for output\n";
	}
#ifdef _FatTree
	outFile<<"\n"<<"FatTree"<<"\n";
#endif
#ifdef _L2VL2
	outFile<<"\n"<<"L2VL2"<<"\n";
#endif
#ifdef _VL2
	outFile<<"\n"<<"VL2"<<"\n";
#endif
#ifdef _Bcube
	outFile<<"\n"<<"Bcube"<<"\n";
#endif
	//outFile<<"\n"<<"Kmax= "<<G.Kmax;
	outFile<<"\n"<<"N=( "<<numOfVM1<<" , "<<numOfVM2<<" )";
	outFile<<"\n"<<"B=( "<<minBw<<" , "<<maxBw<<" )";
	outFile<<"\n"<<"max_Slot= "<<maxSlot;
	
	int i=0;

	i=0;
	outFile<<"\n"<<"sub_tor = ";
    while (i<size) {
         outFile<<sub_tor[i]<<" ";
		 i++;
    }
		i=0;
	outFile<<"\n"<<"sub_agg = ";
    while (i<size) {
         outFile<<sub_agg[i]<<" ";
		 i++;
    }

	i=0;
	outFile<<"\n"<<"process time= ";
	
    while (i<size) {
         outFile<<t[i]<<"ms ";
		 i++;
    }
	

	i=0;
	outFile<<"\n"<<"MLU = ";
    while (i<size) {
         outFile<<max_utilization[i]<<" ";
		 i++;
    }
	i=0;
	outFile<<"\n"<<"success_rate = ";
    while (i<size) {
         outFile<<success_rate[i]<<" ";
		 i++;
    }
	i=0;
	outFile<<"\n"<<"bandwidth_cost = ";
    while (i<size) {
         outFile<<bw[i]<<" ";
		 i++;
    }
	
	outFile<<endl;
	outFile.close();
}


void output2txt(int k_path,int size,float t[],float sub_tor[],float sub_agg[],int numOfVM,float minBw,float maxBw,
				float max_utilization[],float success_rate[],float bw[],char *filename)
{
	// open input file
	ofstream outFile(filename, ios_base::out|ios_base::app);//ofstream outFile(filename);

	if (!outFile) {
	cerr << "Cannot open" <<filename<<" for output\n";
	}
#ifdef _FatTree
	outFile<<"\n"<<"FatTree"<<"\n";
#endif
#ifdef _L2VL2
	outFile<<"\n"<<"L2VL2"<<"\n";
#endif
#ifdef _VL2
	outFile<<"\n"<<"VL2"<<"\n";
#endif
#ifdef _Bcube
	outFile<<"\n"<<"Bcube"<<"\n";
#endif
	outFile<<"\n"<<"Kmax= "<<k_path;
	outFile<<"\n"<<"N= "<<numOfVM;
	outFile<<"\n"<<"B=( "<<minBw<<" , "<<maxBw<<" )";
	outFile<<"\n"<<"max_Slot= "<<maxSlot;
	
	int i=0;

	i=0;
	outFile<<"\n"<<"sub_tor = ";
    while (i<size) {
         outFile<<sub_tor[i]<<" ";
		 i++;
    }
		i=0;
	outFile<<"\n"<<"sub_agg = ";
    while (i<size) {
         outFile<<sub_agg[i]<<" ";
		 i++;
    }

	i=0;
	outFile<<"\n"<<"process time= ";
	
    while (i<size) {
         outFile<<t[i]<<"ms ";
		 i++;
    }
	

	i=0;
	outFile<<"\n"<<"MLU = ";
    while (i<size) {
         outFile<<max_utilization[i]<<" ";
		 i++;
    }
	i=0;
	outFile<<"\n"<<"success_rate = ";
    while (i<size) {
         outFile<<success_rate[i]<<" ";
		 i++;
    }
	i=0;
	outFile<<"\n"<<"bandwidth_cost = ";
    while (i<size) {
         outFile<<bw[i]<<" ";
		 i++;
    }
	
	outFile<<endl;
	outFile.close();
}


#endif