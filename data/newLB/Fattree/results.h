#ifndef results_H
#define results_H
#include "Graphs.h"
#include "Trees.h"
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

#ifdef _Tree
void output2txt(int size,float t[],int numOfVM,float avgBw[],
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
	outFile<<"\n"<<"VM= ";
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

void static_embedding_N(char algorithm)
{
	Tree G;
	G.drawTree();
	//testKSP(G);
	//input
	const int n=6,m=1;
	int numOfreq=1e3,n_test=n;
	// param of physical networks
	float p=0.5,p_minResBw=0,p_maxResBw=1;// with p in [0,1G]
	// param of  VDC
	int numOfVM[n]={2,4,6,8,9,10};
	float minBw=100,maxBw=700;

	//output
	float max_utilization[n]={0},success_rate[n]={0},bandwidth_cost[n]={0};
	// start simulation
	
	time_t Atime,Btime; 
	
	float running_time[n]={0};
	
	for (int i=0;i<n_test;i++)
	{
		time(&Atime); 
		G.SingleRequest(algorithm,numOfreq,
			p,p_minResBw,p_maxResBw,numOfVM[i],minBw,maxBw,
			max_utilization[i],success_rate[i],bandwidth_cost[i]);	
		time(&Btime);
		running_time[i]=1000*(float)(Btime - Atime)/(float)numOfreq;
		cout<<"process time is "<<running_time[i]<<"ms"<<endl;
	}

	switch(algorithm)
	{
	case 'P':
			output2txt(n_test,running_time,numOfVM,minBw,maxBw,max_utilization,success_rate,bandwidth_cost,"s_pertubation.txt");
		break;
	case 'A':		
			output2txt(n_test,running_time,numOfVM,minBw,maxBw,max_utilization,success_rate,bandwidth_cost,"s_allocrange.txt");
		break;

	case 'R':
		
			output2txt(n_test,running_time,numOfVM,minBw,maxBw,max_utilization,success_rate,bandwidth_cost,"s_randomdrop.txt");
		break;
	case 'B':
		
			output2txt(n_test,running_time,numOfVM,minBw,maxBw,max_utilization,success_rate,bandwidth_cost,"s_backtracking.txt");
		break;
	case 'E':
			output2txt(n_test,running_time,numOfVM,minBw,maxBw,max_utilization,success_rate,bandwidth_cost,"s_exhaustive.txt");
		break;
	case 'F':
		 
			output2txt(n_test,running_time,numOfVM,minBw,maxBw,max_utilization,success_rate,bandwidth_cost,"s_firstfit.txt");
		break;
	case 'N':
		 
			output2txt(n_test,running_time,numOfVM,minBw,maxBw,max_utilization,success_rate,bandwidth_cost,"s_nextfit.txt");
		break;
	case 'G':
		 
			output2txt(n_test,running_time,numOfVM,minBw,maxBw,max_utilization,success_rate,bandwidth_cost,"s_bestfit.txt");
		break;
	default:
		break;
	}	
}


void static_embedding_B(char algorithm)
{
	Tree G;
	G.drawTree();
	//testKSP(G);
	//input
	const int n=5;
	int numOfreq=100,n_test=n; 
	// param of physical networks
	float p=1,p_minResBw=0.5,p_maxResBw=1; // with p in [500B,1G]
	// param of  VDC
	int numOfVM=12;
	float devBw=100,avgBw[n]={200,250,300,350,400};
	//output
	float max_utilization[n]={0},success_rate[n]={0},bandwidth_cost[n]={0};
	// start simulation
	
	time_t Atime,Btime; 
	
	float running_time[n]={0};

		
		for (int i=0;i<n_test;i++)
		{
			time(&Atime); 
			G.SingleRequest(algorithm,numOfreq,
				p,p_minResBw,p_maxResBw,numOfVM,avgBw[i]-devBw,avgBw[i]+devBw,
				max_utilization[i],success_rate[i],bandwidth_cost[i]);	

			time(&Btime);
			running_time[i]=1000*(float)(Btime - Atime)/(float)(numOfreq*n_test);
			cout<<"process time is "<<running_time[i]<<"ms"<<endl;
		}

		switch(algorithm)
		{
		case 'P':
				output2txt(n_test,running_time,numOfVM,avgBw,max_utilization,success_rate,bandwidth_cost,"s_pertubation.txt");
			break;
		case 'A':			 
				output2txt(n_test,running_time,numOfVM,avgBw,max_utilization,success_rate,bandwidth_cost,"s_allocrange.txt");
			break;
		case 'R':			 
				output2txt(n_test,running_time,numOfVM,avgBw,max_utilization,success_rate,bandwidth_cost,"s_randomdrop.txt");
			break;

		case 'B':
			
				output2txt(n_test,running_time,numOfVM,avgBw,max_utilization,success_rate,bandwidth_cost,"s_backtracking.txt");
			break;
		case 'E':
				output2txt(n_test,running_time,numOfVM,avgBw,max_utilization,success_rate,bandwidth_cost,"s_exhaustive.txt");
			break;
		case 'F':
			 
				output2txt(n_test,running_time,numOfVM,avgBw,max_utilization,success_rate,bandwidth_cost,"s_firstfit.txt");
			break;
		default:
			break;
		}	
	
}



void dynamic_embedding(char algorithm)
{
	Tree G;
	G.drawTree();
	//input
	int numOfreq=1e3;
	const int n=4,m=2;
	float load[n]={0.2,0.4,0.6,0.8};
	int n_test=4;
	//output
	float max_utilization[n]={0},success_rate[n]={0},bw_cost[n]={0},RC[n]={0};
	// start simulation	
	time_t Atime,Btime; 
	time(&Atime); 
	float running_time=0;
	
	for (int i=0;i<n_test;i++)
		G.ProcessRequest(algorithm,numOfreq,load[i],max_utilization[i],success_rate[i],bw_cost[i],RC[i]);
	
	time(&Btime);
	running_time=1000*1000*(float)(Btime - Atime)/(float)(numOfreq*n_test);
	cout<<"process time is "<<running_time<<"us"<<endl;
	
	switch(algorithm)
	{
		case 'P':
			output2txt(n_test,running_time,load,max_utilization,success_rate,bw_cost,RC,"tree_pertubation.txt");
			break;
		case 'B':
			output2txt(n_test,running_time,load,max_utilization,success_rate,bw_cost,RC,"tree_backtracking.txt");
			break;
		case 'F':			 
			output2txt(n_test,running_time,load,max_utilization,success_rate,bw_cost,RC,"tree_firstfit.txt");
			break;
		default:
			break;
	}	
}
#endif
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