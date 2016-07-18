/* single path - bottleneck & bandwidth cost
   输出数据到xx2.txt 
*/

#include<iostream>
#include<fstream>
#include<cmath>
#include<ctime>
using namespace std;
ofstream fout;

#define cycle 15
#define num_i 4
#define num_s 10

//---------拓扑基本参数--------
int space[num_s]; // node space
int link[num_i][num_s]; // link capacity

//-----中间变量-----
int B,N,times;
int A[num_s+1]; // 选进去的点集
int sum_m; //累计VM数量
double total_band=0;

double max_load;
double bw_cost;
double max_load_array[cycle];
double bw_cost_array[cycle];
int success[cycle]={0};

//----分配结果----
int m[num_s];
int f[num_s][num_s];

void initialize()
{
	int i;
	for(i=0;i<num_s;i++)
	{
		max_load_array[i] = -1;
		bw_cost_array[i] = -1;
	}
}

void generate()
{
	int i,j;
	total_band = 0;

	//srand(10*times);
	//----node space分布---
	for(i=0;i<num_s;i++) //只分布leaf层
		space[i]=rand()%5+1;

	//---link capacity 分布---
	for(i=0;i<num_i;i++)
		for(j=0;j<num_s;j++)
		{
			link[i][j]=rand()%100+1;
			total_band += link[i][j];
		}

	//-----输出拓扑信息-----
	//cout<<"=========================\n";
	//cout<<"round "<<times<<endl;
	/*cout<<"links:\n";
	for(i=0;i<num_i;i++)
	{
		for(j=0;j<num_s;j++)
			cout<<link[i][j]<<" ";
		cout<<"\t";
	}
	cout<<endl;

	cout<<"space:\n";
	for(i=0;i<num_s;i++)
		cout<<space[i]<<" ";
	cout<<endl;*/
}

void clear()
{
	int i,j;

	sum_m = 0;
	max_load = 0.0;
	bw_cost = 0.0;

	for(i=0;i<num_s;i++)
		m[i]=0;

	for(i=0;i<=num_s;i++)
		A[i]=0;

	for(i=0;i<num_s;i++)
		for(j=0;j<num_s;j++)
				f[i][j]=-1;
}
int check_single(int s)
{
	int j,k,i,max_path;
	m[s]=space[s];

	//---带宽分约束---
	for(j=1;j<=A[0];j++)
	{
		k=A[j];
		max_path=0;

		// single path
		for(i=0;i<num_i;i++)
			if( min(link[i][s], link[i][k]) > max_path)
			{
				max_path = min(link[i][s], link[i][k]);
				f[s][k] = i;
				f[k][s] = i;
			}

		while( max_path < B* min(m[s],m[k]) )
			m[s]--;
	}

	sum_m += m[s];

	return 1;
}

int check_all(int i, int j)
{ 
	int s,d,k,sum_n,n,mm;
	double sum_t;
	
	s = A[j];
	sum_n=0; // # of VM
	sum_t=0; // # of traffic

	//--作为判断条件的traffic，（类似取f最大 ，f只能取离散值0,1)
	k=1; n = min(m[s], sum_m-m[s]);
	while(k<=A[0] && sum_n< n)
	{
		if(A[k]==s) k++;
		d = A[k];

		if(f[s][d]==i) //表示path在该switch上
		{
			mm = min( n-sum_n , m[d]); //最坏就 mm = n-0

			sum_n += mm;
			sum_t += B*mm;
		}
		k++;
	}

	if(sum_t > link[i][s])
		return 0;
	
	//
	bw_cost += sum_t;
	max_load = max( max_load, sum_t/link[i][s] );

	return 1;
}

int greedy()
{
	int s,i,j;

	//---greedy算法---
	for(s=0;s<num_s;s++) 
	{
		check_single(s);

		if(m[s]>0)
		{
			A[++A[0]]=s;

			//--check_all
			for(j=1;j<=A[0];j++)
				for(i=0;i<num_i;i++) 
					while(!check_all(i,j))	{
						m[A[j]]--; sum_m--;
					}

			//--找到，收尾
			if(sum_m >= N)
			{
				m[s] -= sum_m - N;

				//re-calculate load (因为m[ss]变动了), 计算bw_cost
				max_load = 0;
				bw_cost = 0;
				for(j=1;j<=A[0];j++)
					for(i=0;i<num_i;i++) 
						check_all(i,j);
				return 1;
			}
		}

		
	} 
	return 0;
}

void record(int flag)
{
	if(flag==1) //success
	{
		success[N-1]++;
		if(bw_cost_array[N-1] == -1)
			bw_cost_array[N-1] = bw_cost/total_band;
		else
			bw_cost_array[N-1] = ( bw_cost_array[N-1]*(success[N-1]-1) + bw_cost/total_band )/success[N-1];

		if(max_load_array[N-1] == -1)
			max_load_array[N-1] = max_load;
		else
			max_load_array[N-1] = ( max_load_array[N-1]*(success[N-1]-1) + max_load )/success[N-1];
	}
}

void figure()
{
	int i;

	fout.open("E:\\data\\N-bn\\bottleneck2.txt");
	fout.clear();
	fout.close();
	fout.open("E:\\data\\N-bw\\bandcost2.txt");
	fout.clear();
	fout.close();
	fout.open("E:\\data\\N-suc\\success2.txt");
	fout.clear();
	fout.close();

	fout.open("E:\\data\\N-bn\\bottleneck2.txt",ios::app);
	for(i=0;i<cycle;i++) //(N=1;N<=20;N++)
		fout<<max_load_array[i]<<" ";
	fout.close();

	fout.open("E:\\data\\N-bw\\bandcost2.txt",ios::app);
	for(i=0;i<cycle;i++) //(N=1;N<=20;N++)
		fout<<bw_cost_array[i]<<" ";
	fout.close();

	fout.open("E:\\data\\N-suc\\success2.txt",ios::app);
	for(i=0;i<cycle;i++) //(N=1;N<=20;N++)
		fout<<double(success[i])/(times-1)<<" ";
	fout.close();
}

void main()
{
	time_t begin,end;
	initialize();

	begin=clock();
	for(times=1;times<=10;times++)
	{
		generate();// hardware network
	
		B=25;//N=10
		for(N=1;N<=cycle;N++) //(B=1;B<=cycle;B++)
		{
			clear();
			record(greedy());
		}
	}
	figure();

	end=clock();
	cout<<"run time: "<<double(end-begin)/CLOCKS_PER_SEC<<" s"<<endl;
	system("pause");
}