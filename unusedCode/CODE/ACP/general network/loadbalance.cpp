/* Perturbation algorithm with load balance on general network - ex.1 */
/*
	1. 本程序中，一般数组，0位存放的是数组大小信息。
	2. 因为只是分配一次，所以直接用link,若是动态仿真，需用res_link 
*/

#include<iostream>
#include<fstream>
#include<cmath>
#include<ctime>
using namespace std;
ofstream fout;

//----仿真规模---
#define K 4

#define B 1
#define N 2

#define n_switch 3
#define n_server 2

#define n_node n_switch+n_server
#define my_MAX 1000

//----硬件资源-----
int switch_set[n_switch];
int server_set[n_server];

int space[n_node]; // server
int weigh[n_node][n_node]; // link length
double link[n_node][n_node]; // link capacity

//-----中间变量-----
int A[n_server+1]; // 选进去的点集
int sum_m; //累计VM数量
int bottleneck_s; //定位的导致bottleneck的server

//Dijstra()
int dis[n_switch+1]; //Dijkstra中点到s的距离
int weigh_t[n_node][n_node]; // Dijksrta中动态weigh
int D[n_switch+1]; //switch set + d
int pre[n_switch+1];

int n_D;
int n_path;
int *pathSet[K][n_node+1]; //当前s,d之间path信息
double capacity[n_switch]; // 当前各path的capacity
int passby_s[n_node][n_node][n_switch+2];
int passby_d[n_node][n_node][n_switch+2];

//----分配结果----
int m[n_node];
double f[n_node][n_node][n_node][n_node];

//----统计参数----
double max_load;
double bw_cost;
int suc;

void Generate()
{
	/* generate harware resource */

	/* server_set,space中从0编到n;
	   weigh(距离),link(带宽)中编号为混编 */

	int i,j=0;
	// node 分布
	for(i=0;i<n_switch;i++)
		switch_set[i] = j++; //0,1
	for(i=0;i<n_server;i++)
		server_set[i] = j++; //2,3

	// VM space 
	for(i=0;i<n_node;i++) //只分布leaf层
		space[i] = 1;
	
	// link length
	for(i=0;i<n_node;i++)
	{
		for(j=0;j<n_node;j++)
			weigh[i][j] = my_MAX;

		weigh[i][i] = 0;
	}

	// link capacity
	for(i=0;i<n_node;i++)
		for(j=0;j<n_node;j++)
			link[i][j] = 0;

	weigh[0][1] = 1; link[0][1] = 1; 
	weigh[0][3] = 1; link[0][3] = 1;
	weigh[1][0] = 1; link[1][0] = 1;
	weigh[1][2] = 1; link[1][2] = 5; 
	weigh[1][3] = 1; link[1][3] = 2;
	weigh[1][4] = 1; link[1][4] = 3; 
	weigh[2][1] = 1; link[2][1] = 5; 
	weigh[2][4] = 1; link[2][4] = 1;
	weigh[3][0] = 1; link[3][0] = 1;
	weigh[3][1] = 1; link[3][1] = 2; 
	weigh[4][1] = 1; link[4][1] = 3;
	weigh[4][2] = 1; link[4][2] = 1; 
}

void Initialize()
{
	/* Initialize 分配结果 and 统计变量 */
	int s,d,i,j;

	sum_m = 0; 
	max_load = 0.0;
	bw_cost = 0.0;
	suc = 0;

	for(i=0;i<=n_server;i++)
		A[i] = 0;

	for(i=0;i<n_node;i++)
		m[i]=0;

	for(s=0;s<n_node;s++)
		for(d=0;d<n_node;d++)
			for(i=0;i<n_node;i++)
				for(j=0;j<n_node;j++)
					f[s][d][i][j] = 0; //f-sd-e(i,j),注：此处的s,d非点下标，而是server下标

	for(i=0;i<n_node;i++)
		for(j=0;j<n_node;j++)
			weigh_t[i][j] = weigh[i][j];

	for(i=0;i<n_node;i++)
		for(j=0;j<n_node;j++)
		{
			passby_s[i][j][0] = 0;
			passby_d[i][j][0] = 0;
		}
}

void Initial(int s, int d)
{
	/* 本函数为Dijksrta算法中weigh情况Initialize */

	int i;

	//Initialize D[]
	n_D = n_switch+1;
	D[0] = d; 
	for(i=1;i<=n_switch;i++)
		D[i] = switch_set[i-1];

	//Initilize dis[]
	for(i=0;i<n_D;i++)
		dis[D[i]] = weigh_t[s][D[i]];

	//Initialize pre[]
	for(i=0;i<n_D;i++)
		pre[i] = s;

	
	
}
int* Dijkstra(int s,int d)
{
	/* s为source节点，d为destination节点 */
	int i,j,k,j_index;

	n_path = 0;
	Initial(s,d);

	//---算法---
	while(n_D > 0)
	{
		cout<<"D: ";
		for(i=0;i<n_D;i++)
			cout<<D[i]<<" ";
		cout<<endl;

		for(i=0;i<n_D;i++)
			cout<<"dis["<<D[i]<<"] = "<<dis[D[i]]<<", ";
		cout<<endl;
	
		//j=d,找到一条路
		if(dis[d]<my_MAX)
		{
			path[n_path][1]=d; //path[][0]存放数量信息
			capacity[n_path]=my_MAX;

			j=d;
			for(i=2;i<=n_switch+2;i++)
			{
				path[n_path][i] = pre[j]; //记录path路径
	
				if( link[pre[j]][j] < capacity[n_path] )
 					capacity[n_path] = link[pre[j]][j]; //计算path.capacity.

				weigh_t[pre[j]][j] = my_MAX; //remove links on current path

				j = pre[j];
				if(j==s) 
				{
					path[n_path][0]=i;
					break;
				}
			}

			cout<<"---find a path: ";
			for(i=1;i<=path[n_path][0];i++)
				cout<<path[n_path][i]<<" ";
			cout<<endl;

			n_path++;
			Initial(s,d);
		}

		//j为中途点
		else 
		{
			//choose min dis[s][j]
			j=d;
			for(i=1;i<n_D;i++) 
				if(dis[D[i]] < dis[j])
					j=D[i];

			//搜索结束，返回
			if(dis[j] >= my_MAX)
				return;

			//借助点j，跟新路由
			for(i=0;i<n_D;i++)
			{
				k = D[i];
				if( dis[j] + weigh_t[j][k] < dis[k])
				{
					dis[k] = dis[j] + weigh_t[j][k];
					pre[k] = j;
				}
			}

			//remove j from D
			for(i=0;i<n_D;i++)
				if(D[i]==j)
					j_index = i;

			for(i=j_index;i<n_D-1;i++)
				D[i]=D[i+1];

			n_D--;
		}
	}
}

void KSP(int s, int d)
{
	k=0;
	pathSet[0] = Dijkstra(s,d);

	while( k < K-1 )
	{
		for(i=0;i<pathSet[k][0];i++) //对path中的每条边依次删去
		{
			spurNode = pathSet[0][i];

			//删去使用过的边
			for(j=0;j<=k;j++) 
			{
				flag = 0;
				for( m=1;m<=i;m++)
					if( pathSet[j][m] != pathSet[k][m] )
						flag = 1;

				if(flag == 0)
					weigh_t[pathSet[j][i]][pathSet[j][i+1]] = my_MAX;
			}

			spurPath = Dijkstra( pathSet[k][i],d ); 

			for(j=0;j<i;j++)
				pathBuffer[size_b][j] = pathSet[k][j];
			for(j=i;j<spur_size;j++)
				pathBuffer[size_b][j] = spurPath[j];

			ResetWeigh_t();
		}

		//选出B中最短的路径
		for(i=0;i<size_B;i++)
			if(pathBuffer[i][0] < min_length)
				min_index = i;

		//加入A
		k++; 
		for(i=1;i<=min_length;i++)
			A[k][i] = pathBuffer[min_index][i];
		
		//从B中删去
		for(i=min_index;i<size_B;i++)
			for(j=0;j<10;j++)
				pathBuffer[i][j] = pathBuffer[i+1][j];
		size_B--;
	}
}

int Routing(int s,int d)
{
	int j,i,p,q,n;
	
	//--- narrow the candidate set ---
	double sum_capacity = 0;

	for(i=0;i<n_path;i++) 
		sum_capacity += capacity[i];

	while(sum_capacity < B* min(m[s],m[d]))
		m[s]--;

	//---route with load balance---

	for(i=0;i<n_path;i++) 
	{
		for(j=1;j<path[i][0];j++)
		{
			q=path[i][j];p=path[i][j+1]; //link e(p,q)

			//f-sd在e(p,q)上的分配比例
			f[s][d][p][q] = capacity[i]/sum_capacity ;
			f[d][s][p][q] = f[s][d][p][q]; // p->q是考虑方向的

			//e(p,q)上存储(s,d)信息
			passby_s[p][q][0]++;   
			passby_d[p][q][0]++;
			n = passby_s[p][q][0];

			passby_s[p][q][n] = s; 
			passby_d[p][q][n] = d;

			//cout<<"passby["<<p<<"]["<<q<<"]: ("<<s<<","<<d<<")\n";
		}
		cout<<"f on path "<<i<<": "<<f[d][s][p][q]<<endl;
	}

	return 1;
}

double max_traffic(int i, int j)
{ 
	int a,b,s,d,source,destination,sa,da,sb,db,temp,k;
	double TS[n_node]={0},TD[n_node]={0};
	double TS_max[n_node]={0},TD_max[n_node]={0};
	double T_left=0,T_add=0;
	double sum_t=0; // sum of traffic on link e(i,j)

	if(passby_s[i][j][0]==0)
		return 0;
	//else
		//cout<<"num of f for ("<<i<<","<<j<<"): "<<passby_s[i][j][0]<<"\n";

	//------ Initialize max concurrent traffic T_max[s] for each server s 
	//source server
	for(a=1;a<=passby_s[i][j][0];a++)
	{
		s=passby_s[i][j][a];
		for(b=1;b<=passby_s[i][j][0];b++)
		{
			source = passby_s[i][j][b];
			d = passby_d[i][j][b];
			if(source == s)
				TS_max[s] += m[d];
		}

		if( TS_max[s] < m[s])
			TS_max[s] = m[s];
		TS_max[s] = TS_max[s] * B;
	}
	//cout<<"TS_max["<<s<<"] = "<<TS_max[s]<<endl;

	//destination server
	for(a=1;a<=passby_d[i][j][0];a++)
	{
		d=passby_d[i][j][a];
		for(b=1;b<=passby_d[i][j][0];b++)
		{
			destination = passby_d[i][j][b];
			s = passby_s[i][j][b];

			if(destination == d)
				TD_max[d] += m[s];
		}

		if( TD_max[d] < m[d])
			TD_max[d] = m[d];
		TD_max[d] = TD_max[d] * B;
	}
	//cout<<"TD_max["<<d<<"] = "<<TD_max[d]<<endl;


	//--------对passby[i][j]内的pair排序
	for(a=1;a<=passby_s[i][j][0];a++)
		for(b=a;b<=passby_s[i][j][0];b++)
		{
			sa=passby_s[i][j][a]; da=passby_d[i][j][a];
			sb=passby_s[i][j][b]; db=passby_d[i][j][b];

			if(f[sa][da][i][j] > f[sb][db][i][j]) //f-sd-e(i,j)
			{
				//交换a,b位置
				temp = passby_s[i][j][a]; 
				passby_s[i][j][a] = passby_s[i][j][b];
				passby_s[i][j][b] = temp;

				temp = passby_d[i][j][a]; 
				passby_d[i][j][a] = passby_d[i][j][b];
				passby_d[i][j][b] = temp;
			}
		}

	//--------计算traffic,并寻找bottleneck_s
	for(k=1;k<=passby_s[i][j][0];k++)
	{
		s=passby_s[i][j][k];	
		d=passby_d[i][j][k];

		//计算traffic
		T_left = min(TS_max[s]-TS[s], TD_max[d]-TD[d]);
		T_add = min( double(min(m[s],m[d]))*B, T_left);
		//cout<<"T_left="<<T_left<<", T_add="<<T_add<<endl;
		sum_t += T_add*f[s][d][i][j]; //max traffic
		TS[s] += T_add;	
		TD[d] += T_add;

		//寻找bottleneck_s
		if(k==1) //选第1个：因为f是从大到小排过序的
		{
			//再选ms,md中较小的那个
			if(m[s] < m[d])
				bottleneck_s = s; 
			else
				bottleneck_s = d;
		}
	}

	//if(passby_s[i][j][0]!=0)
		//cout<<"e(i,j) max_traffic: "<<sum_t<<"\n\n";
	return sum_t;
}

void Perturbation()
{
	int i,j;

	if(A[0]==0) //只有1个点，没经历故Dijkstra(),routing()
		return;

	//cout<<"\nPerturbaion: max_traffic.\n";
	for(i=0;i<n_node;i++) //for each link e(i,j)
		for(j=0;j<n_node;j++)
			while(max_traffic(i,j)>link[i][j]) 
				m[bottleneck_s]--; 
}

void Data()
{
	int i,j;
	double te,load;

	//cout<<"\nData: final traffic.\n";
	for(i=0;i<n_node;i++) //for each link e(i,j)
		for(j=0;j<n_node;j++)
		{
			if(A[0] == 1)
				te = 0;
			else
				te = max_traffic(i,j);

			if(link[i][j] == 0)
				load = 0;
			else
				load = te/link[i][j];

			max_load = max(max_load,load);
			bw_cost += te; //指需要为网络各link准备的最大带宽总和（硬件资源）
		}
}



bool ALG_loadbalance() 
{
	int s,d,i,j,x,y;

	for(i=0;i<n_server;i++) 
	{
		s = server_set[i]; 
		m[s] = space[s];

		for(j=1;j<=A[0];j++)
		{
			d = A[j];	
			KSP(s,d); //返回path[]，该数组对每个s,d保留到routing()结束，下一轮更新

			/*cout<<"\nFind "<<n_path<<" paths for ("<<s<<","<<d<<")\n"; 
			for(x=0;x<n_path;x++)
			{
				cout<<"  path "<<x<<": ";
				for(y=1;y<=path[x][0];y++)
					cout<<path[x][y]<<" ";
				cout<<endl;
			}
			cout<<endl;*/

			Routing(s,d); //返回f-sd-e
		}
	
		Perturbation();

		if(m[s]>0)
		{
			A[++A[0]] = s;
			sum_m += m[s];

			if(sum_m >= N)
			{
				m[s] -= sum_m - N;	
				Data(); //calculate bn,bw_cost
				return 1;
			}
		}	
	} 
	return 0;
}

void Average(bool su)
{
	int i;

	if(su==1)
	{
		cout<<"\nEmbedding success!\n";
		suc++;
	}
	else
		cout<<"\nEmbedding failure!\n";

	for(i=0;i<n_server;i++)
		cout<<"m["<<server_set[i]<<"] = "<<m[server_set[i]]<<endl;

	cout<<"\nmax_load: "<<max_load<<endl;
	cout<<"total bandwidth: "<<bw_cost<<endl;
	cout<<"success rate: "<<suc<<endl;
}

void main()
{
	bool su;

	Generate();
	Initialize();
	su = ALG_loadbalance(); //Embedding算法
	Average(su);

	system("pause");
}