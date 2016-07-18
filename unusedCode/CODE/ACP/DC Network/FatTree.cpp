/* Perturbation algorithm with load balance on general network - ex.1 */
/*
	1. �������У�һ�����飬0λ��ŵ��������С��Ϣ��
	2. ��Ϊֻ�Ƿ���һ�Σ�����ֱ����link,���Ƕ�̬���棬����res_link 
*/

#include<iostream>
#include<fstream>
#include<cmath>
#include<ctime>
using namespace std;
ofstream fout;

//----�����ģ---
#define K 1
#define cycle 15
int B,N=20;

#define n_switch 20
#define n_server 16

#define n_node n_switch+n_server
#define my_MAX 1000

//----Ӳ����Դ-----
int switch_set[n_switch];
int server_set[n_server];

int space[n_node]; // server
int weigh[n_node][n_node]; // link length
double link[n_node][n_node]; // link capacity

//-----�м����-----
int A[n_server+1]; // ѡ��ȥ�ĵ㼯
int sum_m; //�ۼ�VM����
int bn_s; //��λ�ĵ���bottleneck��server
double bw_sum; //�ܴ���
int times;

//Dijstra()
int dis[n_node+1]; //Dijkstra�е㵽s�ľ���
int weigh_t[n_node][n_node]; // Dijksrta�ж�̬weigh
int D[n_switch+1]; //switch set + d
int pre[n_node+1];

int n_D;
int n_path;
int *pathSet[K]; //��ǰs,d֮��path��Ϣ
double capacity[K]; // ��ǰ��path��capacity
int passby_s[n_node][n_node][n_server*n_server];
int passby_d[n_node][n_node][n_server*n_server];

//----������----
int m[n_node];
double f[n_node][n_node][n_node][n_node];

//----ͳ�Ʋ���----
double max_load;
double bw_cost;

double bn_array[cycle+1];
double bw_array[cycle+1];
int suc_array[cycle+1];

void Initialize()
{
	int i;
	for(i=0;i<=cycle;i++)
	{
		bn_array[i] = 0;
		bw_array[i] = 0;
		suc_array[i] = 0;
	}
}
void Generate()
{
	/* generate harware resource */

	/* server_set,space�д�0�ൽn;
	   weigh(����),link(����)�б��Ϊ��� */
	int i,j=0,k;

	// node
	for(i=0;i<n_switch;i++)
		switch_set[i] = j++; //0,1
	for(i=0;i<n_server;i++)
		server_set[i] = j++; //2,3

	// VM space 
	for(i=0;i<n_node;i++) //ֻ�ֲ�leaf��
		space[i] = rand()%5+1;
	
	// link length
	for(i=0;i<n_node;i++)
	{
		for(j=0;j<n_node;j++)
			weigh[i][j] = my_MAX;
		weigh[i][i] = 0;
	}

	for(k=0;k<2;k++)
		for(i=0;i<2;i++)
			for(j=4;j<=10;j+=2)
			{
				weigh[2*k+i][k+j] = 1;
				link[2*k+i][k+j] = rand()%25+25;
			}

	for(i=4;i<=10;i+=2)
		for(j=0;j<=1;j++)
			for(k=i+8;k<=i+9;k++)
			{
				weigh[i+j][k] = 1;
				link[i+j][k] = rand()%25+25;
			}

	for(i=12;i<=19;i++)
		for(j=2*i-4;j<=2*i-3;j++)
		{
			weigh[i][j] = 1;
			link[i][j] = 100;
		}	

	for(i=0;i<n_node;i++)
		for(j=i+1;j<n_node;j++)
		{
			weigh[j][i] = weigh[i][j];
			link[j][i] = link[i][j];
		}

	/*cout<<"weigh matrix:\n";
	for(i=0;i<n_node;i++)
	{
		for(j=0;j<n_node;j++)
		{
			if(weigh[i][j]==my_MAX)
				cout<<"X ";
			else
				cout<<weigh[i][j]<<" ";
		}
		cout<<endl;
	}*/

	/*cout<<"VM space:\n";
	for(i=0;i<n_server;i++) //ֻ�ֲ�leaf��
		cout<<space[server_set[i]]<<" ";
	cout<<"\n";

	cout<<"capacity matrix:\n";
	for(i=0;i<n_node;i++)
	{
		for(j=0;j<n_node;j++)
			cout<<link[i][j]<<" ";
		cout<<endl;
	}*/
}

void print_weigh_t()
{
	int i,j;

	for(i=0;i<n_node;i++)
	{
		cout<<"\t";
		for(j=0;j<n_node;j++)
		{
			if(weigh_t[i][j]==my_MAX)
				cout<<"X ";
			else
				cout<<weigh_t[i][j]<<" ";
		}
		cout<<endl;
	}
}

void clear()
{
	/* Initialize ������ and ͳ�Ʊ��� */
	int s,d,i,j;

	max_load = 0.0;
	bw_cost = 0.0;

	sum_m = 0; 
	bw_sum = 0;

	for(i=0;i<=n_server;i++)
		A[i] = 0;

	for(i=0;i<n_node;i++)
		m[i]=0;

	for(s=0;s<n_node;s++)
		for(d=0;d<n_node;d++)
			for(i=0;i<n_node;i++)
				for(j=0;j<n_node;j++)
					f[s][d][i][j] = 0; //f-sd-e(i,j),ע���˴���s,d�ǵ��±꣬����server�±�

	for(i=0;i<n_node;i++)
		for(j=0;j<n_node;j++)
			weigh_t[i][j] = weigh[i][j];

	for(i=0;i<n_node;i++)
		for(j=0;j<n_node;j++)
		{
			passby_s[i][j][0] = 0;
			passby_d[i][j][0] = 0;
		}

	for(i=0;i<n_node;i++)
		for(j=i+1;j<n_node;j++)
			bw_sum += link[i][j];
}

void Initial(int d, int s)
{
	/* ������ΪDijksrta�㷨��weigh���Initialize */

	int i;

	//Initilize dis[]
	for(i=0;i<n_D;i++)
		dis[D[i]] = weigh_t[s][D[i]];

	//Initialize pre[]
	for(i=0;i<n_D;i++)
		pre[D[i]] = s;
}
void ResetWeigh_t()
{
	int i,j;

	for(i=0;i<n_node;i++)
		for(j=0;j<n_node;j++)
			weigh_t[i][j] = weigh[i][j];
}
int* Dijkstra(int d,int s)
{
	/* s,d֮�����path */
	/* ��¼path�ķ�ʽ����d��s���������� */

	int i,j,k,j_index;
	int *path;
	path = new int[n_node];

	//cout<<"\tDijkstra("<<d<<","<<s<<")\n";
	
	//----------Initialize-----------
	Initial(d,s);

	//---------Dijkstra�㷨----------
	while(n_D > 0)
	{
		/*cout<<"dis[] = ";
		for(i=0;i<n_D;i++)
			cout<<dis[D[i]]<<" ";
		cout<<endl;*/

		//choose min dis[s][j]
		j=d;
		for(i=0;i<n_D;i++) 
			if(dis[D[i]] < dis[j])
				j=D[i];

		//�ҵ����·��������
		if(j==d)
		{
			if(dis[j]==my_MAX)
				return NULL;

			path[1]=d; i=2;
			while(j!=s)
			{
				path[i++]=pre[j];
				j = pre[j];
			}
			path[0]=i-1;

			/*cout<<"\tDijkstra success..\n";
			cout<<"Return path: ";
			for(i=1;i<=path[0];i++)
				cout<<path[i]<<" ";
			cout<<endl;*/

			return path;
		}

		//����·�ɱ�
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

	return NULL;
}

void KSP_my(int d, int s)
{
	/* ��Ϊ�㷨���⣬�Ӻ���ǰ��pair��s,d)��
	   ����ΪDijkstra����path��ǰ����(d->s)���ʴ�ǰ����ִ��KSP�㷨. */

	int i,j,m,n,spurNode,rootnode,flag;
	int k=0,size_B=0,min_length=100;
	int *spurPath; 
	int pathBuffer[1000][n_node];//temp buffer

	//---D[], weigh_t
	for(i=0;i<n_switch;i++)
		D[i] = switch_set[i];
	D[n_switch] = d;
	n_D = n_switch+1;
	ResetWeigh_t();
	//print_weigh_t();

	pathSet[0] = Dijkstra(d,s);
	n_path=1;

	/*cout<<"  Path 1: ";
	for(i=1;i<=pathSet[0][0];i++)
		cout<<pathSet[0][i]<<" ";
	cout<<endl;*/

	while( k < K-1 )
	{
		for(i=1;i<pathSet[k][0];i++) //��path�е�ÿ��������ɾȥ
		{
			spurNode = pathSet[k][i];

			//cout<<"\ti="<<i<<", spurnode="<<spurNode<<endl;

			//---��ʼ��D[]: 
			for(j=0;j<n_switch;j++)
				D[j] = switch_set[j];
			D[n_switch] = d;
			n_D = n_switch+1;
			ResetWeigh_t();

			//---ɾ��rootpath�ϵĵ�
			for(j=2;j<=i;j++) //j==i��Ϊd�㣬ҲҪɾȥ 
			{
				rootnode = pathSet[k][j-1];
				for(m=0;m<n_D;m++)
					if(D[m]==rootnode)
					{
						for(n=m;n<n_D-1;n++)
							D[n] = D[n+1];
						n_D--;
					}
			}

			//---ɾȥ��rootpath��, pathSet������·��ʹ�ù��ı�
			for(j=0;j<=k;j++) 
			{
				flag = 0;
				for( m=1;m<=i;m++)
					if( pathSet[j][m] != pathSet[k][m] )
						flag = 1;

				if(flag == 0)
				{
					weigh_t[pathSet[j][i]][pathSet[j][i+1]] = my_MAX;
					weigh_t[pathSet[j][i+1]][pathSet[j][i]] = my_MAX;
				}
			}

			//---ɾȥ��rootpath��, Buffer������·��ʹ�ù��ı�
			for(j=0;j<size_B;j++) 
			{
				flag = 0;
				for( m=1;m<=i;m++)
					if( pathBuffer[j][m] != pathSet[k][m] )
						flag = 1;

				if(flag == 0)
				{
					weigh_t[pathBuffer[j][i]][pathBuffer[j][i+1]] = my_MAX;
					weigh_t[pathBuffer[j][i+1]][pathBuffer[j][i]] = my_MAX;
				}
			}

			/*cout<<"\tD[] = ";
			for(j=0;j<n_D;j++)
				cout<<D[j]<<" ";
			cout<<endl;

			cout<<"\t----------\n";
			print_weigh_t();
			cout<<"\t----------\n";*/

			//----After initialize D[] and weigh_t
			spurPath = Dijkstra( pathSet[k][i], s); 
			if(spurPath==NULL)
				continue;

			/*cout<<"\tspurPath: ";
			for(j=1;j<=spurPath[0];j++)
				cout<<spurPath[j]<<" ";
			cout<<endl;*/

			//rootPath + spurPath
			if(i-1+spurPath[0] == pathSet[0][0])
			{
				//rootPath + spurPath
				for(j=1;j<i;j++)
					pathBuffer[size_B][j] = pathSet[k][j];
				for(j=1;j<=spurPath[0];j++)
					pathBuffer[size_B][i-1+j] = spurPath[j];
				pathBuffer[size_B][0] = i - 1 + spurPath[0];
				size_B++;
			}

			/*cout<<"  B[]:\n";
			for(m=0;m<size_B;m++)
			{
				cout<<"\t";
				for(n=1;n<=pathBuffer[m][0];n++)
					cout<<pathBuffer[m][n]<<" ";
				cout<<endl;
			}
			cout<<"-------------------------\n";*/

			//ResetWeigh_t();
		}

		if(size_B==0)
			return;

		//����A
		k++; 
		pathSet[k] = new int[n_node];
		for(i=0;i<=pathSet[0][0];i++)
			pathSet[k][i] = pathBuffer[0][i];

		/*cout<<"  Path "<<k+1<<": ";
		for(i=1;i<=pathSet[k][0];i++)
			cout<<pathSet[k][i]<<" ";
		cout<<endl;*/
		n_path++;
		
		//��B��ɾȥ
		for(i=0;i<size_B;i++)
			for(j=0;j<n_node;j++)
				pathBuffer[i][j] = pathBuffer[i+1][j];
		size_B--;
	}
}

int Routing(int d,int s)
{
	int j,i,k,p,q,n,flag;
	double sum_capacity = 0;

	//cout<<"Routing:\n";
	//--- ����path capacity ---
	for(i=0;i<n_path;i++) 
	{
		capacity[i] = my_MAX;
		for(j=1;j<pathSet[i][0];j++)
		{
			q=pathSet[i][j];p=pathSet[i][j+1]; //link e(p,q)
			if( link[p][q] < capacity[i] )
 				capacity[i] = link[p][q];
		}
		sum_capacity += capacity[i];
	}

	//--- route with load balance ---
	for(i=0;i<n_path;i++)  //ÿ��path
	{
		for(j=1;j<pathSet[i][0];j++)  //path��ÿ��link
		{
			q=pathSet[i][j];p=pathSet[i][j+1]; //link e(p,q)

			//f-sd��e(p,q)�ϵķ������
			f[s][d][p][q] += capacity[i]/sum_capacity ; //++++++++++�ɽ����·����Ҫ����������+++++++++
			f[d][s][p][q] = f[s][d][p][q]; // p->q�ǿ��Ƿ����

			//��û�д��棬����e(p,q)�ϴ洢(s,d)��Ϣ
			flag=0;
			for(k=1;k<=passby_s[p][q][0];k++)
				if(passby_s[p][q][k]==s && passby_d[p][q][k]==d)
					flag=1;

			if(flag==0)
			{
				passby_s[p][q][0]++;   
				passby_d[p][q][0]++;
				n = passby_s[p][q][0];

				passby_s[p][q][n] = s; 
				passby_d[p][q][n] = d;
			}

			//cout<<"passby["<<p<<"]["<<q<<"]: ("<<s<<","<<d<<")\n";
		}
		//cout<<"  f on path "<<i+1<<": "<<capacity[i]/sum_capacity<<endl;
	}

	return 1;
}

double max_traffic(int i, int j)
{ 
	int a,b,s,d,sa,da,sb,db,temp,k,bn;
	double TS[n_node]={0},TD[n_node]={0};
	double TS_max[n_node]={0},TD_max[n_node]={0};
	double T_left=0,T_add=0,T_max=0;
	double sum_t=0; // sum of traffic on link e(i,j)

	if(passby_s[i][j][0]==0)
		return 0;

	//---- Initialize max concurrent traffic T_max[s] for each server s 
	//TS_max = min(m[s],sum(m[d]))
	for(a=1;a<=passby_s[i][j][0];a++)
	{
		s = passby_s[i][j][a];
		d = passby_d[i][j][a];
		TS_max[s] += m[d];
	}

	if( TS_max[s] > m[s])
		TS_max[s] = m[s];
	TS_max[s] = TS_max[s] * B;

	//TD_max = min(m[d],sum(m[s]))
	for(a=1;a<=passby_d[i][j][0];a++)
	{
		d = passby_d[i][j][a];
		s = passby_s[i][j][a];
		TD_max[d] += m[s];

		if( TD_max[d] > m[d])
			TD_max[d] = m[d];
		TD_max[d] = TD_max[d] * B;
	}

	//--------��passby[i][j]�ڵ�pair����
	for(a=1;a<=passby_s[i][j][0];a++)
		for(b=a;b<=passby_s[i][j][0];b++)
		{
			sa=passby_s[i][j][a]; da=passby_d[i][j][a];
			sb=passby_s[i][j][b]; db=passby_d[i][j][b];

			if(f[sa][da][i][j] > f[sb][db][i][j]) //f-sd-e(i,j)
			{
				//����a,bλ��
				temp = passby_s[i][j][a]; 
				passby_s[i][j][a] = passby_s[i][j][b];
				passby_s[i][j][b] = temp;

				temp = passby_d[i][j][a]; 
				passby_d[i][j][a] = passby_d[i][j][b];
				passby_d[i][j][b] = temp;
			}
		}

	//--------����traffic,��Ѱ��bn_s
	for(k=1;k<=passby_s[i][j][0];k++)
	{
		s=passby_s[i][j][k];	
		d=passby_d[i][j][k];

		//����traffic
		T_left = min(TS_max[s]-TS[s], TD_max[d]-TD[d]);
		T_add = min( double(min(m[s],m[d]))*B, T_left);
		sum_t += T_add*f[s][d][i][j]; //max traffic
		TS[s] += T_add;	
		TD[d] += T_add;

		//Ѱ��bn_s
		//if(k==1) //ѡ���f

		if(T_add>T_max)
			bn=k;
	}

	if(rand()/double(RAND_MAX) < 0.5)
		bn_s = passby_s[i][j][bn]; 
	else
		bn_s = passby_d[i][j][bn];

	return sum_t;
}

void Perturbation()
{
	int i,j;

	if(A[0]==0) //ֻ��1���㣬û������Dijkstra(),routing()
		return;

	//cout<<"Perturbation:\n";
	for(i=0;i<n_node;i++) //for each link e(i,j)
		for(j=0;j<n_node;j++)
			while(max_traffic(i,j)>link[i][j]) 
			{
				//cout<<"link("<<i<<","<<j<<"), m["<<bn_s<<"]-=1, to"<<m[bn_s]-1<<endl;
				m[bn_s]--; 
				sum_m--;
			}
}

void Data()
{
	/*�������������ͳ������*/
	int i,j;
	double te,load;

	//cout<<"---traffic---\n";
	for(i=0;i<n_node;i++) //for each link e(i,j)
	{
		for(j=i+1;j<n_node;j++)
		{
			if(A[0] == 1)
				te = 0;
			else
				te = max( max_traffic(i,j), max_traffic(j,i) );

			if(link[i][j] == 0)
				load = 0;
			else
				load = te/link[i][j];

			//if(load>0)
			//	cout<<"e("<<i<<","<<j<<") = "<<load<<endl;

			max_load = max(max_load,load);
			bw_cost += te; //ָ��ҪΪ�����link׼�����������ܺͣ�Ӳ����Դ��
		}
		//cout<<endl;
	}
}



bool ALG_loadbalance() 
{
	int s,d,i,j;

	for(i=0;i<n_server;i++) 
	{
		s = server_set[i]; 

		//cout<<"\n=====try server "<<s<<" =====\n";
		m[s] = space[s];
		sum_m += m[s];

		for(j=1;j<=A[0];j++)
		{
			d = A[j];	
			KSP_my(d,s); //����path[]������ֱ����һ�ָ���
			Routing(d,s); //����f-sd-e
		}

		Perturbation();

		if(m[s]>0)
		{
			A[++A[0]] = s;
			//sum_m += m[s];

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
	if(su==1)
	{
		//cout<<"Embedding success!\n\n";
		/*cout<<"VM allocation:\n";
		for(int i=0;i<n_server;i++) //ֻ�ֲ�leaf��
			cout<<m[server_set[i]]<<" ";
		cout<<"\n\n";*/

		suc_array[B]++;
		bn_array[B] = ( bn_array[B]*(suc_array[B]-1) + max_load )/suc_array[B]; 
		bw_array[B] = ( bw_array[B]*(suc_array[B]-1) + bw_cost/bw_sum )/suc_array[B];
	}
		
	//else
		//cout<<"Embedding failure!\n\n";
}

void figure()
{
	int i;

	/*fout.open("E:\\data\\KSP\\Fat-tree\\bn-B\\4.txt");
	fout.clear();
	fout.close();
	fout.open("E:\\data\\KSP\\Fat-tree\\bw-B\\1.txt");
	fout.clear();
	fout.close();*/
	fout.open("E:\\data\\KSP\\Fat-tree\\suc-B\\1.txt");
	fout.clear();
	fout.close();

	/*fout.open("E:\\data\\KSP\\Fat-tree\\bn-B\\4.txt",ios::app);
	for(i=1;i<=cycle;i++)
		fout<<bn_array[i]<<" ";
	fout.close();

	fout.open("E:\\data\\KSP\\Fat-tree\\bw-B\\1.txt",ios::app);
	for(i=1;i<=cycle;i++) 
		fout<<bw_array[i]<<" ";
	fout.close();*/

	fout.open("E:\\data\\KSP\\Fat-tree\\suc-B\\1.txt",ios::app);
	for(i=1;i<=cycle;i++)
		fout<<double(suc_array[i])/(times-1)<<" "; //times��������6
	fout.close();
}
void main()
{
	bool su;

	Initialize();
	for(times=1;times<=100;times++)
	{
		if(times%10==0)
			cout<<"=========="<<times/10<<"===========\n";
		Generate();

		for(B=1;B<=cycle;B++)
		{
			clear();
			su = ALG_loadbalance(); //Embedding�㷨
			Average(su);
		}
	}
	figure();

	system("pause");
}