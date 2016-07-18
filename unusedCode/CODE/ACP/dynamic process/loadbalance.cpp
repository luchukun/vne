/* ��̬ϵͳ֮ - load balance�㷨 */

/* ---���ͳ����---
consume_link:   check_all()
consume_space:  greedy()
revenue & cost: update()
*/

#include<iostream>
#include<fstream>
#include<cmath>
#include<ctime>
using namespace std;
ofstream fout;

//---��ģsize---
#define num_i 4
#define num_s 12
#define N 6   
#define B 5  

//---Poisson arrive---
#define lamda double(1)/1 //���ÿСʱ����lamda�������ʱ�� 1/lamda
double miu;				//�뿪��ÿСʱ�뿪miu����requestʱ�� 1/miu

//---substrateӲ����Դ---
int space[num_s];				// node space
double link[num_i][num_s];		// link capacity
int res_space[num_s];			// ʣ��ռ� residual space
double res_link[num_i][num_s];	// ʣ����� residual link
#define size_r 10*num_s/N+1		// substrate network���������request��

//---VDCռ����Դ---
int consume_space[size_r][num_s]; // ռ�� node space
double consume_link[size_r][num_i][num_s]; // ռ�� link capacity
//int begin_time[size_r];

//----������----
int m[num_s];
double f[num_s][num_s][num_i];

//----ͳ�Ʋ���----
double suc; //accepted rate
double revenue;
double cost;

//-----�м����-----
double total_band;
double bw_cost;
int num_arrive; //������������е�����VDC����
int num_accp;    //�����ܵ�VDC����
int num_VDC;    //��ǰ�����д����VDC����

int A[num_s+1]; // ѡ��ȥ�ĵ㼯
int sum_m; //�ۼ�VM����

void generate()
{
	/* ����Ӳ����Դ */

	int i,j;
	total_band = 0;

	//----node space---
	for(i=0;i<num_s;i++)
		space[i]=5;//rand()%5+1;

	//---link capacity---
	for(i=0;i<num_i;i++)
		for(j=0;j<num_s;j++)
			link[i][j]=rand()%6+5;//10;
}

void print_info()
{
	/* ��ӡ��ǰӲ����Դ & VDCռ����Դ */

	int i,j,k;

	//-----��ӡ ��ǰӲ����Դ----
	cout<<"\ncurrent resource:";
	cout<<"\n   node: ";
	for(i=0;i<num_s;i++)
		cout<<res_space[i]<<" ";

	cout<<"\n   link: ";
	for(i=0;i<num_i;i++)
	{
		for(j=0;j<num_s;j++)
			cout<<res_link[i][j]<<" ";
		cout<<" ";
	}

	//------��ӡ VDCռ����Դ-----
	for(k=0;k<num_VDC;k++)
	{
		cout<<"\nVDC "<<k<<" resource:";
		cout<<"\n   node: ";
		for(i=0;i<num_s;i++)
			cout<<consume_space[k][i]<<" ";

		cout<<"\n   link: ";
		for(i=0;i<num_i;i++)
		{
			for(j=0;j<num_s;j++)
				cout<<consume_link[k][i][j]<<" ";
			cout<<" ";
		}
	}
	cout<<"\n\n";
}

void clear_file()
{
	//-----�ĵ����-----
	fout.open("E:\\data\\paper2\\loadbalance\\suc-load.txt");
	fout.clear();
	fout.close();

	fout.open("E:\\data\\paper2\\loadbalance\\R-load.txt");
	fout.clear();
	fout.close();

	fout.open("E:\\data\\paper2\\loadbalance\\C-load.txt");
	fout.clear();
	fout.close();

	fout.open("E:\\data\\paper2\\loadbalance\\RC-load.txt");
	fout.clear();
	fout.close();
}
void clear_VDC()
{
	/* VDCռ����Դ���� */

	int i,j,k;
	
	//---ͳ�Ʋ���---
	suc = 0;
	revenue = 0;
	cost = 0;
	
	//---��������---
	num_arrive = 0;
	num_accp = 0;
	num_VDC = 0; // num_VDC <= size_r

	for(i=0;i<size_r;i++) 
		for(j=0;j<num_i;j++)
			for(k=0;k<num_s;k++)
				consume_link[i][j][k]=0;

	for(i=0;i<size_r;i++) 
		for(j=0;j<num_s;j++)
			consume_space[i][j]=0;

	/*for(i=0;i<size_r;i++)
		begin_time[i]=0;*/

	for(i=0;i<num_s;i++)
		res_space[i]=space[i];

	for(i=0;i<num_i;i++)
		for(j=0;j<num_s;j++)
			res_link[i][j]=link[i][j];
}

void clear_alloc()
{
	/* ��ջ��������Ϣ */
	int i,j,k;
	sum_m = 0;

	for(i=0;i<num_s;i++) //VM allocation
		m[i]=0;

	for(i=0;i<=num_s;i++) 
		A[i]=0;

	for(i=0;i<num_s;i++)  //�������
		for(j=0;j<num_s;j++)
			for(k=0;k<num_i;k++)
				f[i][j][k]=-1;

}

int check_single(int s)
{
	int j,i,x;
	double sum_link;

	x=res_space[s];
	//---�����Լ��---
	for(j=1;j<=A[0];j++)
	{
		sum_link = 0;
		for(i=0;i<num_i;i++)
			sum_link += min( res_link[i][s], res_link[i][A[j]] );

		while(sum_link < B* min(x,m[A[j]]))
			x--;

		//---load balance---
		for(i=0;i<num_i;i++) 
		{
			f[s][A[j]][i] = double( min(res_link[i][s], res_link[i][A[j]]))/sum_link ;
			f[A[j]][s][i] = f[s][A[j]][i];
		}
	}
	
	//node���� 
	m[s]=x;
	sum_m += x;

	return 1;
}

int check_all(int i, int j, bool final_flag)
{ 
	int a,b,s,k,sum_n,n,temp,mm;
	int sort_f[num_s-1] = {0};
	double sum_t;
	
	s = A[j];
	sum_n=0; // # of VM
	sum_t=0; // # of traffic

	//--��ʼ��sort_f[]
	for(k=1;k<A[0];k++)
	{
		//sort_f[]�����Ҫ���򣬲��ܰ�s�Լ��Ž���
		if(k<j)
			sort_f[k-1] = A[k];
		else
			sort_f[k-1] = A[k+1];
	}
			
	//--����sort_f���水����f�Ӵ�С����
	for(a=0;a<A[0]-1;a++)
		for(b=A[0]-2;b>a;b--)
			if( f[s][sort_f[b]][i] > f[s][sort_f[b-1]][i] )
			{
				temp = sort_f[b];
				sort_f[b] = sort_f[b-1];
				sort_f[b-1] = temp;
			}


	//--��Ϊ�ж�������traffic��ȡf��������ܱ�֤VM����������ʱ��trafficҲ���
	a=0; n = min(m[s], sum_m-m[s]);

	while(a<=A[0]-2 && sum_n< n)
	{
		k = sort_f[a++];
		mm = min( n-sum_n , m[k]); //��� mm = n-0
		
		sum_n += mm;
		sum_t += B* mm *f[s][k][i];
	}

	// �ж�����
	if(sum_t > res_link[i][s])
		return 0;
	
	//------ͳ�Ʋ���-----
	//ֻ��final_flagΪtrue�ż�¼��ǰ��ĵ��ö���ͳ�����²���
	if(final_flag == 1)
	{
		consume_link[num_VDC][i][s] = sum_t;
		res_link[i][s] -= sum_t;
	}
	
	return 1;
}

bool greedy()
{
	/* ����load balance�㷨����VDCǶ�� */

	int s,i,j;

	//---Ƕ���㷨---
	for(s=0;s<num_s;s++) 
	{
		check_single(s); //check single route

		if(m[s]>0)
		{
			A[++A[0]]=s; //�����Լ���Ľڵ�s�ȷ��뼯��

			//--check_all���� link(s,i)���飬�������m(s)��һ
			for(j=1;j<=A[0];j++)
				for(i=0;i<num_i;i++) 
					while(!check_all(i,j,0))	{
						m[A[j]]--; sum_m--;
					}

			//--VM����ɹ�������VMtraffic routing
			if(sum_m >= N)
			{
				//��VDC����node��Դ
				m[s] -= sum_m - N;
				for(i=0;i<num_s;i++)
				{
					consume_space[num_VDC][i] = m[i];
					res_space[i] -= m[i];
				}
				//����bottleneck �� bw_cost
				for(j=1;j<=A[0];j++)
					for(i=0;i<num_i;i++) 
						check_all(i,j,1);
				
				//�ɹ�ͳ����+1
				num_accp++;
				num_VDC++;
				return 1;
			}
		}
		
	} 
	//cout<<"Allocation failed!\n";
	return 0;
}


void update()
{
	/* ÿʱ�̸���revenue,cost */

	int i,s;

	bw_cost = 0;
	for(i=0;i<num_i;i++)
		for(s=0;s<num_s;s++)
			bw_cost += link[i][s]-res_link[i][s];

	revenue += num_VDC*N + num_VDC*N*B; //kv:kb = 1:1
	cost += num_VDC*N + bw_cost;
}

void figure()
{
	/* ����������� */

	double suc = double(num_accp)/num_arrive;
	double RC = revenue/cost;

	fout.open("E:\\data\\paper2\\loadbalance\\suc-load.txt",ios::app);
	fout<<suc<<" ";
	fout.close();

	fout.open("E:\\data\\paper2\\loadbalance\\R-load.txt",ios::app);
	fout<<revenue<<" ";
	fout.close();

	fout.open("E:\\data\\paper2\\loadbalance\\C-load.txt",ios::app);
	fout<<cost<<" ";
	fout.close();

	fout.open("E:\\data\\paper2\\loadbalance\\RC-load.txt",ios::app);
	fout<<RC<<" ";
	fout.close();
}

void free_VDC(int x)
{
	int i,s,y;

	//�ͷ�����ռ�
	for(i=0;i<num_i;i++)
		for(s=0;s<num_s;s++)
			res_link[i][s] += consume_link[x][i][s];

	for(s=0;s<num_s;s++)
		res_space[s] += consume_space[x][s];

	//VDC�б�ɾȥVDC x
	
	for(i=0;i<num_i;i++)
		for(s=0;s<num_s;s++)
		{
			for(y=x;y<num_VDC-1;y++) //VDC�б��λ��0��ʼ����
				consume_link[y][i][s] = consume_link[y+1][i][s];
			consume_link[num_VDC-1][i][s] = 0;
		}

	for(s=0;s<num_s;s++)
	{
		for(y=x;y<num_VDC-1;y++) 
			consume_space[y][s] = consume_space[y+1][s];
		consume_space[num_VDC-1][s] = 0;
	}

	/*for(y=x;y<num_VDC-1;y++) 
		begin_time[y] = begin_time[y+1];
	begin_time[num_VDC-1] = 0;*/

	num_VDC--;
}

void leave()
{
	/* �����û��VDC�����뿪 */

	int i;
	double p2 = 1-exp(-miu); //��ʱ�� VDC�뿪�ĸ���

	for(i=0;i<num_VDC;i++)	 //��������ÿ��VDC����ѯһ��
	{
		double pp = rand()/double(RAND_MAX); //����һ��0-1֮��������
		if( pp < p2) 
		{
			//cout<<"VDC "<<i<<" has leaved\n";
			free_VDC(i); 
			i--; //��Ϊ������洢��VDC�������٣�����VDC������ǰ��
		}
	}
}

bool arrive()
{
	/* һ����VDC���� */

	double p1 = 1-exp(-lamda); //��λʱ���ڣ�����VDC����ĸ���
	double pp = rand()/double(RAND_MAX);

	if( pp < p1)
	{
		//cout<<"New VDC arrived.\n";
		num_arrive++;
		return 1;
	}
	else
		return 0;
}

void main()
{
	clear_file();
	generate();	 // �趨Ӳ����Դ
	for(int lasting=1; lasting<=10; lasting++)
	{

		miu = double(1)/lasting;
		cout<<"====== R/C ratio = "<<lamda/miu<<" =======\n";
		clear_VDC(); // VDCռ����Դ����

		for(int time=1;time<=100000;time++) //����ģ������xxСʱ	
		{
			//cout<<"-----time "<<time<<"------\n";
			leave(); //������û��VDC�뿪
			if(arrive())
			{
				greedy();	  //Ϊ������request����
				clear_alloc(); //��ջ��������Ϣ
			}
			update(); //ÿʱ�̸���R/C��Ϣ
			//print_info();
		}

		print_info(); //��ӡ��ǰӲ����Դ & VDC������Դ
		figure(); //�����������
	}
	system("pause");
}