/* 动态系统之 - load balance算法 */

/* ---相关统计量---
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

//---规模size---
#define num_i 4
#define num_s 12
#define N 6   
#define B 5  

//---Poisson arrive---
#define lamda double(1)/1 //到达，每小时到达lamda个，间隔时长 1/lamda
double miu;				//离开，每小时离开miu个，request时长 1/miu

//---substrate硬件资源---
int space[num_s];				// node space
double link[num_i][num_s];		// link capacity
int res_space[num_s];			// 剩余空间 residual space
double res_link[num_i][num_s];	// 剩余带宽 residual link
#define size_r 10*num_s/N+1		// substrate network最大能容纳request数

//---VDC占用资源---
int consume_space[size_r][num_s]; // 占用 node space
double consume_link[size_r][num_i][num_s]; // 占用 link capacity
//int begin_time[size_r];

//----分配结果----
int m[num_s];
double f[num_s][num_s][num_i];

//----统计参数----
double suc; //accepted rate
double revenue;
double cost;

//-----中间变量-----
double total_band;
double bw_cost;
int num_arrive; //整个仿真过程中到来的VDC数量
int num_accp;    //被接受的VDC数量
int num_VDC;    //当前网络中储存的VDC数量

int A[num_s+1]; // 选进去的点集
int sum_m; //累计VM数量

void generate()
{
	/* 产生硬件资源 */

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
	/* 打印当前硬件资源 & VDC占用资源 */

	int i,j,k;

	//-----打印 当前硬件资源----
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

	//------打印 VDC占用资源-----
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
	//-----文档清空-----
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
	/* VDC占用资源清零 */

	int i,j,k;
	
	//---统计参数---
	suc = 0;
	revenue = 0;
	cost = 0;
	
	//---辅助变量---
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
	/* 清空缓存分配信息 */
	int i,j,k;
	sum_m = 0;

	for(i=0;i<num_s;i++) //VM allocation
		m[i]=0;

	for(i=0;i<=num_s;i++) 
		A[i]=0;

	for(i=0;i<num_s;i++)  //分配比例
		for(j=0;j<num_s;j++)
			for(k=0;k<num_i;k++)
				f[i][j][k]=-1;

}

int check_single(int s)
{
	int j,i,x;
	double sum_link;

	x=res_space[s];
	//---带宽分约束---
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
	
	//node安排 
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

	//--初始化sort_f[]
	for(k=1;k<A[0];k++)
	{
		//sort_f[]待会儿要排序，不能把s自己放进来
		if(k<j)
			sort_f[k-1] = A[k];
		else
			sort_f[k-1] = A[k+1];
	}
			
	//--排序，sort_f里面按比例f从大到小排列
	for(a=0;a<A[0]-1;a++)
		for(b=A[0]-2;b>a;b--)
			if( f[s][sort_f[b]][i] > f[s][sort_f[b-1]][i] )
			{
				temp = sort_f[b];
				sort_f[b] = sort_f[b-1];
				sort_f[b-1] = temp;
			}


	//--作为判断条件的traffic，取f最大，这样能保证VM数量达上限时，traffic也最大。
	a=0; n = min(m[s], sum_m-m[s]);

	while(a<=A[0]-2 && sum_n< n)
	{
		k = sort_f[a++];
		mm = min( n-sum_n , m[k]); //最坏就 mm = n-0
		
		sum_n += mm;
		sum_t += B* mm *f[s][k][i];
	}

	// 判断条件
	if(sum_t > res_link[i][s])
		return 0;
	
	//------统计参数-----
	//只有final_flag为true才记录，前面的调用都不统计以下参数
	if(final_flag == 1)
	{
		consume_link[num_VDC][i][s] = sum_t;
		res_link[i][s] -= sum_t;
	}
	
	return 1;
}

bool greedy()
{
	/* 采用load balance算法进行VDC嵌入 */

	int s,i,j;

	//---嵌入算法---
	for(s=0;s<num_s;s++) 
	{
		check_single(s); //check single route

		if(m[s]>0)
		{
			A[++A[0]]=s; //满足分约束的节点s先放入集合

			//--check_all：对 link(s,i)检验，超载则对m(s)减一
			for(j=1;j<=A[0];j++)
				for(i=0;i<num_i;i++) 
					while(!check_all(i,j,0))	{
						m[A[j]]--; sum_m--;
					}

			//--VM分配成功，进行VMtraffic routing
			if(sum_m >= N)
			{
				//对VDC分配node资源
				m[s] -= sum_m - N;
				for(i=0;i<num_s;i++)
				{
					consume_space[num_VDC][i] = m[i];
					res_space[i] -= m[i];
				}
				//计算bottleneck 和 bw_cost
				for(j=1;j<=A[0];j++)
					for(i=0;i<num_i;i++) 
						check_all(i,j,1);
				
				//成功统计量+1
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
	/* 每时刻更新revenue,cost */

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
	/* 输出曲线数据 */

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

	//释放网络空间
	for(i=0;i<num_i;i++)
		for(s=0;s<num_s;s++)
			res_link[i][s] += consume_link[x][i][s];

	for(s=0;s<num_s;s++)
		res_space[s] += consume_space[x][s];

	//VDC列表删去VDC x
	
	for(i=0;i<num_i;i++)
		for(s=0;s<num_s;s++)
		{
			for(y=x;y<num_VDC-1;y++) //VDC列表从位置0开始储存
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
	/* 检查有没有VDC到期离开 */

	int i;
	double p2 = 1-exp(-miu); //该时刻 VDC离开的概率

	for(i=0;i<num_VDC;i++)	 //对网络中每个VDC都轮询一遍
	{
		double pp = rand()/double(RAND_MAX); //产生一个0-1之间的随机数
		if( pp < p2) 
		{
			//cout<<"VDC "<<i<<" has leaved\n";
			free_VDC(i); 
			i--; //因为是数组存储，VDC总数减少，后面VDC会整体前移
		}
	}
}

bool arrive()
{
	/* 一个新VDC到达 */

	double p1 = 1-exp(-lamda); //单位时刻内，有新VDC到达的概率
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
	generate();	 // 设定硬件资源
	for(int lasting=1; lasting<=10; lasting++)
	{

		miu = double(1)/lasting;
		cout<<"====== R/C ratio = "<<lamda/miu<<" =======\n";
		clear_VDC(); // VDC占用资源清零

		for(int time=1;time<=100000;time++) //仿真模拟运行xx小时	
		{
			//cout<<"-----time "<<time<<"------\n";
			leave(); //看看有没有VDC离开
			if(arrive())
			{
				greedy();	  //为新来的request分配
				clear_alloc(); //清空缓存分配信息
			}
			update(); //每时刻更新R/C信息
			//print_info();
		}

		print_info(); //打印当前硬件资源 & VDC消耗资源
		figure(); //输出曲线数据
	}
	system("pause");
}