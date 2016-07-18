/* 动态系统之 - single path算法 */

/* ---相关统计量---
consume_link,consume_space,res_link; solve_LP()
revenue & cost: update()
*/

#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

#include<ctime>
using namespace std;
ofstream fout;

typedef IloArray<IloNumVarArray> NumVarMatrix;
typedef IloArray<NumVarMatrix> NumVarCubic;

//---常量参数---
#define num_i 4
#define num_s 12
#define N 6 
#define B 5  

//---Poisson arrive---
#define lamda double(1)/1 //到达，每小时到达lamda个，间隔时长 1/lamda
double miu;				  //离开，每小时离开miu个，request时长 1/miu

//---substrate硬件资源---
int space[num_s];				// node space
double link[num_i][num_s];		// link capacity
int res_space[num_s];			// 剩余空间 residual space
double res_link[num_i][num_s];	// 剩余带宽 residual link
#define size_r 10*num_s/N+1		// substrate network最大能容纳request数
#define C 50					// 动态规划中 Mv缓存空间,与space大小有关

//---VDC占用资源---
int consume_space[size_r][num_s]; // 占用 node space
double consume_link[size_r][num_i][num_s]; // 占用 link capacity

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

int y[num_s];
int Mv[num_s+1][C]; // Mv:第k个node的容量, k=1,..,s_num
int Lv[num_s+1][C*num_s]; // 前k个node能容纳的所有可能值
int Dv[num_s+1][C*num_s][C*num_s]; // 分配在第k个node上的值

void generate()
{
	/* 产生硬件资源 */

	int i,j;

	//----node space分布---
	for(i=0;i<num_s;i++)
		space[i]=5;//rand()%10+1;

	//---link capacity 分布---
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
	fout.open("E:\\data\\paper2\\exhaustive\\suc-load.txt");
	fout.clear();
	fout.close();

	fout.open("E:\\data\\paper2\\exhaustive\\R-load.txt");
	fout.clear();
	fout.close();

	fout.open("E:\\data\\paper2\\exhaustive\\C-load.txt");
	fout.clear();
	fout.close();

	fout.open("E:\\data\\paper2\\exhaustive\\RC-load.txt");
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

	for(i=0;i<num_s;i++)
			res_space[i]=space[i];

	for(i=0;i<num_i;i++)
		for(j=0;j<num_s;j++)
			res_link[i][j]=link[i][j];
}

void clear_alloc()
{
	/* 清空缓存分配信息 */

	for(int s=0;s<num_s;s++)
		y[s] = 0;
}

int solve_LP()
{
   //--------基本环境设置--------
   IloEnv env; //环境environment
   IloModel model(env); //model
   int i,s,j,k;

   //-------常量变量设置--------
   IloNumArray a1 (env, num_i); //[1,1,...,1]
   for(i=0;i<num_i;i++)
	   a1[i] = 1;

   IloNumArray a2 (env, num_s); //[1,1,...,1]
   for(s=0;s<num_s;s++)
	   a2[s] = 1;

   IloIntArray m (env,num_s); //VM allocation
   for(s=0;s<num_s;s++)
	   m[s] = y[s];

   //-------待定变量设置------
   NumVarCubic f(env, num_s); //f[s][j][i], f[4][4][2]
   for(s=0;s<num_s;s++)
   {
	   f[s] = NumVarMatrix(env,num_s);
	   for(j=0;j<num_s;j++)
		   f[s][j] = IloNumVarArray(env, num_i, 0, 1 );
   }

   NumVarCubic x(env, num_i); //x[i][s][j], x[2][4][4]
   for(i=0;i<num_i;i++)
   {
	   x[i] = NumVarMatrix(env,num_s);  
	   for(s=0;s<num_s;s++)
		   x[i][s] = IloNumVarArray(env, num_s, 0, IloInfinity );
   }   
   
   IloNumVar bandcost (env); //目标变量

   IloNumVarArray b(env, num_i, 0, IloInfinity); //中间变量
   NumVarMatrix c(env, num_i); 

   for(i=0;i<num_i;i++)
	   c[i] = IloNumVarArray(env, num_s, 0, IloInfinity);

   //===============[分割界: 上用=,定义类型; 下用==,赋予值.]============

    //----优化目标----
   model.add(IloMinimize(env,bandcost));

   //--------cplex列方程-------
   //--bandcost--
   for(i=0;i<num_i;i++)
		    model.add( b[i] == IloScalProd( a2, c[i]) );

   model.add( bandcost == IloScalProd(a1, b) );

   //--f--
   for(i=0;i<num_i;i++)  //f[i][s][s] =0
	   for(s=0;s<num_s;s++)
		   model.add( f[s][s][i] == 0);

   for(i=0;i<num_i;i++)  //f[i][j][k] = f[i][k][j]
	   for(j=0;j<num_s;j++)
		   for(k=j+1;k<num_s;k++)
			   model.add( f[j][k][i] - f[k][j][i] == 0 );

   for(j=0;j<num_s;j++)  //比例和为一
	   for(k=j+1;k<num_s;k++)
		   model.add(  IloScalProd( a1,f[j][k] ) == 1);

   //--x--
   for(i=0;i<num_i;i++)
	   for(s=0;s<num_s;s++)
		   model.add(  B* IloScalProd( m, x[i][s]) <= c[i][s]);

   for(i=0;i<num_i;i++)
	   for(s=0;s<num_s;s++)
		   model.add( c[i][s] <= res_link[i][s]);

	for(i=0;i<num_i;i++)
	   for(s=0;s<num_s;s++)
		   for(j=0;j<num_s;j++)
		   {
			   if(j==s)
				   continue;
			   model.add( x[i][s][s] + x[i][s][j] - f[s][j][i] >= 0 );
		   }
	
	//---------Cplex解方程---------	   
	IloCplex cplex(model); //依model建的cplex
	cplex.setOut(env.getNullStream());

	if(cplex.solve())
	{		
		for(s=0;s<num_s;s++)
		{
			res_space[s] -= m[s];
			consume_space[num_VDC][s] = m[s];
		}

		for(i=0;i<num_i;i++)
			for(s=0;s<num_s;s++)
			{
				res_link[i][s] -= cplex.getValue(c[i][s]); 
				consume_link[num_VDC][i][s] = cplex.getValue(c[i][s]);
			}

		bw_cost = cplex.getObjValue(); //dynamic中要作为判断依据

		env.end();
		return 1;
	}
	
	else
	{
		//cout<<"No Solution!\n\n";
		env.end();
		return 0;
	}
}


void alloc( int s, int n ) //node, amount, for循环
{
	//alloc的递归，是不断打印出不同的组合。在本题中是有一个组合成功即可了。
	int i;

	if(bw_cost != -1)
		return;

	if(n==0 || s==0)
	{
		solve_LP();
		
		for(i=0;i<=s;i++)
			y[i]=0;
		return;
	}

	for(i=1;i<=Dv[s][n][0];i++)
	{
		y[s-1] = Dv[s][n][i];
		alloc(s-1,n-y[s-1]);
	}
}

int dynamic()
{
	int i,s,k,j,e_i,h_i,eh,a,flag;
	double link_sum;

	//----变量初始化----
	for(i=0;i<=num_s;i++)
	{
		for(j=0;j<C;j++)
			Mv[i][j]=0;

		for(j=0;j<C*num_s;j++)
		{
			Lv[i][j]=0;
			for(k=0;k<C*num_s;k++)
				Dv[i][j][k]=0;
		}
	}
	
	//-------储备叶子节点的Mv--------
	for(s=0;s<num_s;s++)
	{		
		link_sum=0;
		
		//check simple bandwidth constrain
		for(i=0;i<num_i;i++)
			link_sum += res_link[i][s];//这里可以再修改

		k=0;
		for(j=0;j<=res_space[s];j++)
		{			
			if( link_sum >= B* min(j,N-j) ) 
				Mv[s+1][++k]=j; 
		}
		Mv[s+1][0]=k; //Mv中可取值的数量
	}
	
	//------------列出所有 m1+m2+...+mn =N 的case-----------
	Lv[0][0]=1; Lv[0][1]=0;

	for(k=1;k<=num_s;k++) // 前 k 个node
	{
		for( e_i=1; e_i<=Mv[k][0]; e_i++ ) //取遍Mv[k]所有值, 从2开始因为：[0]表示数量, [1]放0个无意义
		{
			for( h_i=1; h_i<=Lv[k-1][0]; h_i++ ) //取遍lv[k-1]所有值				
			{			
				flag = 0;

				eh = Mv[k][e_i] + Lv[k-1][h_i]; //计算e+h
				a = ++Dv[k][eh][0];//eh最大size所有mi加起来
				Dv[k][eh][a] = Mv[k][e_i];

				//找是不是已经有该元素了
				for(j=1;j<=Lv[k][0];j++)
					if(Lv[k][j]==eh)
						flag =1;

				if(flag==0)
					Lv[k][++Lv[k][0]] = eh;
			}
		}
	} 

	bw_cost = -1;
	for(j=1;j<=Lv[num_s][0];j++)
		if(Lv[num_s][j]==N)
			alloc(num_s,N);

	if(bw_cost != -1)
	{
		//cout<<"Allocation for VDC "<<num_VDC<<" succeed!\n";
		num_VDC++;
		num_accp++;
	}
	//else
		//cout<<"Allocation Fails!\n";

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

	fout.open("E:\\data\\paper2\\exhaustive\\suc-load.txt",ios::app);
	fout<<suc<<" ";
	fout.close();

	fout.open("E:\\data\\paper2\\exhaustive\\R-load.txt",ios::app);
	fout<<revenue<<" ";
	fout.close();

	fout.open("E:\\data\\paper2\\exhaustive\\C-load.txt",ios::app);
	fout<<cost<<" ";
	fout.close();

	fout.open("E:\\data\\paper2\\exhaustive\\RC-load.txt",ios::app);
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
	generate();  // 设定硬件资源
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
				dynamic();	   //为新来的request分配
				clear_alloc(); //清空缓存分配信息
			}
			update(); //每时刻更新R/C信息
			//print_info();
		}
	
		print_info();
		figure();
	}
	system("pause");
}