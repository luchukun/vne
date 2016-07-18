#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

#include<ctime>
using namespace std;
ofstream fout;

typedef IloArray<IloNumVarArray> NumVarMatrix;
typedef IloArray<NumVarMatrix> NumVarCubic;

#define cycle 15
#define num_i 4
#define num_s 10
#define C 20 //缓存空间

//---------拓扑基本参数--------
int space[num_s]; // node space
int link[num_i][num_s]; // link capacity

//------中间参数------
int B,N,times;
double bn[cycle];
int success[cycle]={0};
double bn_best;

int y[num_s];
int Mv[num_s+1][C]; // Mv:第k个node的容量, k=1,..,s_num
int Lv[num_s+1][C*num_s]; // 前k个node能容纳的所有可能值
int Dv[num_s+1][C*num_s][C*num_s]; // 分配在第k个node上的值

void initialize()
{
	int i;
	for(i=0;i<cycle;i++)
		bn[i] = -1;
}

void generate()
{
	int i,j;

	//srand(10*times);

	//----node space分布---
	for(i=0;i<num_s;i++) //只分布leaf层
		space[i]=rand()%5+1;

	//---link capacity 分布---
	for(i=0;i<num_i;i++)
		for(j=0;j<num_s;j++)
			link[i][j]=rand()%100+1;

	cout<<"============= round "<<times<<" =============\n";
	//-----输出拓扑信息-----
	/*cout<<"space:\n";
	for(i=0;i<num_s;i++) //只分布leaf层
		cout<<space[i]<<" ";

	cout<<"\nlinks:\n";
	for(i=0;i<num_i;i++)
	{
		for(j=0;j<num_s;j++)
			cout<<link[i][j]<<" ";
		cout<<"\t";
	}
	cout<<"\n-------------------------\n";*/
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

   IloIntArray m (env,num_s); //VM allocation
   for(s=0;s<num_s;s++)
	   m[s] = y[s];

   /*cout<<"VM: ";
   for(s=0;s<num_s;s++)
	   cout<<m[s]<<" ";
   cout<<endl;*/

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

   IloNumVar bottleneck (env);
   
   //----优化目标----
   model.add(IloMinimize(env,bottleneck));

   //----------cplex列方程---------
   //--bottleneck--
   for(i=0;i<num_i;i++)
	   for(s=0;s<num_s;s++)
		   model.add(  B* IloScalProd( m, x[i][s]) /link[i][s] <= bottleneck);

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
   /*for(i=0;i<num_i;i++)
	   for(s=0;s<num_s;s++)
		   model.add(  B* IloScalProd( m, x[i][s]) <= link[i][s]);*/
   model.add( bottleneck <= 1);

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

	double tar;
	if(cplex.solve())
	{
		tar = cplex.getObjValue();
		if(tar<bn_best) 
			bn_best = tar; 
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
	int i;

	if(bn_best!=2)
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
	int i,s,k,j,e_i,h_i,eh,a,link_sum,flag;

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
			link_sum += link[i][s];//这里可以再修改

		k=0;
		for(j=0;j<=space[s];j++)
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

	bn_best = 2;
	for(j=1;j<=Lv[num_s][0];j++)
		if(Lv[num_s][j]==N)
			alloc(num_s,N);

	if(bn_best != 2)
	{
		success[N-1]++;
		//cout<<"N="<<N<<", bn_best = "<<bn_best<<endl;

		if(bn[N-1] == -1)
			bn[N-1] = bn_best;
		else
			bn[N-1] = ( bn[N-1]*(success[N-1]-1) +bn_best )/success[N-1];
	}
	//else
		//cout<<"N="<<N<<", No solution!\n";

	return 0;
}

void figure()
{
	int i;

	fout.open("E:\\data\\N-bn\\bottleneck3.txt");
	fout.clear();
	fout.close();
	fout.open("E:\\data\\N-suc\\success3.txt");
	fout.clear();
	fout.close();

	fout.open("E:\\data\\N-bn\\bottleneck3.txt",ios::app);
	for(i=0;i<cycle;i++) //(N=1;N<=20;N++)
		fout<<bn[i]<<" ";
	fout.close();

	fout.open("E:\\data\\N-suc\\success3.txt",ios::app);
	for(i=0;i<cycle;i++) //(N=1;N<=20;N++)
		fout<<double(success[i])/(times-1)<<" "; //times-1是因为for循环多加了1
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
	
		B=25;
		for(N=1;N<=cycle;N++) 
			dynamic();
	}

	figure();
	end=clock();
	cout<<"run time: "<<double(end-begin)/CLOCKS_PER_SEC<<" s"<<endl;
	system("pause");
}