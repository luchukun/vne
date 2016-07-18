/* ��̬ϵͳ֮ - single path�㷨 */

/* ---���ͳ����---
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

//---��������---
#define num_i 4
#define num_s 12
#define N 6 
#define B 5  

//---Poisson arrive---
#define lamda double(1)/1 //���ÿСʱ����lamda�������ʱ�� 1/lamda
double miu;				  //�뿪��ÿСʱ�뿪miu����requestʱ�� 1/miu

//---substrateӲ����Դ---
int space[num_s];				// node space
double link[num_i][num_s];		// link capacity
int res_space[num_s];			// ʣ��ռ� residual space
double res_link[num_i][num_s];	// ʣ����� residual link
#define size_r 10*num_s/N+1		// substrate network���������request��
#define C 50					// ��̬�滮�� Mv����ռ�,��space��С�й�

//---VDCռ����Դ---
int consume_space[size_r][num_s]; // ռ�� node space
double consume_link[size_r][num_i][num_s]; // ռ�� link capacity

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

int y[num_s];
int Mv[num_s+1][C]; // Mv:��k��node������, k=1,..,s_num
int Lv[num_s+1][C*num_s]; // ǰk��node�����ɵ����п���ֵ
int Dv[num_s+1][C*num_s][C*num_s]; // �����ڵ�k��node�ϵ�ֵ

void generate()
{
	/* ����Ӳ����Դ */

	int i,j;

	//----node space�ֲ�---
	for(i=0;i<num_s;i++)
		space[i]=5;//rand()%10+1;

	//---link capacity �ֲ�---
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

	for(i=0;i<num_s;i++)
			res_space[i]=space[i];

	for(i=0;i<num_i;i++)
		for(j=0;j<num_s;j++)
			res_link[i][j]=link[i][j];
}

void clear_alloc()
{
	/* ��ջ��������Ϣ */

	for(int s=0;s<num_s;s++)
		y[s] = 0;
}

int solve_LP()
{
   //--------������������--------
   IloEnv env; //����environment
   IloModel model(env); //model
   int i,s,j,k;

   //-------������������--------
   IloNumArray a1 (env, num_i); //[1,1,...,1]
   for(i=0;i<num_i;i++)
	   a1[i] = 1;

   IloNumArray a2 (env, num_s); //[1,1,...,1]
   for(s=0;s<num_s;s++)
	   a2[s] = 1;

   IloIntArray m (env,num_s); //VM allocation
   for(s=0;s<num_s;s++)
	   m[s] = y[s];

   //-------������������------
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
   
   IloNumVar bandcost (env); //Ŀ�����

   IloNumVarArray b(env, num_i, 0, IloInfinity); //�м����
   NumVarMatrix c(env, num_i); 

   for(i=0;i<num_i;i++)
	   c[i] = IloNumVarArray(env, num_s, 0, IloInfinity);

   //===============[�ָ��: ����=,��������; ����==,����ֵ.]============

    //----�Ż�Ŀ��----
   model.add(IloMinimize(env,bandcost));

   //--------cplex�з���-------
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

   for(j=0;j<num_s;j++)  //������Ϊһ
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
	
	//---------Cplex�ⷽ��---------	   
	IloCplex cplex(model); //��model����cplex
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

		bw_cost = cplex.getObjValue(); //dynamic��Ҫ��Ϊ�ж�����

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


void alloc( int s, int n ) //node, amount, forѭ��
{
	//alloc�ĵݹ飬�ǲ��ϴ�ӡ����ͬ����ϡ��ڱ���������һ����ϳɹ������ˡ�
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

	//----������ʼ��----
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
	
	//-------����Ҷ�ӽڵ��Mv--------
	for(s=0;s<num_s;s++)
	{		
		link_sum=0;
		
		//check simple bandwidth constrain
		for(i=0;i<num_i;i++)
			link_sum += res_link[i][s];//����������޸�

		k=0;
		for(j=0;j<=res_space[s];j++)
		{			
			if( link_sum >= B* min(j,N-j) ) 
				Mv[s+1][++k]=j; 
		}
		Mv[s+1][0]=k; //Mv�п�ȡֵ������
	}
	
	//------------�г����� m1+m2+...+mn =N ��case-----------
	Lv[0][0]=1; Lv[0][1]=0;

	for(k=1;k<=num_s;k++) // ǰ k ��node
	{
		for( e_i=1; e_i<=Mv[k][0]; e_i++ ) //ȡ��Mv[k]����ֵ, ��2��ʼ��Ϊ��[0]��ʾ����, [1]��0��������
		{
			for( h_i=1; h_i<=Lv[k-1][0]; h_i++ ) //ȡ��lv[k-1]����ֵ				
			{			
				flag = 0;

				eh = Mv[k][e_i] + Lv[k-1][h_i]; //����e+h
				a = ++Dv[k][eh][0];//eh���size����mi������
				Dv[k][eh][a] = Mv[k][e_i];

				//���ǲ����Ѿ��и�Ԫ����
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
	generate();  // �趨Ӳ����Դ
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
				dynamic();	   //Ϊ������request����
				clear_alloc(); //��ջ��������Ϣ
			}
			update(); //ÿʱ�̸���R/C��Ϣ
			//print_info();
		}
	
		print_info();
		figure();
	}
	system("pause");
}