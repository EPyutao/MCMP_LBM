//Author: Yu Tao ,Huazhong University of Science and Technology, Wuhan, China 2016.6.1
//Email: lucciola2015x@gmail.com
//battery porous media
#include <iostream>
#include<cmath>
#include<cstdlib>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>
#include<ctime>

using namespace std;

const int Q=9;
const int NX=224;
const int NY=124;

int e[Q][2]={{0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}};
int r[Q]={0,3,4,1,2,7,8,5,6};
double w[Q]={4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36};
int flag[NX+1][NY+1];


double au[NX+1][NY+1][2],bu[NX+1][NY+1][2],arho[NX+1][NY+1],brho[NX+1][NY+1],af[NX+1][NY+1][Q],bf[NX+1][NY+1][Q],aF[NX+1][NY+1][Q],bF[NX+1][NY+1][Q],u[NX+1][NY+1][2];

double atu = 1,btu = 1;
double gff = 0.49;
double gfsa = -0.03, gfsb = -0.03;
double feq(int k,double rho,double u[2]);
double xforceFluid(int i, int j, int component);
double yforceFluid(int i, int j, int component);
double xforceSolid(int i, int j, int component);
double yforceSolid(int i, int j, int component);
double aphi(int i, int j);
double bphi(int i, int j);
int s(int i, int j);

int n,k,i,j,ip,jp;

void init();
void evolution();
void marco();
void error();
void output(int m );

int main()
{
    init();
    for(n=1;n<=30;n++)
    {
        evolution();
        marco();
        if(n%5==0)
        {
          output(n);
        }
        cout<<n<<endl;
        cout<<"velocity_45,45_x "<<u[45][45][0]<<endl;
        cout<<"velocity_45,45_y "<<u[45][45][1]<<endl;
        cout<<arho[46][45]<<endl;
        cout<<endl;
    }
    return 0;
}

void init()
{
    //flag define
    ifstream porous;
    porous.open("test2.txt");
    for(j=0;j<125;j++)
    {
        for(i=0;i<225;i++)
        {
            porous >> flag[i][j];
            //flag == 1 => porous structure
        }
    }
    //initial rho
    for(i=0;i<=NX;i++)
    {
        for(j=0;j<=NY;j++)
        {
            if(flag[i][j]==1)
            {
                arho[i][j] = 0;
                brho[i][j] = 0;
            }
            else{
                arho[i][j] = 0;
                brho[i][j] = 1;
            }
            au[i][j][0] = 0;
            au[i][j][1] = 0;
            bu[i][j][0] = 0;
            bu[i][j][1] = 0;
            u[i][j][0] = 0;
            u[i][j][1] = 0;
        }
    }
    for(i=50;i<=100;i++)
    {
        for(j=30;j<=80;j++)
        {
            if(flag[i][j]==1)
            {
                arho[i][j] = 0;
                brho[i][j] = 0;
            }
            else{
                arho[i][j] = 1;
                brho[i][j] = 0;
            }
        }
    }
    //initial the distribution function
    for(i=0;i<=NX;i++)
    {
        for(j=0;j<=NY;j++)
        {
            for(k=0;k<Q;k++)
                {
                    af[i][j][k] = feq(k,arho[i][j],au[i][j]);
                    bf[i][j][k] = feq(k,brho[i][j],bu[i][j]);
                }
        }
    }
    ostringstream name;
	name<<"MCMP_"<<".dat";
	ofstream out(name.str().c_str());
	out<<"Title= \"LBM Lid Driven Flow\"\n"<<"VARIABLES=\"X\",\"Y\",\"Rho1\",\"Rho2\",\"Flag\"\n"<<"ZONE T=\"BOX\",I="<<NX+1<<",J="<<NY+1<<",F=POINT"<<endl;
	for(j=0;j<=NY;j++)
	{
		for(i=0;i<=NX;i++)
		{
			out<<i<<" "<<j<<" "<<(double)arho[i][j]<<" "<<(double)brho[i][j]<<" "<<flag[i][j]<<endl;
		}
	}
}

void evolution()
{
    //collision
    for(i=0;i<=NX;i++)
    {
        for(j=0;j<=NY;j++)
        {
            if(flag[i][j]==0)//fluid zone
            {
                for(k=0;k<Q;k++)
                {
                    aF[i][j][k] = af[i][j][k] + (feq(k,arho[i][j],au[i][j]) - af[i][j][k])/atu;
                    bF[i][j][k] = bf[i][j][k] + (feq(k,brho[i][j],bu[i][j]) - bf[i][j][k])/btu;
                }
            }

        }
    }
    //streaming
    for(i=0;i<=NX;i++)
    {
        for(j=0;j<=NY;j++)
        {
            if(flag[i][j]==0)//fluid zone
            {
                for(k=0;k<Q;k++)
                {
                    ip = i - e[k][0];
                    jp = j - e[k][1];
                    if(ip>=0 && ip<=NX && jp>=0 && jp<=NY)
                    {
                        if(flag[ip][jp]==0)//fluid zone
                        {
                            af[i][j][k] = aF[ip][jp][k];
                            bf[i][j][k] = bF[ip][jp][k];
                        }
                        else{
                            af[i][j][k] = aF[i][j][r[k]];
                            bf[i][j][k] = bF[i][j][r[k]];
                        }
                    }

                }
            }
        }
    }
    //boundary condition periodic  af = AF
    for(j=1;j<NY;j++)//left
    {
        af[0][j][5]=aF[NX][j][5];
        af[0][j][1]=aF[NX][j][1];
        af[0][j][8]=aF[NX][j][8];
        bf[0][j][5]=bF[NX][j][5];
        bf[0][j][1]=bF[NX][j][1];
        bf[0][j][8]=bF[NX][j][8];
    }
    for(j=1;j<NY;j++)//right
    {
        af[NX][j][6]=aF[0][j][6];
        af[NX][j][3]=aF[0][j][3];
        af[NX][j][7]=aF[0][j][7];
        bf[NX][j][6]=bF[0][j][6];
        bf[NX][j][3]=bF[0][j][3];
        bf[NX][j][7]=bF[0][j][7];
    }
    for(i=1;i<NX;i++)//TOP
    {
        af[i][0][6]=aF[i][NY][6];
        af[i][0][5]=aF[i][NY][5];
        af[i][0][4]=aF[i][NY][4];
        bf[i][0][6]=bF[i][NY][6];
        bf[i][0][5]=bF[i][NY][5];
        bf[i][0][4]=bF[i][NY][4];
    }
    for(i=1;i<NX;i++)//BOTTOM
    {
        af[i][NY][7]=aF[i][0][7];
        af[i][NY][4]=aF[i][0][4];
        af[i][NY][8]=aF[i][0][8];
        bf[i][NY][7]=bF[i][0][7];
        bf[i][NY][4]=bF[i][0][4];
        bf[i][NY][8]=bF[i][0][8];
    }
    // corner

    // 0,NX
    af[0][NY][4]=aF[0][0][4];
    af[0][NY][7]=aF[0][0][7];
    af[0][NY][5]=aF[NX][NY][5];
    af[0][NY][1]=aF[NX][NY][1];
    af[0][NY][8]=aF[NX][0][8];
    // 0,0
    af[0][0][6]=aF[0][NY][6];
    af[0][0][2]=aF[0][NY][2];
    af[0][0][1]=aF[NX][0][1];
    af[0][0][8]=aF[NX][0][8];
    af[0][0][5]=aF[NX][NY][5];
    // NX,NY
    af[NX][NY][6]=aF[0][NY][6];
    af[NX][NY][3]=aF[0][NY][3];
    af[NX][NY][4]=aF[NX][0][4];
    af[NX][NY][8]=aF[NX][0][8];
    af[NX][NY][7]=aF[0][0][7];
    // NX,0
    af[NX][0][3]=aF[0][0][3];
    af[NX][0][7]=aF[0][0][7];
    af[NX][0][2]=aF[NX][NY][2];
    af[NX][0][5]=aF[NX][NY][5];
    af[NX][0][6]=aF[0][NY][6];
    // component B
    // 0,NX
    bf[0][NY][4]=bF[0][0][4];
    bf[0][NY][7]=bF[0][0][7];
    bf[0][NY][5]=bF[NX][NY][5];
    bf[0][NY][1]=bF[NX][NY][1];
    bf[0][NY][8]=bF[NX][0][8];
    // 0,0
    bf[0][0][6]=bF[0][NY][6];
    bf[0][0][2]=bF[0][NY][2];
    bf[0][0][1]=bF[NX][0][1];
    bf[0][0][8]=bF[NX][0][8];
    bf[0][0][5]=bF[NX][NY][5];
    // NX,NY
    bf[NX][NY][6]=bF[0][NY][6];
    bf[NX][NY][3]=bF[0][NY][3];
    bf[NX][NY][4]=bF[NX][0][4];
    bf[NX][NY][8]=bF[NX][0][8];
    bf[NX][NY][7]=bF[0][0][7];
    // NX,0
    bf[NX][0][3]=bF[0][0][3];
    bf[NX][0][7]=bF[0][0][7];
    bf[NX][0][2]=bF[NX][NY][2];
    bf[NX][0][5]=bF[NX][NY][5];
    bf[NX][0][6]=bF[0][NY][6];
}

void marco()
{
    for(i=0;i<=NX;i++)
    {
        for(j=0;j<=NY;j++)
        {
            au[i][j][0] = 0;
            au[i][j][1] = 0;
            bu[i][j][0] = 0;
            bu[i][j][1] = 0;
            arho[i][j] = 0;
            brho[i][j] = 0;
            for(k=0;k<Q;k++)
            {
                arho[i][j]+=af[i][j][k];
                brho[i][j]+=bf[i][j][k];
                //
                au[i][j][0]+=(af[i][j][k]*e[k][0]);
                au[i][j][1]+=(af[i][j][k]*e[k][1]);
                bu[i][j][0]+=(bf[i][j][k]*e[k][0]);
                bu[i][j][1]+=(bf[i][j][k]*e[k][1]);
            }
            if(arho[i][j]==0)
            {
                 au[i][j][0]=0;
                 au[i][j][1]=0;
            }
            else{
                au[i][j][0]/=(arho[i][j]);
                au[i][j][1]/=(arho[i][j]);
            }
            if(brho[i][j]==0)
            {
                bu[i][j][0]=0;
                bu[i][j][0]=0;
            }
            else{
                bu[i][j][0]/=(brho[i][j]);
                bu[i][j][1]/=(brho[i][j]);
            }

            // macro speed
            if(arho[i][j]==0&&brho[i][j]==0)
            {
                u[i][j][0]=0;
                u[i][j][1]=0;
            }
            else{
                u[i][j][0] = (arho[i][j]*au[i][j][0]/atu + brho[i][j]*bu[i][j][0]/btu) / (arho[i][j]/atu + brho[i][j]/btu );
                u[i][j][1] = (arho[i][j]*au[i][j][1]/atu + brho[i][j]*bu[i][j][1]/btu) / (arho[i][j]/atu + brho[i][j]/btu );
            }
            // each component speed
            if(arho[i][j]==0)
            {
                au[i][j][0]=0;
                au[i][j][1]=0;
            }
            else{
                au[i][j][0] = u[i][j][0] + atu * ( xforceFluid(i,j,1) + xforceSolid(i,j,1) ) / (arho[i][j]);
                au[i][j][1] = u[i][j][1] + atu * ( yforceFluid(i,j,1) + yforceSolid(i,j,1) ) / (arho[i][j]);
            }
            if(brho[i][j]==0)
            {
                bu[i][j][0]=0;
                bu[i][j][1]=0;
            }
            else{
                bu[i][j][0] = u[i][j][0] + btu * (xforceFluid(i,j,2) + xforceSolid(i,j,2)) / (brho[i][j]);
                bu[i][j][1] = u[i][j][1] + btu * (yforceFluid(i,j,2) + yforceSolid(i,j,2)) / (brho[i][j]);
            }
            if(flag[i][j]==1)//solid
            {
                arho[i][j]=0;
                brho[i][j]=0;
                for(k=0;k<Q;k++)
                {
                    af[i][j][k]=0;
                    bf[i][j][k]=0;
                }
            }
        }
    }
}

void error()
{
}


void output(int m)
{
    ostringstream name;
	name<<"MCMP_"<<m<<".dat";
	ofstream out(name.str().c_str());
	out<<"Title= \"LBM Lid Driven Flow\"\n"<<"VARIABLES=\"X\",\"Y\",\"U\",\"V\",\"Rho1\",,\"Rho2\",\"Flag\"\n"<<"ZONE T=\"BOX\",I="<<NX+1<<",J="<<NY+1<<",F=POINT"<<endl;
	for(j=0;j<=NY;j++)
	{
		for(i=0;i<=NX;i++)
		{
			out<<i<<" "<<j<<" "<<au[i][j][0]<<" "<<au[i][j][1]<<" "<<(double)arho[i][j]<<" "<<(double)brho[i][j]<<" "<<flag[i][j]<<endl;
		}
	}
}

double feq(int k,double rho,double u[2])
{
    double eu,uv,feq;
    eu=(e[k][0]*u[0]+e[k][1]*u[1]);
    uv=(u[0]*u[0]+u[1]*u[1]);
    feq=w[k]*rho*(1.0+3.0*eu+4.5*eu*eu-1.5*uv);
    return feq;
}

double xforceFluid(int i, int j, int component)
{
    double xforcefluidnumber;
    switch(component)
    {
        case 1 : xforcefluidnumber = -gff*aphi(i,j)*(1/3.0*(bphi(i+1,j)-bphi(i-1,j))+1/12.0*(bphi(i+1,j+1)-bphi(i-1,j+1))+1/12.0*(bphi(i+1,j-1)-bphi(i-1,j-1)));break;
        case 2 : xforcefluidnumber = -gff*bphi(i,j)*(1/3.0*(aphi(i+1,j)-aphi(i-1,j))+1/12.0*(aphi(i+1,j+1)-aphi(i-1,j+1))+1/12.0*(aphi(i+1,j-1)-aphi(i-1,j-1)));break;
        default: cout<<"xforceFluid is wrong"<<endl;
    }
    return xforcefluidnumber;
}

double yforceFluid(int i, int j, int component)
{
    double yforcefluidnumber;
    switch(component)
    {
        case 1 : yforcefluidnumber = -gff*aphi(i,j)*(1/3.0*(bphi(i,j+1)-bphi(i,j-1))+1/12.0*(bphi(i+1,j+1)-bphi(i+1,j-1))+1/12.0*(bphi(i-1,j+1)-bphi(i-1,j-1)));break;
        case 2 : yforcefluidnumber = -gff*bphi(i,j)*(1/3.0*(aphi(i,j+1)-aphi(i,j-1))+1/12.0*(aphi(i+1,j+1)-aphi(i+1,j-1))+1/12.0*(aphi(i-1,j+1)-aphi(i-1,j-1)));break;
        default: cout<<"yforceFluid is wrong"<<endl;
    }
    return yforcefluidnumber;
}

double xforceSolid(int i, int j, int component)
{
    double xforcesolidnubmebr;
    switch(component)
    {
        case 1 : xforcesolidnubmebr = -gfsa*arho[i][j]*(1.0*(s(i+1,j)-s(i-1,j)) + 1/4.0*(s(i+1,j+1) - s(i-1,j+1) + s(i+1,j-1) - s(i-1,j-1)));break;
        case 2 : xforcesolidnubmebr = -gfsb*brho[i][j]*(1.0*(s(i+1,j)-s(i-1,j)) + 1/4.0*(s(i+1,j+1) - s(i-1,j+1) + s(i+1,j-1) - s(i-1,j-1)));break;
        default: cout<<"xforceSolid is wrong"<<endl;
    }
    return xforcesolidnubmebr;
}

double yforceSolid(int i, int j, int component)
{
    double yforcesolidnumber;
    switch(component)
    {
        case 1 : yforcesolidnumber = -gfsa*arho[i][j]*(1.0*(s(i,j+1)-s(i,j-1)) + 1/4.0*(s(i+1,j+1) + s(i-1,j+1) - s(i-1,j-1) - s(i+1,j-1)));break;
        case 2 : yforcesolidnumber = -gfsb*brho[i][j]*(1.0*(s(i,j+1)-s(i,j-1)) + 1/4.0*(s(i+1,j+1) + s(i-1,j+1) - s(i-1,j-1) - s(i+1,j-1)));break;
        default: cout<<"yforceSolid is wrong"<<endl;
    }
    return yforcesolidnumber;
}

double aphi(int i, int j)
{
    double aphinumber;
    if(i<0||i>NX||j<0||j>NY)
    {
        aphinumber = 0;
    }
    else
        aphinumber = 1.0 - exp(-arho[i][j]);
    return aphinumber;
}

double bphi(int i, int j)
{
    double bphinumber;
    if(i<0||i>NX||j<0||j>NY)
    {
        bphinumber = 0;
    }
    else
        bphinumber = 1.0 - exp(brho[i][j]);
    return bphinumber;
}
int s(int i, int j)
{
    int snumber;
    if(i<0||i>NX||j<0||j>NY)
    {
        snumber = 0;
    }
    else
        snumber = flag[i][j];
    return snumber;
}
