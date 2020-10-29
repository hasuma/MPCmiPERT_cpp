//SIMULACAO PARA CONTROLE DE ATITUDE(THETA4)
//MPC com malha interna(LQR) e ESTIMATIVA de PERTUBACAO 
//**********************************************************
// Bibliotecas
#include <bits/stdc++.h>//All basic librarys
#include <armadillo>	//std::arma
#include <boost/numeric/odeint.hpp>//std::integrate
#include "include/Controlador/CONTROL.hpp"
#include "include/libBANCADA.h"
//**********************************************************		
using namespace arma;
using namespace std;
using namespace boost::numeric::odeint;
using namespace libBANCADA;
//**********************************************************
//Variaveis Globais
#define pi 3.14159265
#define rad2deg (180/3.14159265358979323846)
#define deg2rad (3.14159265358979323846/180)

int kmax;
double Ts,Td;
double dt_BOOST = Td/40;
string file2save;

int numStates = 4;
int numInput = 2;
int numRef = 2;
int numStatesEX = numStates + numInput;
double up0[2];
mat posIni(2,1);
mat yref(2,1);

//Coeficientes F - PWM MOTOR1
#define aF1 0.000458547428571
#define bF1 0.036316433142860
#define cF1 -0.944062080000104
//Coeficientes F - PWM MOTOR2
#define aF2 0.000796758857143
#define bF2 0.009037906285718
#define cF2 -0.285379440000110


//**********************************************************		
//Dinamica da atitude do bicoptero
class fnlin{
	mat u;
	//Parametros da bancada
	double g = 9.81;
    	double l2 = 0.173;
    	double Ix4 = 0.01;
    	double Iy2 = 0.306954;
    	double Iz1 = 0.584807507;
    	double b1 = 0.217;
   	double b2 = 0.215;

    	//Parametros de theta4
    	double mi_d4 = 0.0959 * 1.0;
    	double CG = 0.6376 * 1.0;
    	//Parametros de theta2
    	double mi_d2 = 1.598815;
    	//Parametros do contra-peso
    	double mcp = 1.709;
    	double lcp = 0.33;
    	//Condiçoes de contorno/iniciais
    	mat u0 = {1.65,1.7};
	//Parametros do motor
    	mat p = {{8.016, 0},
                 {0, 9.072}};
     
    	mat coef1 = {-0.000063260930216, 0.013139466306494, -0.175069018957576};
    	mat coef2 = {-0.000021708654545, 0.010473514036364, -0.077660122909090};
	mat q = zeros<mat>(2,2);
	//Parametros do tilt
    	mat pt = {{7.8413, 0},
               	  {0, 7.5815}};
    	mat qt = eye<mat>(2,2);
public:
	fnlin(mat param) : u(param) {}
	void operator()(const vector<double> &x, vector<double> &dxdt, const double t){
		//**********************************************************		
		//Dinamica dos motores
		q(0,0) = coef1(0)*pow(u(0),2) + coef1(1)*u(0) + coef1(2);
		q(1,1) = coef2(0)*pow(u(1),2) + coef2(1)*u(1) + coef2(2);
		mat uF = {u(0),u(1)};
		uF = uF.t();
		mat Fm = {x[2],x[3]};
		Fm = Fm.t();
		mat Fmp = -p * Fm + q * uF;

		//**********************************************************
		//Dinamica dos tilts
		mat alfa = {0,0};
		//alfa = alfa.t();
		//mat alfap = -pt * alfa + qt * u.rows(2,3);

		//**********************************************************
		//Dinamica Bancada
		double x3 = -x[0]-pi*0.5;
		mat q = {0,x[0],x3,x[1]};
		mat qdot = {0,x[4],-x[4],x[5]};
	
		mat Dq,Cq,Gq,nJp;
        	Dq = dinamica::Dq(q(0),q(1),q(2),q(3),qdot(0),qdot(1),qdot(2),qdot(3));
        	Cq = dinamica::Cq(q(0),q(1),q(2),q(3),qdot(0),qdot(1),qdot(2),qdot(3));
        	Gq = dinamica::Gq(q(0),q(1),q(2),q(3),qdot(0),qdot(1),qdot(2),qdot(3));
        	nJp = dinamica::nJp(q(0),q(1),q(2),q(3),qdot(0),qdot(1),qdot(2),qdot(3));	
		
		//Definindo termo Gamma b
		mat tilt = {{  -cos(alfa(0)),   -cos(alfa(1))},
             		    { 	sin(alfa(0)),   -sin(alfa(1))},
                       	    {		   0,               0},
                            {		   0,               0},
                            {		   0,               0},
          	            {b1*cos(alfa(0)),-b2*cos(alfa(1))}};
		mat gamma_b = nJp * tilt * Fm;

		mat Fdot = {0,lcp*mcp*9.81*cos(q(1)) - mi_d2*qdot(1),0,-mi_d4*qdot(3) - CG*q(3)};
		Fdot = Fdot.t();
		
		mat qddot = inv(Dq) * (-Cq + Fdot - Gq + gamma_b); 
		
		//**********************************************************
		//Explicit Function
		dxdt[0] = x[4];
		dxdt[1] = x[5];
		dxdt[2] = Fmp(0);
		dxdt[3] = Fmp(1);
		dxdt[4] = qddot(1);
		dxdt[5] = qddot(3);
	}
};

//**********************************************************
// Rotina principal		
int main(){
	ifstream parametros_sistema;
	parametros_sistema.open("parametros_sistema.txt",ifstream::in);

	if (parametros_sistema.is_open())
	{
		parametros_sistema>> Ts;         	   //time sample
		parametros_sistema>> Td;		
		parametros_sistema>> kmax;      	   //number of iterations
		parametros_sistema>> up0[0];       // PWM value at hoover condition
		parametros_sistema>> up0[1];       // PWM value at hoover condition
		parametros_sistema>> file2save;            // file to save data name
		parametros_sistema>> posIni(0,0);
		parametros_sistema>> posIni(1,0);
		parametros_sistema>> yref(0,0);		
		parametros_sistema>> yref(1,0);             // Theta4 reference value
	}
	else{
		cerr << "Arquivo nao existente ou nao encontrado! "<< endl;
		return 1;
	}
	parametros_sistema.close();
	
	cout << "Parametros Iniciais da SIMULACAO"<< endl;
	cout << "Nome do ensaio: "<< file2save <<  endl;
	cout << "Periodo de amostragem: "<< Ts <<  endl;
	cout << "Atraso de transporte: "<< Td <<  endl;
	cout << "Posicao inicial: "<< posIni <<  endl;
	cout << "Referencia: "<< yref<<  endl;
	cout << "Forcas de trimagem "<< up0[0] << "\t" << up0[1]<<  endl;
	yref = yref * deg2rad;
	posIni = posIni * deg2rad;
	//**********************************************************
	//Inicializacao do controle
	cout <<"INICIALIZACAO DO CONTROLE"<<endl;
	Control controle(numStates, numInput, numRef, posIni,yref, Ts, up0);
	
	//**********************************************************
	//Variaveis de estado e controle
	// Estado x
	mat x = zeros<mat>(numStatesEX,kmax+1);
	// Entrada
	mat up = zeros<mat>(numInput,kmax);
	mat pwm = zeros<mat>(numInput,kmax);
	// Saida
	mat y = zeros<mat>(numRef,kmax);
	mat ym = zeros<mat>(numRef,kmax);
	//CPU time
	mat cpu_time = zeros<mat>(kmax,1);
	clock_t tictoc;

	//**********************************************************
	//Condicoes iniciais
	vector<double> xini(numStatesEX);
	x(0,0) = posIni(1);
	x(1,0) = posIni(1);
	x(2,0) = up0[0];
	x(3,0) = up0[1];
	mat pwmkm1 = zeros<mat>(numInput,1);
	mat xAux(numStates,1);
	mat u0 = {up0[0],up0[1]};
	u0 = u0.t();
	//**********************************************************
	//Ruido posicao angular - Range de 0.05 rad
	//Ruido velocidade angular - Range de 0.4 rad/s
	mat ruido(4,kmax);
	for(int i = 0; i<kmax;i++){
		ruido(0,i) = ((double)((rand()%5001)-2500))/((double)100000);
		ruido(1,i) = ((double)((rand()%5001)-2500))/((double)100000);
		ruido(2,i) = ((double)((rand()%28001)-14000))/((double)100000);
		ruido(3,i) = ((double)((rand()%28001)-14000))/((double)100000);
	}
	
	//**********************************************************
	// Evolucao da dinamica da planta
	for(int k = 0; k<kmax;k++){
		//**********************************************************
		tictoc = clock();

		//**********************************************************
		xAux = {x(0,k),
			x(1,k),
			x(4,k),
		     	x(5,k)};
		xAux = xAux.t();
		y.col(k) = xAux.rows(0,1);
		xAux = xAux + ruido.col(k);
		ym.col(k) = xAux.rows(0,1);
		//**********************************************************
		//up.col(k) = u0;
		up.col(k) = controle.INPUT(xAux,ym.col(k));
		pwm(0,k) = ((-bF1 + sqrt(bF1*bF1-4*aF1*(cF1 - up(0,k))))/(2*aF1));
		pwm(1,k) = ((-bF2 + sqrt(bF2*bF2-4*aF2*(cF2 - up(1,k))))/(2*aF2));
		//**********************************************************
		tictoc = clock()-tictoc;		
		cpu_time(k) = ((float)tictoc)/CLOCKS_PER_SEC*1000000;
		//**********************************************************
		
		//**********************************************************
		// Simulacao da dinamica da planta
		for(int j = 0;j < numStatesEX;j++){xini[j] = x(j,k);}
		// 0 - Td
		
		fnlin fnlinTd(pwmkm1);
		integrate(fnlinTd,xini,0.0,Td,0.00001);
		//Td - Ts
		fnlin fnlinTs(pwm.col(k));
		integrate(fnlinTs,xini,0.0,Ts-Td,0.00001);
		
		for(int z = 0;z < numStatesEX;z++){x(z,k+1) = xini[z];}

		//**********************************************************
		pwmkm1 = pwm.col(k);
	}

	//**********************************************************
	//Armazenamento dos dados
	y = y.t() * 180/pi;
	ym = ym.t() * 180/pi;
	up = up.t();
	pwm = pwm.t();
	//du = du.t();
	x = x.cols(0,kmax-1).t();

	mat k = linspace<mat>(0,kmax-1,kmax);
	mat data;

	data = join_rows(k,y);
	data = join_rows(data,ym);
	data = join_rows(data,pwm);
	data = join_rows(data,x);
	data = join_rows(data,up);
	//data = join_rows(data,du);
	data = join_rows(data,cpu_time);
	data.save("Output/data.txt", raw_ascii);

	//**********************************************************
	//Definindo cpu time medio para calcular a acao de controle
	mat time = mean(cpu_time,0);
	mat dtime = stddev(cpu_time,0);
	cout<<"CPU time: "<<time(0)<<"+/-"<<2*dtime(0)<<" microsegundos, para uma confianca de 95%."<<endl;

	return 0;
}
