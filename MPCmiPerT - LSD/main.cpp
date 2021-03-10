//SIMULACAO PARA CONTROLE DE ATITUDE(THETA4)
//MPC com malha interna(LQR) e ESTIMATIVA de PERTUBACAO 
//**********************************************************
// Bibliotecas
#include <bits/stdc++.h>//All basic librarys
#include <armadillo>	//std::arma
#include <boost/numeric/odeint.hpp>//std::integrate
#include "include/Controlador/CONTROL.hpp"
//**********************************************************		
using namespace arma;
using namespace std;
using namespace boost::numeric::odeint;
//**********************************************************
//Variaveis Globais
#define pi 3.14159265
#define rad2deg (180/3.14159265358979323846)
#define deg2rad (3.14159265358979323846/180)

int kmax;
double Ts,Td;
double dt_BOOST = Td/40;
string file2save;

int numStates = 6;
int numInput = 4;
int numRef = 3;
int numStatesEX = numStates + numInput;
double up0[2];
mat posIni(3,1);
mat yref(3,1);
double ref1;
double ref2;
double ref4;
//Coeficientes F - PWM MOTOR1
#define aF1 0.000386794285714
#define bF1 0.044941011428573
#define cF1 -1.149182640000060
//Coeficientes F - PWM MOTOR2
#define aF2 0.000617749714286
#define bF2 0.020354348571431
#define cF2 -0.502036560000080

//Coefientes T - PWM TILT1
double pt1 = 7.8413;
//Coefientes T - PWM TILT2
double pt2 = 7.5815;
//**********************************************************		
//Dinamica da atitude do bicoptero
class fnlin{
	mat u;
	//Parametros da bancada
	double g = 9.81;
    	double l2 = 0.173;
    	double Ix4 = 0.01 * 1;
    	double Iy2 = 0.306954 * 1;
    	double Iz1 = 0.584807507 * 1;
    	double b1 = 0.215;
   	double b2 = 0.215;

    	//Parametros de theta4
    	double mi_d4 = 0.0959 * 1;
    	double CG = 0.6376 * 1;
    	//Parametros de theta2
    	double mi_d2 = 1.598815 * 1;
    	//Parametros do contra-peso
    	double mcp = 1.709;
    	double lcp = 0.33;
    	//Condiçoes de contorno/iniciais
    	mat u0 = {1.7,1.7};
	//Parametros do motor
    	mat p = {{8.016, 0},
                 {0, 9.072}};
     
    	mat coef1 = {-0.000063260930216, 0.013139466306494, -0.175069018957576};
    	mat coef2 = {-0.000021708654545, 0.010473514036364, -0.077660122909090};
	mat qF = zeros<mat>(2,2);
	//Parametros do tilt
    	mat pt = {{7.8413, 0},
               	  {0, 7.5815}};
    	mat qt = eye<mat>(2,2);
public:
	fnlin(mat param) : u(param) {}
	void operator()(const vector<double> &x, vector<double> &dxdt, const double t){
		//**********************************************************		
		//Dinamica dos motores
		qF(0,0) = coef1(0)*pow(u(0),2) + coef1(1)*u(0) + coef1(2);
		qF(1,1) = coef2(0)*pow(u(1),2) + coef2(1)*u(1) + coef2(2);
		//mat uF = {u(0),u(1)};
		mat uF = {0,0};		
		uF = uF.t();
		mat Fm = {x[3],x[4]};
		Fm = Fm.t();
		mat Fmp = -p * Fm + qF * uF;
		//**********************************************************
		//Dinamica dos tilts
		//mat uA = {u(2),u(3)};
		mat uA = {0,0};		
		uA = uA.t();
		mat alfa = {x[5],x[6]};
		alfa = alfa.t();
		mat alfap = -pt * alfa + qt * uA;

		//**********************************************************
		//Dinamica Bancada
		double x3 = -x[1]-pi*0.5;
		mat q = {x[0],x[1],x3,x[2]};
		mat qdot = {x[7],x[8],-x[8],x[9]};
	
		//double qddot1 = (1/Iz1) * l2 * (u0(0)*alfa(0) - u0(1)*alfa(1) - q(3)*(u0(0)+u0(1)));
    		//double qddot2 = (1/Iy2) * (l2*(Fm(0) - u0(0)) + l2*(Fm(1) - u0(1)) - mi_d2*qdot(1));
    		//double qddot4 = (1/Ix4) * (b1*(Fm(0) - u0(0)) - b2*(Fm(1) - u0(1)) - mi_d4*qdot(3) - CG*q(3));
    		
		double qddot1 = (1/Iz1) * l2 * (u0(0)*u(2) - u0(1)*u(3) - q(3)*(u0(0)+u0(1)));
    		double qddot2 = (1/Iy2) * (l2*(u(0) - u0(0)) + l2*(u(1) - u0(1)) - mi_d2*qdot(1));
    		double qddot4 = (1/Ix4) * (b1*(u(0) - u0(0)) - b2*(u(1) - u0(1)) - mi_d4*qdot(3) - CG*q(3));
		//**********************************************************
		//Explicit Function
		dxdt[0] = x[7];
		dxdt[1] = x[8];
		dxdt[2] = x[9];
		dxdt[3] = Fmp(0);
		dxdt[4] = Fmp(1);
		dxdt[5] = alfap(0);
		dxdt[6] = alfap(1);
		dxdt[7] = qddot1;
		dxdt[8] = qddot2;
		dxdt[9] = qddot4;
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
		parametros_sistema>> posIni(2,0);
		parametros_sistema>> ref1;
		parametros_sistema>> ref2;	
		parametros_sistema>> ref4;             // Theta4 reference value
	}
	else{
		cerr << "Arquivo nao existente ou nao encontrado! "<< endl;
		return 1;
	}
	parametros_sistema.close();
	yref(0,0) = 0;
	yref(1,0) = 0;
	yref(2,0) = 0;

	cout << "Parametros Iniciais da SIMULACAO"<< endl;
	cout << "Nome do ensaio: "<< file2save <<  endl;
	cout << "Periodo de amostragem: "<< Ts <<  endl;
	cout << "Atraso de transporte: "<< Td <<  endl;
	cout << "Posicao inicial: "<< endl << posIni <<  endl;
	cout << "Forcas de trimagem "<< up0[0] << "\t" << up0[1]<<  endl;
	cout << "k = 0"<< endl;
	cout << "Referencia Theta 1 [graus] "<< yref(0,0)<<  endl;
	cout << "Referencia Theta 2 [graus] "<< yref(1,0)<<  endl;
	cout << "Referencia Theta 3 [graus] "<< yref(2,0) <<  endl;
	cout << "============="<< endl;
	cout << "k = 100"<< endl;
	cout << "Referencia Theta 1 [graus] "<< yref(0,0)<<  endl;
	cout << "Referencia Theta 2 [graus] "<< ref2 <<  endl;
	cout << "Referencia Theta 3 [graus] "<< yref(2,0) <<  endl;
	cout << "============="<< endl;
	cout << "k = 700"<< endl;
	cout << "Referencia Theta 1 [graus] "<< ref1 <<  endl;
	cout << "Referencia Theta 2 [graus] "<< ref2 <<  endl;
	cout << "Referencia Theta 4 [graus] "<< yref(2,0) <<  endl;
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
	mat u = zeros<mat>(numInput,kmax);
	mat up = zeros<mat>(numInput,kmax);
	mat ulqr = zeros<mat>(numInput,kmax);
	mat umpc = zeros<mat>(numInput,kmax);
	mat pwm = zeros<mat>(numInput,kmax);
	mat du = zeros<mat>(numInput,kmax);
	// Saida
	mat y = zeros<mat>(numRef,kmax);
	mat ym = zeros<mat>(numRef,kmax);
	//CPU time
	mat cpu_time = zeros<mat>(kmax,1);
	clock_t tictoc;

	//**********************************************************
	//Condicoes iniciais
	vector<double> xini(numStatesEX);
	x(0,0) = posIni(0);
	x(1,0) = posIni(1);
	x(2,0) = posIni(2);
	x(3,0) = up0[0];
	x(4,0) = up0[1];
	mat uAux = zeros<mat>(numInput,2);
	mat ukm1 = zeros<mat>(numInput,1);
	mat pwmkm1 = zeros<mat>(numInput,1);
	mat xAux(numStates,kmax);
	mat xAux2(numStates,1);
	mat u0 = {up0[0],up0[1],0,0};
	u0 = u0.t();
	//**********************************************************
	//Ruido posicao angular - Range de 0.05 rad
	//Ruido velocidade angular - Range de 0.4 rad/s
	mat ruido(6,kmax);
	for(int i = 0; i<kmax;i++){
		ruido(0,i) = ((double)((rand()%501)-250))/((double)100000);
		ruido(1,i) = ((double)((rand()%501)-250))/((double)100000);
		ruido(2,i) = ((double)((rand()%501)-250))/((double)100000);
		ruido(3,i) = ((double)((rand()%2801)-1400))/((double)100000);
		ruido(4,i) = ((double)((rand()%2801)-1400))/((double)100000);
		ruido(5,i) = ((double)((rand()%2801)-1400))/((double)100000);
	}
	
	//**********************************************************
	// Evolucao da dinamica da planta
	for(int k = 0; k<kmax;k++){
		//**********************************************************
		tictoc = clock();
		//**********************************************************
		if(k==0){
			//controle.setRefTheta1inGraus(ref1);
			controle.setRefTheta2inGraus(ref2);
			//controle.setRefTheta4inGraus(ref4);
		}
		if(k==200){
			//controle.setRefTheta1inGraus(ref1);
			//controle.setRefTheta2inGraus(ref2);
			controle.setRefTheta4inGraus(ref4);
		}

		if(k==400){
			controle.setRefTheta1inGraus(ref1);
			//controle.setRefTheta2inGraus(ref2);
			//controle.setRefTheta4inGraus(ref4);
		}
		//**********************************************************
		xAux2 = {x(0,k),
			x(1,k),
			x(2,k),
			x(7,k),
			x(8,k),
		     	x(9,k)};
		xAux.col(k) = xAux2.t();
		xAux.col(k) = xAux.col(k);// + ruido.col(k);
		ym.col(k) = xAux.submat(0,k,2,k);
		//**********************************************************
		//up.col(k) = u0;
		uAux = controle.INPUT(xAux.col(k),ym.col(k));
		ulqr.col(k) = uAux.col(0);
		umpc.col(k) = uAux.col(1);
		u.col(k) = ulqr.col(k) + umpc.col(k);
		du.col(k) = u.col(k) - ukm1;
		up.col(k) = u.col(k) + u0;		
		
		pwm(0,k) = ((-bF1 + sqrt(bF1*bF1-4*aF1*(cF1 - up(0,k))))/(2*aF1));
		pwm(1,k) = ((-bF2 + sqrt(bF2*bF2-4*aF2*(cF2 - up(1,k))))/(2*aF2));
		pwm(2,k) = pt1 * up(2,k);
		pwm(3,k) = pt2 * up(3,k);
		//**********************************************************
		tictoc = clock()-tictoc;		
		cpu_time(k) = ((float)tictoc)/CLOCKS_PER_SEC*1000000;
		//**********************************************************
		
		//**********************************************************
		// Simulacao da dinamica da planta
		for(int j = 0;j < numStatesEX;j++){xini[j] = x(j,k);}
		// 0 - Td

		//fnlin fnlinTd(pwmkm1);
		fnlin fnlinTd(ukm1);
		integrate(fnlinTd,xini,0.0,Td,0.00001);
		//Td - Ts
		//fnlin fnlinTs(pwm.col(k));
		fnlin fnlinTs(up.col(k));
		integrate(fnlinTs,xini,0.0,Ts-Td,0.00001);
		
		for(int z = 0;z < numStatesEX;z++){x(z,k+1) = xini[z];}
		//**********************************************************
		ukm1 = u.col(k);
		pwmkm1 = pwm.col(k);
	}

	//**********************************************************
	//Armazenamento dos dados
	ulqr = ulqr.t();	
	umpc = umpc.t();
	u = u.t();
	up = up.t();
	up.col(2) = up.col(2) * 180/pi;
	up.col(3) = up.col(3) * 180/pi;
	du = du.t();
	xAux = xAux.t() * 180/pi;

	mat k = linspace<mat>(0,kmax-1,kmax);
	mat data;

	data = join_rows(k,xAux);
	data = join_rows(data,up);
	data = join_rows(data,u);
	data = join_rows(data,ulqr);
	data = join_rows(data,umpc);
	data = join_rows(data,du);	
	data = join_rows(data,cpu_time);
	data.save("Output/data.txt", raw_ascii);

	//**********************************************************
	//Definindo cpu time medio para calcular a acao de controle
	mat time = mean(cpu_time,0);
	mat dtime = stddev(cpu_time,0);
	cout<<"CPU time: "<<time(0)<<"+/-"<<2*dtime(0)<<" microsegundos, para uma confianca de 95%."<<endl;

	return 0;
}