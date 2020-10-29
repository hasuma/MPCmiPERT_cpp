#include "libBANCADA.h"
#include <armadillo>

using namespace arma;

namespace libBANCADA{
	mat dinamica::Dq(double q1,double q2,double q3,double q4,double qdot1,double qdot2,double qdot3,double qdot4){
		mat A = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
		A(0,0) = (-pow(l(3),2) * (m(2) + 4 * m(3) + 4 * m(4)) * pow(sin(q3),2) + pow(l(3),2) * (m(2) + 4 * m(3) + 4 * m(4)) * pow(cos(q3),2) + 0.4e1 * l(2) * l(3) * (m(2) + 2 * m(3) + 2 * m(4)) * cos(q3) + 0.4e1 * pow(l(2),2) * (m(2) + m(3) + m(4))) * pow(cos(q2),2) / 0.4e1 - sin(q2) * l(3) * sin(q3) * (l(3) * (m(2) + 4 * m(3) + 4 * m(4)) * cos(q3) + 0.2e1 * l(2) * (m(2) + 2 * m(3) + 2 * m(4))) * cos(q2) / 0.2e1 + pow(l(3),2) * (m(2) + 4 * m(3) + 4 * m(4)) * pow(sin(q3),2) / 0.4e1 + Iz(4) + Iz(1) + Iz(2) + Iz(3);
		A(0,1) = 0.0e0;
		A(0,2) = 0.0e0;
		A(0,3) = (-sin(q2) * sin(q3) + cos(q2) * cos(q3)) * Iz(4);
		A(1,0) = 0.0e0;
		A(1,1) = (-4 * Ix(2) - 4 * Ix(3) - 4 * Ix(4) + 4 * Iy(2) + 4 * Iy(3) + 4 * Iy(4)) * pow(cos(q1),2) / 0.4e1 + pow(l(3),2) * (m(2) + 4 * m(3) + 4 * m(4)) * pow(cos(q3),2) / 0.4e1 + l(2) * l(3) * (m(2) + 2 * m(3) + 2 * m(4)) * cos(q3) + pow(l(3),2) * (m(2) + 4 * m(3) + 4 * m(4)) * pow(sin(q3),2) / 0.4e1 + (4 * m(2) + 4 * m(3) + 4 * m(4)) * pow(l(2),2) / 0.4e1 + Ix(2) + Ix(3) + Ix(4);
		A(1,2) = (-Ix(3) - Ix(4) + Iy(3) + Iy(4)) * pow(cos(q1),2) + pow(l(3),2) * (m(3) + m(4)) * pow(cos(q3),2) + l(2) * l(3) * (m(3) + m(4)) * cos(q3) + pow(l(3),2) * (m(3) + m(4)) * pow(sin(q3),2) + Ix(3) + Ix(4);
		A(1,3) = -cos(q1) * sin(q1) * (Ix(4) - Iy(4)) * (sin(q2) * cos(q3) + cos(q2) * sin(q3));
		A(2,0) = 0.0e0;
		A(2,1) = A(1,2);
		A(2,2) = (-Ix(3) - Ix(4) + Iy(3) + Iy(4)) * pow(cos(q1),2) + pow(l(3),2) * (m(3) + m(4)) * pow(cos(q3),2) + pow(l(3),2) * (m(3) + m(4)) * pow(sin(q3),2) + Ix(3) + Ix(4);
		A(2,3) = -cos(q1) * sin(q1) * (Ix(4) - Iy(4)) * (sin(q2) * cos(q3) + cos(q2) * sin(q3));
		A(3,0) = A(0,3);
		A(3,1) = A(1,3);
		A(3,2) = A(2,3);
		A(3,3) = -((Ix(4) - Iy(4)) * pow(cos(q1),2) + Iy(4) - Iz(4)) * (cos(q3) - sin(q3)) * (cos(q3) + sin(q3)) * pow(cos(q2),2) + 0.2e1 * ((Ix(4) - Iy(4)) * pow(cos(q1),2) + Iy(4) - Iz(4)) * sin(q2) * sin(q3) * cos(q3) * cos(q2) + ((Ix(4) - Iy(4)) * pow(cos(q1),2) + Iy(4)) * pow(cos(q3),2) + pow(sin(q3),2) * Iz(4);		
		return A;
	}
	mat dinamica::Cq(double q1,double q2,double q3,double q4,double qdot1,double qdot2,double qdot3,double qdot4){
		mat A = {0,0,0,0};
  		A(0) = -0.2e1 * qdot2 * pow(l(2),2) * (m(2) + m(3) + m(4)) * cos(q2) * sin(q2) * qdot1 - (-qdot2 * (-4 * Ix(2) - 4 * Ix(3) - 4 * Ix(4) + 4 * Iy(2) + 4 * Iy(3) + 4 * Iy(4)) * cos(q1) * sin(q1) / 0.2e1 - 0.2e1 * qdot3 * (-Ix(3) - Ix(4) + Iy(3) + Iy(4)) * cos(q1) * sin(q1) + qdot4 * (-(Ix(4) - Iy(4)) * pow(sin(q1),2) + (Ix(4) - Iy(4)) * pow(cos(q1),2))) * qdot2 / 0.2e1 - (-0.2e1 * qdot2 * (-Ix(3) - Ix(4) + Iy(3) + Iy(4)) * cos(q1) * sin(q1) - 0.2e1 * qdot3 * (-Ix(3) - Ix(4) + Iy(3) + Iy(4)) * cos(q1) * sin(q1) + qdot4 * (-(Ix(4) - Iy(4)) * pow(sin(q1),2) + (Ix(4) - Iy(4)) * pow(cos(q1),2))) * qdot3 / 0.2e1 - (qdot2 * (-(Ix(4) - Iy(4)) * pow(sin(q1),2) + (Ix(4) - Iy(4)) * pow(cos(q1),2)) + qdot3 * (-(Ix(4) - Iy(4)) * pow(sin(q1),2) + (Ix(4) - Iy(4)) * pow(cos(q1),2)) - 0.2e1 * qdot4 * cos(q1) * sin(q1) * (Ix(4) - Iy(4))) * qdot4 / 0.2e1;
  		A(1) = (-qdot1 * (-4 * Ix(2) - 4 * Ix(3) - 4 * Ix(4) + 4 * Iy(2) + 4 * Iy(3) + 4 * Iy(4)) * cos(q1) * sin(q1) / 0.2e1 - qdot2 * l(2) * l(3) * (m(2) + 2 * m(3) + 2 * m(4)) * cos(q2)) * qdot2 + (-0.2e1 * qdot1 * (-Ix(3) - Ix(4) + Iy(3) + Iy(4)) * cos(q1) * sin(q1) - qdot2 * l(2) * l(3) * (m(3) + m(4)) * cos(q2)) * qdot3 + qdot1 * (-(Ix(4) - Iy(4)) * pow(sin(q1),2) + (Ix(4) - Iy(4)) * pow(cos(q1),2)) * qdot4 + pow(qdot1,2) * pow(l(2),2) * (m(2) + m(3) + m(4)) * cos(q2) * sin(q2) - (-qdot2 * l(2) * l(3) * (m(2) + 2 * m(3) + 2 * m(4)) * cos(q2) - qdot3 * l(2) * l(3) * (m(3) + m(4)) * cos(q2)) * qdot2 / 0.2e1 + qdot2 * l(2) * l(3) * (m(3) + m(4)) * cos(q2) * qdot3 / 0.2e1;
  		A(2) = (-0.2e1 * qdot1 * (-Ix(3) - Ix(4) + Iy(3) + Iy(4)) * cos(q1) * sin(q1) - qdot2 * l(2) * l(3) * (m(3) + m(4)) * cos(q2)) * qdot2 - 0.2e1 * qdot1 * (-Ix(3) - Ix(4) + Iy(3) + Iy(4)) * cos(q1) * sin(q1) * qdot3 + qdot1 * (-(Ix(4) - Iy(4)) * pow(sin(q1),2) + (Ix(4) - Iy(4)) * pow(cos(q1),2)) * qdot4;
  		A(3) = qdot1 * (-(Ix(4) - Iy(4)) * pow(sin(q1),2) + (Ix(4) - Iy(4)) * pow(cos(q1),2)) * qdot2 + qdot1 * (-(Ix(4) - Iy(4)) * pow(sin(q1),2) + (Ix(4) - Iy(4)) * pow(cos(q1),2)) * qdot3 - 0.2e1 * qdot1 * cos(q1) * sin(q1) * (Ix(4) - Iy(4)) * qdot4;
  		return A.t();
	}
	mat dinamica::Gq(double q1,double q2,double q3,double q4,double qdot1,double qdot2,double qdot3,double qdot4){
		mat A = {0,0,0,0};
	        A(0) = 0;
	        A(1) = 0*m(2) * g * l(2) * cos(q2) +0* m(3) * g * l(2) * cos(q2) + m(4) * g * l(2) * cos(q2);
	        A(2) = 0;
	        A(3) = 0;
	        return A.t();
	}
	mat dinamica::nJp(double q1,double q2,double q3,double q4,double qdot1,double qdot2,double qdot3,double qdot4){
		mat A = {{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0}};
	  	A(0,0) = sin(q4) * l(2) * cos(q2);
	  	A(0,1) = cos(q4) * l(2) * cos(q2);
	  	A(0,2) = 0.0e0;
	  	A(0,3) = -cos(q4);
	  	A(0,4) = sin(q4);
	  	A(0,5) = 0.0e0;
	  	A(1,0) = -A(0,1);
	  	A(1,1) = sin(q4) * l(2) * cos(q2);
	  	A(1,2) = -l(2) * sin(q2) + l(3);
	  	A(1,3) = -sin(q4);
	  	A(1,4) = -cos(q4);
	  	A(1,5) = 0.0e0;
	  	A(2,0) = 0.0e0;
	  	A(2,1) = 0.0e0;
	  	A(2,2) = l(3);
	  	A(2,3) = -sin(q4);
	  	A(2,4) = -cos(q4);
	  	A(2,5) = 0.0e0;
	  	A(3,0) = 0.0e0;
	  	A(3,1) = 0.0e0;
	  	A(3,2) = 0.0e0;
	  	A(3,3) = 0.0e0;
	  	A(3,4) = 0.0e0;
	  	A(3,5) = 0.1e1;
	   	return A;
	}
	
}
