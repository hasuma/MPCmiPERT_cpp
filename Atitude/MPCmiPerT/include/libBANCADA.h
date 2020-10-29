#pragma once
#include <armadillo>
using namespace arma;

namespace libBANCADA{
	double g = 9.81;
  	mat m = {0,2.761, 0.154, 0.115, 0.559};
  	mat l = {0,    0, 0.595, 0.173,     0};
  	mat Ix= {0,0.004165125, 0.001006095, 0.00001723333, 0.01};
  	mat Iy= {0,0.289994709, 0.306954, 0.00008259913, 0.00051408};
  	mat Iz= {0,0.28761702, 0.00001835292, 0.00007608364, 0.01329993};
	class dinamica{
		public:
			static mat Dq(double q1,double q2,double q3,double q4,double qdot1,double qdot2,double qdot3,double qdot4);
			static mat Cq(double q1,double q2,double q3,double q4,double qdot1,double qdot2,double qdot3,double qdot4);
			static mat Gq(double q1,double q2,double q3,double q4,double qdot1,double qdot2,double qdot3,double qdot4);
			static mat nJp(double q1,double q2,double q3,double q4,double qdot1,double qdot2,double qdot3,double qdot4);
	};
}
