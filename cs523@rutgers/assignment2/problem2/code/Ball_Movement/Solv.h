#pragma once
#include"RigidBody.h"
#include<iostream>
#define THRESHOLD 1e-3
class Solv
{
	struct Contact
	{
		RigidBody	*a, *b;
		
		Eigen::Vector3d p, n, ea, eb;

		bool vf;
	};
	
public:
	Solv();
	~Solv();
	void State_to_Array(RigidBody *rb, double *y);
	void Array_to_State(RigidBody *rb, double *y);
	void Array_to_Bodies(RigidBody bodies[],double y[]);
	void Bodies_to_Array(RigidBody bodies[], double y[]);
	
	void Compute_Force_and_Torque(double t, RigidBody *rb);
	void RunSimulation(RigidBody bodies[], double ori_vertex_data_[], double vertex_date_[]);
	void ODE(double y0_[], double yfinal_[], int len_,
		double t_, double t1, RigidBody bodies[]); \

	void PartialUpdate(double y0_[], double yfinal_[], RigidBody bodies_[], double dt_, int i_);
	void VertexPosition(RigidBody bodies[], double ori_vertex_data_[], double vertex_data_[], unsigned int n_);
	Matrix3d Star(Vector3d a_);



	double* TouchGround(double x_[], double p_[], double l_[], double R[], double ta_, double tb_, double epsilon_, RigidBody bodies_[]);
	bool	CollisionCheck(double x_[], double p_[], double l_[], double R[], double t_, RigidBody body_);
	//double TouchGround(double x_[], double p_[], double l_[], double R[], double ta_, double dt_, double epsilon_, RigidBody bodies_[]);
	void TemUpdate(double x_[], double p_[], double l_[], double R[], double y0_[], int i_);



	void Collision(double y0_[], double yfinal_[], RigidBody bodies[], int i_, double* marker_);
	double current_vertex_data[NBODIES * NVertex * 3];

	Contact contacts[NVertex];
private:
	Vector3d force;
	Vector3d torque;
	
};

