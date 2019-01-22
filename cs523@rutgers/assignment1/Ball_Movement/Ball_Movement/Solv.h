#pragma once
#include"RigidBody.h"
#include<iostream>
class Solv
{
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
		double t_, double t1,RigidBody bodies[]);
	void VertexPosition(RigidBody bodies[], double ori_vertex_data_[], double vertex_data_[], unsigned int n_);
	Matrix3d Star(Vector3d a_);

private:
	Vector3d force;
	Vector3d torque;
};

