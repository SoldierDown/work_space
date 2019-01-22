#pragma once
#include<Eigen/Dense>
#include<Eigen/SVD>
#include<iostream>
#define STATE_SIZE 18
#define GRAVATIY_CONST	9.8
#define NBODIES	2
#define MASS 2
#define NVertex	8
#define TOTAL_RUNNING_TIME 3000
#define TIME_STEP 1
using namespace Eigen;
enum RigidBodyType
{
	STATIC,
	DYNAMIC
};
class RigidBody
{
public:

	RigidBody();
	RigidBody(double px_, double py_, double pz_, double lx_, double ly_, double lz_);
	//RigidBody(double px_, double py_, double pz_, double lx_, double ly_, double lz_, RigidBodyType type_);
	~RigidBody();
	void SetMass(double mass_);
	void SetEdge(double edge_);
	void SetWidth(double width_);
	void SetHeight(double height_);
	void SetDepth(double depth_);

	void SetX(double x_, double y_, double z_);
	void SetRij(int i_,int j_,double r_);
	void SetP(double x_,double y_,double z_);

	void SetL(double x_, double y_, double z_);
	void SetV();
	
	void SetOmega();

	void SetForce();
	void SetTorque();
	

	//Whole
	void SetI_body(Matrix3d I_body_);
	void SetI_body_inv();
	void SetI();
	void SetI_inv();
	

	double		GetMass();

	Vector3d	GetX();
	Vector3d	GetV();
	Matrix3d	GetR();
	Vector3d	GetP();
	Vector3d	GetL();
	Vector3d	GetOmega();
	Vector3d	GetForce();
	Vector3d	GetTorque();
	RigidBodyType GetType();
	Matrix3d	GetI_body();
	Matrix3d	GetI_body_inv();
	Matrix3d	GetI();
	Matrix3d    GetI_inv();
	double*		GetVerteices();
	double		Phi(double x_, double y_, double z_);
private:
	RigidBodyType	type = DYNAMIC;	// rigid body type
	double		mass;				// mass
	double		edge;
	double		width;
	double		height;
	double		depth;
	Vector3d	x;					//x(t)
	Matrix3d	R;					//R(t)
	Vector3d	P, L;				//P(t), L(t)
	Vector3d	v;					//v(t), omega(t)

	Vector3d	force, torque;		//F(t), tau(t)

	// whole
	Matrix3d	I_body, I_body_inv, I, I_inv;
	Vector3d	omega;
	double		vertices[24];
	double		worldspace_vertices[24];
};

