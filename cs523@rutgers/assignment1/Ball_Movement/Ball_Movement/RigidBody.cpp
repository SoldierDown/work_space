#include "RigidBody.h"


RigidBody::RigidBody()
{
}

RigidBody::RigidBody(double px_, double py_, double pz_, double lx_, double ly_, double lz_)
{
	SetMass(1.0);
	SetEdge(1.0);
	SetX(0, 0, 0);

	SetRij(0, 0, 1), SetRij(0, 1, 0), SetRij(0, 2, 0);
	SetRij(1, 0, 0), SetRij(1, 1, 1), SetRij(1, 2, 0);
	SetRij(2, 0, 0), SetRij(2, 1, 0), SetRij(2, 2, 1);

	SetP(px_, py_, pz_);
	SetV();

	SetL(lx_,ly_,lz_);


	Matrix3d I_body_;
	I_body_(0, 0) = mass * edge*edge / 6.0;		I_body_(0, 1) = 0.0;						I_body_(0, 2) = 0.0;
	I_body_(1, 0) = 0.0;						I_body_(1, 1) = mass * edge*edge / 6.0;		I_body_(1, 2) = 0.0;
	I_body_(2, 0) = 0.0;						I_body_(2, 1) = 0.0;						I_body_(2, 2) = mass * edge*edge / 6.0;


	// Set Inertia Matrix
	SetI_body(I_body_);
	SetI_body_inv();
	SetI_inv();
	SetI();


	
	SetForce();
	SetTorque();


	SetOmega();
	
}


RigidBody::~RigidBody()
{

}


//**********************************   Set   *****************************************************************************************
void RigidBody::SetMass(double mass_)
{
	mass = mass_;
}

void RigidBody::SetEdge(double edge_)
{
	edge = edge_;
}

void RigidBody::SetX(double x_, double y_, double z_)
{
	x(0) = x_;
	x(1) = y_;
	x(2) = z_;
}
void RigidBody::SetRij(int i_, int j_, double r_)
{
	R(i_, j_) = r_;

}
void RigidBody::SetP(double x_, double y_, double z_)
{
	P(0) = x_;
	P(1) = y_;
	P(2) = z_;
	
	
}
void RigidBody::SetV()
{
	v(0) = 1.0 / mass * P(0);
	v(1) = 1.0 / mass * P(1);
	v(2) = 1.0 / mass * P(2);
}

void RigidBody::SetL(double x_, double y_, double z_)
{
	L(0) = x_;
	L(1) = y_;
	L(2) = z_;

}

void RigidBody::SetOmega()
{
	omega = I_inv * L;
}

void RigidBody::SetForce()
{
	force(0) = 0;
	force(1) = -mass * GRAVATIY_CONST;
	force(2) = 0;
}
void RigidBody::SetTorque()
{
	torque(0) = 0;
	torque(1) = 0;
	torque(2) = 0;
}

void RigidBody::SetI_body(Matrix3d I_body_)
{
	I_body = I_body_;
}
void RigidBody::SetI_body_inv()
{
	I_body_inv = I_body.inverse();
}
void RigidBody::SetI()
{
	I = I_inv.inverse();
}
void RigidBody::SetI_inv()
{
	I_inv = R * I_body_inv*R.transpose();
}



//**********************************   Get   *****************************************************************************************
double RigidBody::GetMass()
{
	return mass;
}



Vector3d RigidBody::GetX()
{
	return x;
}
Vector3d RigidBody::GetV()
{
	return v;
}
Matrix3d RigidBody::GetR()
{
	return R;
}
Vector3d RigidBody::GetP()
{
	return P;
}
Vector3d RigidBody::GetL()
{
	return L;
}
Vector3d RigidBody::GetOmega()
{
	return omega;
}
Vector3d RigidBody::GetForce()
{
	return force;
}
Vector3d RigidBody::GetTorque()
{
	return torque;
}

Matrix3d RigidBody::GetI()
{
	return I;
}

Matrix3d RigidBody::GetI_inv()
{
	return I_inv;
}

Matrix3d RigidBody::GetI_body()
{
	return I_body;
}
Matrix3d RigidBody::GetI_body_inv()
{
	return I_body_inv;
}