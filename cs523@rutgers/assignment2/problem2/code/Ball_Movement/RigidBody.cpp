#include "RigidBody.h"


RigidBody::RigidBody()
{
	type = STATIC;
	SetMass(9999999999999999999);
	SetWidth(2);
	SetDepth(2);
	SetHeight(0.001);
	SetX(0, 0, 0);
	
	vertices[0] = 0.5*width;	vertices[1] = 0.5*height;		vertices[2] = 0.5*depth;
	vertices[3] = -0.5*width;	vertices[4] = 0.5*height;		vertices[5] = 0.5*depth;
	vertices[6] = -0.5*width;	vertices[7] = -0.5*height;		vertices[8] = 0.5*depth;
	vertices[9] = 0.5*width;	vertices[10] = -0.5*height;		vertices[11] = 0.5*depth;


	vertices[12] = 0.5*width;	vertices[13] = 0.5*height;	vertices[14] = -0.5*depth;
	vertices[15] = -0.5*width;	vertices[16] = 0.5*height;	vertices[17] = -0.5*depth;
	vertices[18] = -0.5*width;	vertices[19] = -0.5*height; vertices[20] = -0.5*depth;
	vertices[21] = 0.5*width;	vertices[22] = -0.5*height; vertices[23] = -0.5*depth;


	worldspace_vertices[0] = x(0) + 0.5*width;		worldspace_vertices[1] = x(1) + 0.5*height;			worldspace_vertices[2] = x(2) + 0.5*depth;
	worldspace_vertices[3] = x(0) - 0.5*width;		worldspace_vertices[4] = x(1) + 0.5*height;			worldspace_vertices[5] = x(2) + 0.5*depth;
	worldspace_vertices[6] = x(0) - 0.5*width;		worldspace_vertices[7] = x(1) - 0.5*height;			worldspace_vertices[8] = x(2) + 0.5*depth;
	worldspace_vertices[9] = x(0) + 0.5*width;		worldspace_vertices[10] = x(1) - 0.5*height;		worldspace_vertices[11] = x(2) + 0.5*depth;


	worldspace_vertices[12] = x(0) + 0.5*width;	worldspace_vertices[13] = x(1) + 0.5*height;		worldspace_vertices[14] = x(2) -0.5*depth;
	worldspace_vertices[15] = x(0) -0.5*width;	worldspace_vertices[16] = x(1) + 0.5*height;		worldspace_vertices[17] = x(2) -0.5*depth;
	worldspace_vertices[18] = x(0) -0.5*width;	worldspace_vertices[19] = x(1) -0.5*height;			worldspace_vertices[20] = x(2) -0.5*depth;
	worldspace_vertices[21] = x(0) + 0.5*width;	worldspace_vertices[22] = x(1) -0.5*height;			worldspace_vertices[23] = x(2) -0.5*depth;


	SetRij(0, 0, 1), SetRij(0, 1, 0), SetRij(0, 2, 0);
	SetRij(1, 0, 0), SetRij(1, 1, 1), SetRij(1, 2, 0);
	SetRij(2, 0, 0), SetRij(2, 1, 0), SetRij(2, 2, 1);

	SetP(0, 0, 0);
	SetV();

	SetL(0, 0, 0);
	//SetL(0, 0, 0);

	Matrix3d I_body_;
	I_body_(0, 0) = mass * height*depth / 6.0;		I_body_(0, 1) = 0.0;						I_body_(0, 2) = 0.0;
	I_body_(1, 0) = 0.0;							I_body_(1, 1) = mass * width*depth / 6.0;		I_body_(1, 2) = 0.0;
	I_body_(2, 0) = 0.0;							I_body_(2, 1) = 0.0;						I_body_(2, 2) = mass * width*height / 6.0;


	// Set Inertia Matrix
	SetI_body(I_body_);
	SetI_body_inv();
	SetI_inv();
	SetI();



	SetForce();
	SetTorque();


	SetOmega();
}

RigidBody::RigidBody(double px_, double py_, double pz_, double lx_, double ly_, double lz_)
{
	type = DYNAMIC;
	SetMass(100);
	SetEdge(0.1);
	SetX(0, .9, 0);
	width = edge;
	height = edge;
	depth = edge;
	vertices[0] = 0.5*width;	vertices[1] = 0.5*height;		vertices[2] = 0.5*depth;
	vertices[3] = -0.5*width;	vertices[4] = 0.5*height;		vertices[5] = 0.5*depth;
	vertices[6] = -0.5*width;	vertices[7] = -0.5*height;		vertices[8] = 0.5*depth;
	vertices[9] = 0.5*width;	vertices[10] = -0.5*height;		vertices[11] = 0.5*depth;


	vertices[12] = 0.5*width;	vertices[13] = 0.5*height;	vertices[14] = -0.5*depth;
	vertices[15] = -0.5*width;	vertices[16] = 0.5*height;	vertices[17] = -0.5*depth;
	vertices[18] = -0.5*width;	vertices[19] = -0.5*height; vertices[20] = -0.5*depth;
	vertices[21] = 0.5*width;	vertices[22] = -0.5*height; vertices[23] = -0.5*depth;


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

void RigidBody::SetWidth(double width_)
{
	width = width_;
}

void RigidBody::SetHeight(double height_)
{
	height = height_;
}

void RigidBody::SetDepth(double depth_)
{
	depth = depth_;
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

RigidBodyType RigidBody::GetType()
{
	return type;
}
double*		RigidBody::GetVerteices()
{
	return vertices;
}

double RigidBody::Phi(double x_, double y_, double z_)
{

		return y_;
	

}