#include "Solv.h"

Solv::Solv()
{
	force(0) = 0;
	force(1) = -MASS * GRAVATIY_CONST;
	force(2) = 0;
	torque(0) = 0;
	torque(1) = 0;
	torque(2) = 0;	
}


Solv::~Solv()
{

}

void Solv::State_to_Array(RigidBody *rb, double *y)
{
	*y++ = rb->GetX()(0);
	*y++ = rb->GetX()(1);
	*y++ = rb->GetX()(2);

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
			*y++ = rb->GetR()(i,j);
	}
	
	*y++ = rb->GetP()(0);
	*y++ = rb->GetP()(1);
	*y++ = rb->GetP()(2);



	*y++ = rb->GetL()(0);
	*y++ = rb->GetL()(1);
	*y++ = rb->GetL()(2);

}


void Solv::Array_to_State(RigidBody *rb, double *y)
{	
	double x_, y_, z_, r_;
	x_ = *y++;
	y_ = *y++;
	z_ = *y++;
	rb->SetX(x_, y_, z_);
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			r_ = *y++;
			rb->SetRij(i, j, r_);
		}
	}
	x_ = *y++;
	y_ = *y++;
	z_ = *y++;
	rb->SetP(x_, y_, z_);
	x_ = *y++;
	y_ = *y++;
	z_ = *y++;
	rb->SetL(x_, y_, z_);

	
	rb->SetI_inv();
	rb->SetI();
	rb->SetForce();
	rb->SetTorque();
	rb->SetOmega();
}


void Solv::Array_to_Bodies(RigidBody bodies[], double y[])
{
	for (int i = 0; i < NBODIES; i++)
	{
		Array_to_State(&bodies[i], &y[i*STATE_SIZE]);
	}
}

void Solv::Bodies_to_Array(RigidBody bodies[], double y[])
{
	for (int i = 0; i < NBODIES; i++)
	{
		State_to_Array(&bodies[i], &y[i*STATE_SIZE]);
	}
}



void Solv::Compute_Force_and_Torque(double t, RigidBody *rb)
{
	force = rb->GetForce();
	torque = rb->GetTorque();
}


Matrix3d Solv::Star(Vector3d a_)
{
	Matrix3d a;
	a(0, 0) = 0,		a(0, 1) =-a_[2],		a(0, 2) = a_[1];
	a(1, 0) = a_[2],	a(1, 1) = 0,			a(1, 2) =-a_[0];
	a(2, 0) =-a_[1];	a(2, 1) = a_[0];		a(2, 2) = 0;
	return a;
}

void Solv::ODE(double y0_[], double yfinal_[], int len_,
	double t_, double t1_, RigidBody bodies[])
{
	double dt = 1.0 / 300.0;
	
	for (size_t i = 0; i < NBODIES; i++)
	{
		// update P: 12 13 14
		yfinal_[12 + i * NBODIES] = y0_[12 + i * NBODIES] + dt * force(0);
		yfinal_[13 + i * NBODIES] = y0_[13 + i * NBODIES] + dt * force(1);
		yfinal_[14 + i * NBODIES] = y0_[14 + i * NBODIES] + dt * force(2);

		// update L: 15 16 17
		yfinal_[15 + i * NBODIES] = y0_[15 + i * NBODIES] + dt * torque(0);
		yfinal_[16 + i * NBODIES] = y0_[16 + i * NBODIES] + dt * torque(1);
		yfinal_[17 + i * NBODIES] = y0_[17 + i * NBODIES] + dt * torque(2);
		
		//std::cout << yfinal_[15 + i * NBODIES] << " " << yfinal_[16 + i * NBODIES] << " " << yfinal_[17 + i * NBODIES] << std::endl;


		// update X: 0 1 2
		yfinal_[0 + i * NBODIES] = y0_[0 + i * NBODIES] + dt * yfinal_[12 + i * NBODIES] / bodies[i].GetMass();
		yfinal_[1 + i * NBODIES] = y0_[1 + i * NBODIES] + dt * yfinal_[13 + i * NBODIES] / bodies[i].GetMass();
		yfinal_[2 + i * NBODIES] = y0_[2 + i * NBODIES] + dt * yfinal_[14 + i * NBODIES] / bodies[i].GetMass();
		
	

		// update R: 3 4 5; 6 7 8; 9 10 11;
		Matrix3d Rnm1,Rn;
		//double norm;
		Rnm1(0, 0) = y0_[3], Rnm1(0, 1) = y0_[4], Rnm1(0, 2) = y0_[5];
		Rnm1(1, 0) = y0_[6], Rnm1(1, 1) = y0_[7], Rnm1(1, 2) = y0_[8];
		Rnm1(2, 0) = y0_[9], Rnm1(2, 1) = y0_[10], Rnm1(2, 2) = y0_[11];

		Rn = (Matrix3d::Identity() + dt * Star(bodies[i].GetOmega()))*Rnm1;
		
		//normalize Rn
		Quaterniond q(Rn);
		q.normalize();
		Rn = q.toRotationMatrix();
		
		

		yfinal_[3 + i * NBODIES] = Rn(0, 0);
		yfinal_[6 + i * NBODIES] = Rn(1, 0);	
		yfinal_[9 + i * NBODIES] = Rn(2, 0);
			


		yfinal_[4 + i * NBODIES] = Rn(0, 1);
		
		yfinal_[7 + i * NBODIES] = Rn(1, 1);
		
		yfinal_[10 + i * NBODIES] = Rn(2, 1);
	

		
		yfinal_[5 + i * NBODIES] = Rn(0, 2);
		yfinal_[8 + i * NBODIES] = Rn(1, 2);
		yfinal_[11 + i * NBODIES] = Rn(2, 2);
		

	}
}

void Solv::RunSimulation(RigidBody bodies[], double ori_vertex_data_[], double vertex_data_[])
{
	
	double y0[STATE_SIZE*NBODIES];
	double yfinal[STATE_SIZE*NBODIES];
	Bodies_to_Array(bodies,yfinal);

	//std::cout << TIME_STEP;
	for (int t = 0; t < TOTAL_RUNNING_TIME; t += TIME_STEP)
	{
		for (int i = 0; i < STATE_SIZE*NBODIES; i++)
		{
			y0[i] = yfinal[i];
		}
		
		//std::cout << t << std::endl;
		VertexPosition(bodies, ori_vertex_data_, vertex_data_, t);
			

		ODE(y0, yfinal, STATE_SIZE*NBODIES, t, t + TIME_STEP, bodies);
		

		
		Array_to_Bodies(bodies, yfinal);


		//std::cout.precision(4);
		//std::cout << bodies[0].GetI_body() << std::endl;
		//std::cout << bodies[0].GetI() << std::endl;
		//printf("%lf", bodies[0].GetR());
		//std::cout << std::endl;
		//std::cout << bodies[0].GetR() << std::endl;

		
	}
}
void Solv::VertexPosition(RigidBody bodies[],double ori_vertex_date_[], double vertex_data_[],unsigned int n_)
{
	Vector3d ori_pos;
	Vector3d tra_pos;
	
	for (int i = 0; i < NBODIES; i++) //rigid body #
	{
		//std::cout << bodies[i].GetR() << std::endl;
		for (int j = 0; j < NVertex; j++)// vertex #
		{
			ori_pos = Vector3d(ori_vertex_date_[i*NVertex * 3 + j * 3], ori_vertex_date_[i*NVertex * 3 + j * 3 + 1], ori_vertex_date_[i*NVertex * 3 + j * 3 + 2]);
			
			/*
			std::cout << "Origin Vertex " << j << ": " << std::endl;
			std::cout << "X: " << ori_pos(0) << std::endl;
			std::cout << "Y: " << ori_pos(1) << std::endl;
			std::cout << "Z: " << ori_pos(2) << std::endl;
			*/
			
			tra_pos = bodies[i].GetR()*ori_pos + bodies[i].GetX();
				
			
			/*
			std::cout << "Origin Vertex " << j << ": " << std::endl;
			std::cout << "X: " << ori_pos(0) << std::endl;
			std::cout << "Y: " << ori_pos(1) << std::endl;
			std::cout << "Z: " << ori_pos(2) << std::endl;
			
			std::cout << "Translated Vertex " << j << ": " << std::endl;
			std::cout << "X: " << tra_pos(0) << std::endl;
			std::cout << "Y: " << tra_pos(1) << std::endl;
			std::cout << "Z: " << tra_pos(2) << std::endl;
			*/

			vertex_data_[n_*NBODIES*NVertex * 3 + i * NVertex * 3 + j * 3] = tra_pos(0);
			vertex_data_[n_*NBODIES*NVertex * 3 + i * NVertex * 3 + j * 3 + 1] = tra_pos(1);
			vertex_data_[n_*NBODIES*NVertex * 3 + i * NVertex * 3 + j * 3 + 2] = tra_pos(2);
			
			
			/*check deformation
			if ((j + 1) % 8 == 0)
			{
				
				Vector3d edge1 = { vertex_data_[n_*NBODIES*NVertex * 3 + i * NVertex * 3 + j * 3] - vertex_data_[n_ * NBODIES*NVertex * 3 + i * NVertex * 3 + (j - 2) * 3],
				vertex_data_[n_*NBODIES*NVertex * 3 + i * NVertex * 3 + j * 3 + 1] - vertex_data_[n_ * NBODIES*NVertex * 3 + i * NVertex * 3 + (j - 2) * 3 + 1],
				vertex_data_[n_*NBODIES*NVertex * 3 + i * NVertex * 3 + j * 3 + 2] - vertex_data_[n_ * NBODIES*NVertex * 3 + i * NVertex * 3 + (j - 2) * 3 + 2] };
				
				Vector3d edge2 = { vertex_data_[n_*NBODIES*NVertex * 3 + i * NVertex * 3 + (j-1) * 3] - vertex_data_[n_ * NBODIES*NVertex * 3 + i * NVertex * 3 + (j - 3) * 3],
				vertex_data_[n_*NBODIES*NVertex * 3 + i * NVertex * 3 + (j-1) * 3 + 1] - vertex_data_[n_ * NBODIES*NVertex * 3 + i * NVertex * 3 + (j - 3) * 3 + 1],
				vertex_data_[n_*NBODIES*NVertex * 3 + i * NVertex * 3 + (j-1) * 3 + 2] - vertex_data_[n_ * NBODIES*NVertex * 3 + i * NVertex * 3 + (j - 3) * 3 + 2] };

					
				
				std::cout << "Length: " << sqrt(edge1.dot(edge1)) << " " << sqrt(edge2.dot(edge2)) << std::endl;
			} 
			*/
			
		}

	}


}


