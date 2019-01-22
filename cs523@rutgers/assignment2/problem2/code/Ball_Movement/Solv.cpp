#include "Solv.h"
Solv::Solv()
{
	force(0) = 0;
	force(1) = -GRAVATIY_CONST;
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
void Solv::Bodies_to_Array(RigidBody bodies[], double y[])
{
	for (int i = 0; i < NBODIES; i++)
	{
		State_to_Array(&bodies[i], &y[i*STATE_SIZE]);
	}
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


double* Solv::TouchGround(double x_ta_[], double p_ta_[], double L_ta_[], double R_ta_[], double ta_, double tb_, double epsilon_, RigidBody bodies_[])
{

	double ta = ta_, tb = tb_;
	double x_tb[3], x_tmid[3];
	double p_tb[3], p_tmid[3];
	double L_tb[3], L_tmid[3];
	Matrix3d R_ta, R_tb, R_tmid;
	R_ta(0, 0) = R_ta_[0],	R_ta(0, 1) = R_ta_[1],	R_ta(0, 2) = R_ta_[2];
	R_ta(1, 0) = R_ta_[3],	R_ta(1, 1) = R_ta_[4],	R_ta(1, 2) = R_ta_[5];
	R_ta(2, 0) = R_ta_[6],	R_ta(2, 1) = R_ta_[7],	R_ta(2, 2) = R_ta_[8];
	//double vertex_x_current[24];
	double vertex_x_tb[24];
	//initialize r
	double *r = new double[10]{ 0,0,0,0,0,0,0,0,tb_ - ta_,0 };

	for (size_t i = 0; i < NBODIES; i++)
	{
		if (bodies_[i].GetType() == DYNAMIC)
		{
//------------------------ FIRST STEP: store the VIRTUAL properties of the rigid body ---------------------------------------
			


			// virtual X at tb
			x_tb[0] = x_ta_[0] + tb * p_ta_[0] / bodies_[i].GetMass();
			x_tb[1] = x_ta_[1] + tb * p_ta_[1] / bodies_[i].GetMass();
			x_tb[2] = x_ta_[2] + tb * p_ta_[2] / bodies_[i].GetMass();
			// virtual P:
			p_tb[0] = p_ta_[0] + tb * force(0)* bodies_[i].GetMass();
			p_tb[1] = p_ta_[1] + tb * force(1)* bodies_[i].GetMass();
			p_tb[2] = p_ta_[2] + tb * force(2)* bodies_[i].GetMass();
			// virtual L
			L_tb[0] = L_ta_[0] + tb * torque(0);
			L_tb[1] = L_ta_[1] + tb * torque(1);
			L_tb[2] = L_ta_[2] + tb * torque(2);
			// virtual R: 3 4 5; 6 7 8; 9 10 11;
			
			Eigen::Vector3d L_tb = Eigen::Vector3d(L_ta_[0], L_ta_[1], L_ta_[2]);
			Eigen::Vector3d omega = R_ta * bodies_[i].GetI_body_inv()*R_ta.transpose()*L_tb;

			R_tb = (Matrix3d::Identity() + tb * Star(omega))*R_ta;
			//std::cout << "current R: " << R_current << std::endl;
			//normalize Rn
			Quaterniond q_tb(R_tb);
			q_tb.normalize();
			R_tb = q_tb.toRotationMatrix();




			double min = 10000;
			int marker=-1;
			Eigen::Vector3d ori_pos, tra_pos;
			Eigen::Vector3d center_x_tb = Eigen::Vector3d(x_tb[0], x_tb[1], x_tb[2]);
			int duplicate = 0;
			// find out the vertex with a lowest y component
			for (size_t j = 0; j < NVertex; j++)
			{
				ori_pos = Eigen::Vector3d(bodies_[i].GetVerteices()[j * 3 + 0], bodies_[i].GetVerteices()[j * 3 + 1], bodies_[i].GetVerteices()[j * 3 + 2]);
				tra_pos = R_tb * ori_pos + center_x_tb;
				vertex_x_tb[j * 3 + 0] = tra_pos(0);
				vertex_x_tb[j * 3 + 1] = tra_pos(1);
				vertex_x_tb[j * 3 + 2] = tra_pos(2);

				//std::cout << "Origin Vertex " << j << ": " << std::endl;
				//std::cout << "X: " << ori_pos(0) << std::endl;
				//std::cout << "HHHHHHHHHHHHHHHHHHHHHHH" << std::endl;
				//std::cout << "Y at" << j << ": " << tra_pos(1) << std::endl;
				//std::cout << "Z: " << ori_pos(2) << std::endl;

				if (vertex_x_tb[j * 3 + 1] <= min)
				{
					// more than one vertex: store all the vertex
					if (fabs((vertex_x_tb[j*3+1]-min)/min)<1e-2)
					{
						//std::cout << "MORE THAN ONE" << std::endl;
						duplicate += 1;
						r[j] = 1;
					}
					// just one vertex: update
					else
					{

						duplicate = 1;
						//duplicate += 1;
						if (marker == -1) // first vertex
						{
							//std::cout << "First One" << std::endl;
							min = vertex_x_tb[j * 3 + 1];
							marker = j;
							r[j] = 1;
						}
						else
						{
							//std::cout << "Smaller One" << std::endl;
							min = vertex_x_tb[j * 3 + 1];
							for (size_t k = 0; k < j; k++)
							{
								r[k] = 0;
							}
							marker = j; 
							r[j] = 1;
							
						}

						
						//std::cout << "min is at " << j << ": " << min << std::endl;
					}

				}


			}

			//std::cout << "duplicate: " << duplicate << std::endl;
			r[9] = duplicate;
			//for (size_t i = 0; i < NVertex; i++)
			//{
			//	std::cout << "r[" << i << "]" << r[i] << std::endl;
			//}
			//std::cout << "marker: " << marker << std::endl;
			//std::cout << "min:" << min << std::endl;
			double tmid = 0.5*(ta + tb);

			ori_pos = Eigen::Vector3d(bodies_[i].GetVerteices()[marker * 3 + 0], bodies_[i].GetVerteices()[marker * 3 + 1], bodies_[i].GetVerteices()[marker * 3 + 2]);
			Eigen::Vector3d center_x_ta = Eigen::Vector3d(x_ta_[0], x_ta_[1], x_ta_[2]);
			tra_pos = R_ta * ori_pos + center_x_ta;



			//std::cout << tmid << std::endl;
			//double fai0 = bodies_[i].Phi();
			double fai_ta = bodies_[i].Phi(tra_pos(0), tra_pos(1), tra_pos(2));
			//std::cout << "fai_ta: " << fai_ta << std::endl;

			double fai_tb = bodies_[i].Phi(vertex_x_tb[marker * 3 + 0], vertex_x_tb[marker * 3 + 1], vertex_x_tb[marker * 3 + 2]);
			//std::cout << "fai_tb: " << fai_tb << std::endl;
			//std::cout << "marker: " << marker << std::endl;
			//std::cout << "phi_ta: " << fai_ta << std::endl;
			//std::cout << "phi_tb: " << fai_tb << std::endl;

			if (fai_tb < 0 and fai_ta > 0)
			{

				while (fabs(tb - ta) > epsilon_)
				{
					// update X
					x_tmid[0] = x_ta_[0] + tmid * p_ta_[0] / bodies_[i].GetMass();
					x_tmid[1] = x_ta_[1] + tmid * p_ta_[1] / bodies_[i].GetMass();
					x_tmid[2] = x_ta_[2] + tmid * p_ta_[2] / bodies_[i].GetMass();
					// update P: 12 13 14
					p_tmid[0] = p_ta_[0] + tmid * force(0)* bodies_[i].GetMass();
					p_tmid[1] = p_ta_[1] + tmid * force(1)* bodies_[i].GetMass();
					p_tmid[2] = p_ta_[2] + tmid * force(2)* bodies_[i].GetMass();
					// update L: 15 16 17
					L_tmid[0] = L_ta_[0] + tmid * torque(0);
					L_tmid[1] = L_ta_[1] + tmid * torque(1);
					L_tmid[2] = L_ta_[2] + tmid * torque(2);
					// update R: 3 4 5; 6 7 8; 9 10 11;
					
					//double norm;
					
					R_tmid = (Matrix3d::Identity() + tmid * Star(omega))*R_ta;
					//normalize Rn
					Quaterniond q_tmid(R_tmid);
					q_tmid.normalize();
					R_tmid = q_tmid.toRotationMatrix();

					ori_pos = Eigen::Vector3d(bodies_[i].GetVerteices()[marker * 3 + 0], bodies_[i].GetVerteices()[marker * 3 + 1], bodies_[i].GetVerteices()[marker * 3 + 2]);
					Eigen::Vector3d center_x_tmid = Eigen::Vector3d(x_tb[0], x_tb[1], x_tb[2]);
					tra_pos = R_tb * ori_pos + center_x_tmid;




					if (bodies_[i].Phi(tra_pos(0), tra_pos(1), tra_pos(2)) > 0)
					{
						ta = tmid;
						tmid = 0.5*(ta + tb);
					}
					if (bodies_[i].Phi(tra_pos(0), tra_pos(1), tra_pos(2)) < 0)
					{
						tb = tmid;
						tmid = 0.5*(ta + tb);
					}
					if (bodies_[i].Phi(tra_pos(0), tra_pos(1), tra_pos(2)) == 0)
					{
						r[8] = tmid;
						return r;
					}

					//std::cout << tb - ta << std::endl;


				}
				r[8] = tmid;
				return r;
			}

			else if (fai_ta < 0 && fai_tb < 0)
			{
				std::cout << "WRONG!" << std::endl;
				
				
				return r;
			}


		}

	}

	return r;
}


bool Solv::CollisionCheck(double x0_[], double p0_[], double L0_[], double R0_[], double t_, RigidBody body_)
{
	double dt = t_;
	double x_t[3];
	double p_t[3];
	double L_t[3];
	Matrix3d R0, R_t;
	// R at 0
	R0(0, 0) = R0_[0], R0(0, 1) = R0_[1], R0(0, 2) = R0_[2];
	R0(1, 0) = R0_[3], R0(1, 1) = R0_[4], R0(1, 2) = R0_[5];
	R0(2, 0) = R0_[6], R0(2, 1) = R0_[7], R0(2, 2) = R0_[8];
	//double vertex_x_current[24];
	double min_position[3];
	//initialize r


			//------------------------ FIRST STEP: store the VIRTUAL properties of the rigid body ---------------------------------------
	//std::cout << "dt: " << dt << std::endl;

	// X at t
	x_t[0] = x0_[0] + dt * p0_[0] / body_.GetMass();
	x_t[1] = x0_[1] + dt * p0_[1] / body_.GetMass();
	x_t[2] = x0_[2] + dt * p0_[2] / body_.GetMass();
	// P at t
	p_t[0] = p0_[0] + dt * force(0)* body_.GetMass();
	p_t[1] = p0_[1] + dt * force(1)* body_.GetMass();
	p_t[2] = p0_[2] + dt * force(2)* body_.GetMass();
	// L at t
	L_t[0] = L0_[0] + dt * torque(0);
	L_t[1] = L0_[1] + dt * torque(1);
	L_t[2] = L0_[2] + dt * torque(2);
	
	Eigen::Vector3d L0=Eigen::Vector3d(L0_[0], L0_[1], L0_[2]);
	Eigen::Vector3d omega = R0 * body_.GetI_body_inv()*R0.transpose()*L0;
	

	R_t = (Matrix3d::Identity() + dt * Star(omega))*R0;
	//std::cout << "current R: " << R_current << std::endl;
	//normalize Rn
	Quaterniond q_t(R_t);
	q_t.normalize();
	R_t = q_t.toRotationMatrix();



	double min = 10000;
	Eigen::Vector3d ori_pos, tra_pos;
	Eigen::Vector3d center_x_t = Eigen::Vector3d(x_t[0], x_t[1], x_t[2]);

	// find out the vertex with a lowest y component
	for (size_t j = 0; j < NVertex; j++)
	{
		ori_pos = Eigen::Vector3d(body_.GetVerteices()[j * 3 + 0], body_.GetVerteices()[j * 3 + 1], body_.GetVerteices()[j * 3 + 2]);
		tra_pos = R_t * ori_pos + center_x_t;

		//std::cout << "Origin Vertex " << j << ": " << std::endl;
		//std::cout << "X: " << ori_pos(0) << std::endl;
		//std::cout << "HHHHHHHHHHHHHHHHHHHHHHH" << std::endl;
		//std::cout << "Y at" << j << ": " << tra_pos(1) << std::endl;
		//std::cout << "Z: " << ori_pos(2) << std::endl;

		if (tra_pos[1] < min)
		{
			min = tra_pos(1);
			min_position[0] = tra_pos(0);
			min_position[1] = tra_pos(1);
			min_position[2] = tra_pos(2);
			//std::cout << "min is at " << j << ": " << min << std::endl;
		}


		//for (size_t i = 0; i < NVertex; i++)
		//{
		//	std::cout << "r[" << i << "]" << r[i] << std::endl;
		//}
		//std::cout << "marker: " << marker << std::endl;
		//std::cout << "min:" << min << std::endl;


	}


	//std::cout << tmid << std::endl;
	//double fai0 = bodies_[i].Phi();
	double fai_t = body_.Phi(min_position[0], min_position[1], min_position[2]);

	//std::cout << "fai_t: " << fai_t << std::endl;
	return fai_t < 0 ? true : false;

}


void Solv::ODE(double y0_[], double yfinal_[], int len_,
	double t_, double t1_, RigidBody bodies[])
{

	double dt = t1_ - t_;
	double x_tmp[3];
	double p_tmp[3];
	double l_tmp[3];
	double R_tmp[9];
	for (size_t i = 0; i < NBODIES; i++)
	{
		if (bodies[i].GetType() == DYNAMIC)
		{
			double ddt = 0;
			double dddt = 0;
			double ta = 0.0;
			double tb = dt;
			//std::cout << "y0_[0 + i * STATE_SIZE]: " << y0_[0 + i * STATE_SIZE] << std::endl;
			//std::cout << "y0_[1 + i * STATE_SIZE]: " << y0_[1 + i * STATE_SIZE] << std::endl;
			//std::cout << "y0_[2 + i * STATE_SIZE]: " << y0_[2 + i * STATE_SIZE] << std::endl;

			// we assign the value of each property to xxx_tmp at first
			TemUpdate(x_tmp, p_tmp, l_tmp, R_tmp, y0_, i);
			//for (size_t i = 0; i < 3; i++)
			//{
			//	std::cout << "x_tmp: " << x_tmp[0] << "," << x_tmp[1] << "," << x_tmp[2] << std::endl;
			//	std::cout << "p_tmp: " << p_tmp[0] << "," << p_tmp[1] << "," << p_tmp[2] << std::endl;
			//	std::cout << "l_tmp: " << l_tmp[0] << "," << l_tmp[1] << "," << l_tmp[2] << std::endl;
			//}


			if (!CollisionCheck(x_tmp, p_tmp, l_tmp, R_tmp, tb, bodies[i]))
			{
				PartialUpdate(y0_, yfinal_, bodies, tb, i);

				// Sync
				for (size_t j = 0; j < STATE_SIZE; j++)
				{
					y0_[j + i * STATE_SIZE] = yfinal_[j + i * STATE_SIZE];
				}

			}
			else
			{
				while (CollisionCheck(x_tmp, p_tmp, l_tmp, R_tmp, dt, bodies[i]))
				{
					double* para = TouchGround(x_tmp, p_tmp, l_tmp, R_tmp, ta, tb, 1e-5, bodies);
			//	for (size_t i = 0; i < 9; i++)
			//{
			//	std::cout << "x_tmp: " << x_tmp[0] << "," << x_tmp[1] << "," << x_tmp[2] << std::endl;
			//	std::cout << "para[" << i << "]: " << para[i] << std::endl;
			//	std::cout << "l_tmp: " << l_tmp[0] << "," << l_tmp[1] << "," << l_tmp[2] << std::endl;
			//}
					ddt = para[8];
					dddt = dt - ddt;
					PartialUpdate(y0_, yfinal_, bodies, ddt, i);

					// Sync
					for (size_t j = 0; j < STATE_SIZE; j++)
					
					{
						y0_[j + i * STATE_SIZE] = yfinal_[j + i * STATE_SIZE];
					}

					Collision(y0_, yfinal_, bodies, i, para);
					TemUpdate(x_tmp, p_tmp, l_tmp, R_tmp, y0_, i);
					
					// Sync
					for (size_t j = 0; j < STATE_SIZE; j++)
					{
						y0_[j + i * STATE_SIZE] = yfinal_[j + i * STATE_SIZE];
					}
		
					
					dt = dddt;
					tb = dddt;

				}
				PartialUpdate(y0_, yfinal_, bodies, dddt, i);
			}
			
			

			

			//std::cout << "Before: " << std::endl;
			//for (size_t j = 0; j < STATE_SIZE; j++)
			//{
			//	std::cout << "y0_[" << j << "]: \t" << y0_[j] << std::endl;
			//	std::cout << "yfinal_[" << j << "]:\t" << yfinal_[j] << std::endl;
			//}

			

			//std::cout << "After: " << std::endl;
			//for (size_t j = 0; j < STATE_SIZE; j++)
			//{
			//	std::cout << "y0_[" << j << "]: \t" << y0_[j] << std::endl;
			//	std::cout << "yfinal_[" << j << "]:\t" << yfinal_[j] << std::endl;
			//}


			// assign exact collision time to ddt
			
			//for (size_t i = 0; i < NVertex+1; i++)
			//{
			//	std::cout << "para[" << i << "]: " << para[i] << std::endl;
			//}
			//ddt = dt;
			
			//std::cout << "dt: " << dt << ", ddt: " << ddt << ", dddt: " << dddt << std::endl;


			







		}



		
		if (bodies[i].GetType() == STATIC)
		{

			//PartialUpdate(y0_, yfinal_, bodies, 0, i);
			for (size_t j = 0; j < STATE_SIZE; j++)
			{
				y0_[j + i * STATE_SIZE] = yfinal_[j + i * STATE_SIZE];
			}
			
			//for (size_t i = 0; i < NBODIES* STATE_SIZE; i++)
			//{
			//	std::cout << "y0_[" << i << "]: " << y0_[i] << " yfinal_[" << i << "]: " << yfinal_[i] << std::endl;
			//}
			
			// update P: 12 13 14
			yfinal_[12 + i * STATE_SIZE] = y0_[12 + i * STATE_SIZE];
			yfinal_[13 + i * STATE_SIZE] = y0_[13 + i * STATE_SIZE];
			yfinal_[14 + i * STATE_SIZE] = y0_[14 + i * STATE_SIZE];

			// update L: 15 16 17
			yfinal_[15 + i * STATE_SIZE] = y0_[15 + i * STATE_SIZE];
			yfinal_[16 + i * STATE_SIZE] = y0_[16 + i * STATE_SIZE];
			yfinal_[17 + i * STATE_SIZE] = y0_[17 + i * STATE_SIZE];

			//std::cout << yfinal_[15 + i * STATE_SIZE] << " " << yfinal_[16 + i * STATE_SIZE] << " " << yfinal_[17 + i * STATE_SIZE] << std::endl;


			// update X: 0 1 2
			yfinal_[0 + i * STATE_SIZE] = y0_[0 + i * STATE_SIZE];
			yfinal_[1 + i * STATE_SIZE] = y0_[1 + i * STATE_SIZE];
			yfinal_[2 + i * STATE_SIZE] = y0_[2 + i * STATE_SIZE];


			// update R
			//std::cout << "R" << std::endl;
			//for (int j = 3; j <= 11; j++)
			//{
			//	std::cout << y0_[i*STATE_SIZE + j] << "   " << std::endl;
			//}
			yfinal_[3 + i * STATE_SIZE] = y0_[3 + i * STATE_SIZE];
			yfinal_[6 + i * STATE_SIZE] = y0_[6 + i * STATE_SIZE];
			yfinal_[9 + i * STATE_SIZE] = y0_[9 + i * STATE_SIZE];

			yfinal_[4 + i * STATE_SIZE] = y0_[4 + i * STATE_SIZE];
			yfinal_[7 + i * STATE_SIZE] = y0_[7 + i * STATE_SIZE];
			yfinal_[10 + i * STATE_SIZE] = y0_[10 + i * STATE_SIZE];

			yfinal_[5 + i * STATE_SIZE] = y0_[5 + i * STATE_SIZE];
			yfinal_[8 + i * STATE_SIZE] = y0_[8 + i * STATE_SIZE];
			yfinal_[11 + i * STATE_SIZE] = y0_[11 + i * STATE_SIZE];
			
		}
	}
}


// backup version
/*
void Solv::ODE(double y0_[], double yfinal_[], int len_,
	double t_, double t1_, RigidBody bodies[])
{

	bool is_collide = false;
	double dt = t1_ - t_;
	double x_tmp[3];
	double p_tmp[3];
	double l_tmp[3];
	double R_tmp[9];
	for (size_t i = 0; i < NBODIES; i++)
	{
		if (bodies[i].GetType() == DYNAMIC)
		{
			double ddt = 0;

			//std::cout << "y0_[0 + i * STATE_SIZE]: " << y0_[0 + i * STATE_SIZE] << std::endl;
			//std::cout << "y0_[1 + i * STATE_SIZE]: " << y0_[1 + i * STATE_SIZE] << std::endl;
			//std::cout << "y0_[2 + i * STATE_SIZE]: " << y0_[2 + i * STATE_SIZE] << std::endl;

			x_tmp[0] = y0_[0 + i * STATE_SIZE];
			x_tmp[1] = y0_[1 + i * STATE_SIZE];
			x_tmp[2] = y0_[2 + i * STATE_SIZE];

			p_tmp[0] = y0_[12 + i * STATE_SIZE];
			p_tmp[1] = y0_[13 + i * STATE_SIZE];
			p_tmp[2] = y0_[14 + i * STATE_SIZE];

			l_tmp[0] = y0_[15 + i * STATE_SIZE];
			l_tmp[1] = y0_[16 + i * STATE_SIZE];
			l_tmp[2] = y0_[17 + i * STATE_SIZE];

			R_tmp[0] = y0_[3 + i * STATE_SIZE];
			R_tmp[1] = y0_[4 + i * STATE_SIZE];
			R_tmp[2] = y0_[5 + i * STATE_SIZE];
			R_tmp[3] = y0_[6 + i * STATE_SIZE];
			R_tmp[4] = y0_[7 + i * STATE_SIZE];
			R_tmp[5] = y0_[8 + i * STATE_SIZE];
			R_tmp[6] = y0_[9 + i * STATE_SIZE];
			R_tmp[7] = y0_[10 + i * STATE_SIZE];
			R_tmp[8] = y0_[11 + i * STATE_SIZE];

			double ta = 0.0;
			double tb = dt;


			double* para = TouchGround(x_tmp, p_tmp, l_tmp, R_tmp, 0, dt, 1e-5, bodies);

			//std::cout << "Before: " << std::endl;
			//for (size_t j = 0; j < STATE_SIZE; j++)
			//{
			//	std::cout << "y0_[" << j << "]: \t" << y0_[j] << std::endl;
			//	std::cout << "yfinal_[" << j << "]:\t" << yfinal_[j] << std::endl;
			//}



			//std::cout << "After: " << std::endl;
			//for (size_t j = 0; j < STATE_SIZE; j++)
			//{
			//	std::cout << "y0_[" << j << "]: \t" << y0_[j] << std::endl;
			//	std::cout << "yfinal_[" << j << "]:\t" << yfinal_[j] << std::endl;
			//}


			// assign exact collision time to ddt
			ddt = para[8];
			//for (size_t i = 0; i < NVertex+1; i++)
			//{
			//	std::cout << "para[" << i << "]: " << para[i] << std::endl;
			//}
			//ddt = dt;
			double dddt = dt - ddt;
			//std::cout << "dt: " << dt << ", ddt: " << ddt << ", dddt: " << dddt << std::endl;


			PartialUpdate(y0_, yfinal_, bodies, ddt, i);

			// Sync
			for (size_t j = 0; j < STATE_SIZE; j++)
			{
				y0_[j + i * STATE_SIZE] = yfinal_[j + i * STATE_SIZE];
			}






			// Collision part

			//if (ddt != dt)
			if (CollisionCheck(x_tmp, p_tmp, l_tmp, R_tmp, dt, bodies[i]))
			{

				Collision(y0_, yfinal_, bodies, i, para);
				// Sync
				for (size_t j = 0; j < STATE_SIZE; j++)
				{
					y0_[j + i * STATE_SIZE] = yfinal_[j + i * STATE_SIZE];
				}

				// after collision
				// update X

				PartialUpdate(y0_, yfinal_, bodies, dddt, i);



			}



		}
		if (bodies[i].GetType() == STATIC)
		{
			// update P: 12 13 14
			yfinal_[12 + i * STATE_SIZE] = y0_[12 + i * STATE_SIZE];
			yfinal_[13 + i * STATE_SIZE] = y0_[13 + i * STATE_SIZE];
			yfinal_[14 + i * STATE_SIZE] = y0_[14 + i * STATE_SIZE];

			// update L: 15 16 17
			yfinal_[15 + i * STATE_SIZE] = y0_[15 + i * STATE_SIZE];
			yfinal_[16 + i * STATE_SIZE] = y0_[16 + i * STATE_SIZE];
			yfinal_[17 + i * STATE_SIZE] = y0_[17 + i * STATE_SIZE];

			//std::cout << yfinal_[15 + i * STATE_SIZE] << " " << yfinal_[16 + i * STATE_SIZE] << " " << yfinal_[17 + i * STATE_SIZE] << std::endl;


			// update X: 0 1 2
			yfinal_[0 + i * STATE_SIZE] = y0_[0 + i * STATE_SIZE];
			yfinal_[1 + i * STATE_SIZE] = y0_[1 + i * STATE_SIZE];
			yfinal_[2 + i * STATE_SIZE] = y0_[2 + i * STATE_SIZE];


			// update R
			//std::cout << "R" << std::endl;
			//for (int j = 3; j <= 11; j++)
			//{
			//	std::cout << y0_[i*STATE_SIZE + j] << "   " << std::endl;
			//}
			yfinal_[3 + i * STATE_SIZE] = y0_[3 + i * STATE_SIZE];
			yfinal_[6 + i * STATE_SIZE] = y0_[6 + i * STATE_SIZE];
			yfinal_[9 + i * STATE_SIZE] = y0_[9 + i * STATE_SIZE];

			yfinal_[4 + i * STATE_SIZE] = y0_[4 + i * STATE_SIZE];
			yfinal_[7 + i * STATE_SIZE] = y0_[7 + i * STATE_SIZE];
			yfinal_[10 + i * STATE_SIZE] = y0_[10 + i * STATE_SIZE];

			yfinal_[5 + i * STATE_SIZE] = y0_[5 + i * STATE_SIZE];
			yfinal_[8 + i * STATE_SIZE] = y0_[8 + i * STATE_SIZE];
			yfinal_[11 + i * STATE_SIZE] = y0_[11 + i * STATE_SIZE];
		}
	}
}

*/
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
			//std::cout << "y0_["<<i<<"]: " << y0[i] << std::endl;
			
		}

		//std::cout << t << std::endl;
		VertexPosition(bodies, ori_vertex_data_, vertex_data_, t);
			
		
		ODE(y0, yfinal, STATE_SIZE*NBODIES, 0, 1.0 / 300.0, bodies);
		


		
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
				
			
			
			//std::cout << "Origin Vertex " << j << ": " << std::endl;
			//std::cout << "X: " << ori_pos(0) << std::endl;
			//std::cout << "Y: " << ori_pos(1) << std::endl;
			//std::cout << "Z: " << ori_pos(2) << std::endl;
			//
			//std::cout << "Translated Vertex " << j << ": " << std::endl;
			//std::cout << "X: " << tra_pos(0) << std::endl;
			//std::cout << "Y: " << tra_pos(1) << std::endl;
			//std::cout << "Z: " << tra_pos(2) << std::endl;
			

			vertex_data_[n_*NBODIES*NVertex * 3 + i * NVertex * 3 + j * 3] = tra_pos(0);
			vertex_data_[n_*NBODIES*NVertex * 3 + i * NVertex * 3 + j * 3 + 1] = tra_pos(1);
			vertex_data_[n_*NBODIES*NVertex * 3 + i * NVertex * 3 + j * 3 + 2] = tra_pos(2);
			
			current_vertex_data[i * NVertex * 3 + j * 3] = tra_pos(0);
			current_vertex_data[i * NVertex * 3 + j * 3 + 1] = tra_pos(1);
			current_vertex_data[i * NVertex * 3 + j * 3 + 2] = tra_pos(2);
			
		}

	}


}

void Solv::Collision(double y0_[], double yfinal_[], RigidBody bodies_[], int i_, double* marker_)
{

	double epsilon = 0.05;
	Eigen::Vector3d n = Eigen::Vector3d(0, 1, 0);
	Eigen::Vector3d ra;
	Eigen::Vector3d dl;
	Eigen::Vector3d vrel;
	double	numerator;
	double	term1, term2, jtotal;

	if (marker_[9] == 4) // face
	{
		yfinal_[13 + i_ * STATE_SIZE] = -epsilon * y0_[13 + i_ * STATE_SIZE];
		y0_[13 + i_ * STATE_SIZE] = yfinal_[13 + i_ * STATE_SIZE];
	}
	else   // edge or vertex
	{
		for (size_t j = 0; j < NVertex; j++)
		{
			if (marker_[j] == 1)
			{
				//std::cout << "Apply to index " << j << std::endl;
				ra = bodies_[i_].GetR()*Eigen::Vector3d(bodies_[i_].GetVerteices()[j * 3 + 0], bodies_[i_].GetVerteices()[j * 3 + 1], bodies_[i_].GetVerteices()[j * 3 + 2]);

				//std::cout << "relative r: " << (bodies_[i_].GetR()*Eigen::Vector3d(bodies_[i_].GetVerteices()[marker_ * 3 + 0], bodies_[i_].GetVerteices()[marker_ * 3 + 1], bodies_[i_].GetVerteices()[marker_ * 3 + 2])).norm() << std::endl;
				//Eigen::Vector3d padot = pt_velocity(&bodies[i_], bodies[i_].GetX());

					//ra = c->p - c->a->GetX(),
					//rb = c->p - c->b->GetX();


				//double vrel = yfinal_[13 + i_ * STATE_SIZE] / bodies_[i_].GetMass(),

				vrel = bodies_[i_].GetV() + bodies_[i_].GetOmega().cross(ra);
				//std::cout << "vrel.y: " << vrel(1) << std::endl;

				numerator = -(1 + epsilon)*vrel(1);
				term1 = 1 / bodies_[i_].GetMass();
				term2 = n.dot(((bodies_[i_].GetI_inv()*(ra.cross(n))).cross(ra)));
				jtotal = numerator / (term1 + term2);
				//double	jtotal = numerator / term1;

				//std::cout << "jtotal: " << jtotal << std::endl;
				//std::cout << "ra: " << ra.norm() << std::endl;
				//std::cout << "numerator: " << numerator << std::endl;
				//std::cout << "term1: " << term1 << std::endl;
				//std::cout << "term2: " << term2 << std::endl;

				Eigen::Vector3d impulse_tau = jtotal * n;
				dl = ra.cross(impulse_tau);
				//std::cout << "jtotal: " << jtotal << ", dl: (" << dl(0) << "," << dl(1) << "," << dl(2) << ")" << std::endl;

				//std::cout << "py before: " << yfinal_[13 + i_ * STATE_SIZE] << " or " << y0_[13 + i_ * STATE_SIZE] << std::endl;
				yfinal_[13 + i_ * STATE_SIZE] += jtotal;
				y0_[13 + i_ * STATE_SIZE] += jtotal;
				//std::cout << "py after: " << yfinal_[13 + i_ * STATE_SIZE] << " or " << y0_[13 + i_ * STATE_SIZE] << std::endl;
				yfinal_[15 + i_ * STATE_SIZE] += dl(0);
				y0_[15 + i_ * STATE_SIZE] += dl(0);
				yfinal_[16 + i_ * STATE_SIZE] += dl(1);
				y0_[16 + i_ * STATE_SIZE] += dl(1);
				yfinal_[17 + i_ * STATE_SIZE] += dl(2);
				y0_[17 + i_ * STATE_SIZE] += dl(2);
			}

		}


	}

	

	/*
	Eigen::Vector3d padot = pt_velocity(c->a, c->p),
		pbdot = pt_velocity(c->b, c->p),
		n = Eigen::Vector3d()
		ra = c->p - c->a->GetX(),
		rb = c->p - c->b->GetX();


	double vrel = c->n.dot(padot - pbdot),
		numerator = -(1 + epsilon)*vrel;

	double term1 = 1 / c->a->GetMass(),
		term2 = n.dot(((c->a->GetI_inv()*(ra.cross(n))).cross(ra)));



	double j = numerator / (term1 + term2);
	Eigen::Vector3d force = j * n;
	


 ; // = y0_[13 + i * STATE_SIZE] + dddt * force(1)* bodies[i].GetMass();


	// update L: 15 16 17
	*/




}


void Solv::PartialUpdate(double y0_[], double yfinal_[],RigidBody bodies_[], double dt_, int i_)
{
	// update X: 0 1 2
	yfinal_[0 + i_ * STATE_SIZE] = y0_[0 + i_ * STATE_SIZE] + dt_ * y0_[12 + i_ * STATE_SIZE] / bodies_[i_].GetMass();
	yfinal_[1 + i_ * STATE_SIZE] = y0_[1 + i_ * STATE_SIZE] + dt_ * y0_[13 + i_ * STATE_SIZE] / bodies_[i_].GetMass();
	yfinal_[2 + i_ * STATE_SIZE] = y0_[2 + i_ * STATE_SIZE] + dt_ * y0_[14 + i_ * STATE_SIZE] / bodies_[i_].GetMass();
	// update P: 12 13 14
	yfinal_[12 + i_ * STATE_SIZE] = y0_[12 + i_ * STATE_SIZE] + dt_ * force(0)* bodies_[i_].GetMass();
	yfinal_[13 + i_ * STATE_SIZE] = y0_[13 + i_ * STATE_SIZE] + dt_ * force(1)* bodies_[i_].GetMass();
	yfinal_[14 + i_ * STATE_SIZE] = y0_[14 + i_ * STATE_SIZE] + dt_ * force(2)* bodies_[i_].GetMass();
	// update L: 15 16 17
	yfinal_[15 + i_ * STATE_SIZE] = y0_[15 + i_ * STATE_SIZE] + dt_ * torque(0);
	yfinal_[16 + i_ * STATE_SIZE] = y0_[16 + i_ * STATE_SIZE] + dt_ * torque(1);
	yfinal_[17 + i_ * STATE_SIZE] = y0_[17 + i_ * STATE_SIZE] + dt_ * torque(2);




	// update R: 3 4 5; 6 7 8; 9 10 11;
	Matrix3d Rnm1, Rn;
	//double norm;
	Rnm1(0, 0) = y0_[3],	Rnm1(0, 1) = y0_[4],	 Rnm1(0, 2) = y0_[5];
	Rnm1(1, 0) = y0_[6],	Rnm1(1, 1) = y0_[7],	Rnm1(1, 2) = y0_[8];
	Rnm1(2, 0) = y0_[9],	Rnm1(2, 1) = y0_[10],	Rnm1(2, 2) = y0_[11];

	Rn = (Matrix3d::Identity() + dt_ * Star(bodies_[i_].GetOmega()))*Rnm1;

	//normalize Rn
	Quaterniond q(Rn);
	q.normalize();
	Rn = q.toRotationMatrix();



	yfinal_[3 + i_ * STATE_SIZE] = Rn(0, 0);
	yfinal_[6 + i_ * STATE_SIZE] = Rn(1, 0);
	yfinal_[9 + i_ * STATE_SIZE] = Rn(2, 0);



	yfinal_[4 + i_ * STATE_SIZE] = Rn(0, 1);

	yfinal_[7 + i_ * STATE_SIZE] = Rn(1, 1);

	yfinal_[10 + i_ * STATE_SIZE] = Rn(2, 1);



	yfinal_[5 + i_ * STATE_SIZE] = Rn(0, 2);
	yfinal_[8 + i_ * STATE_SIZE] = Rn(1, 2);
	yfinal_[11 + i_ * STATE_SIZE] = Rn(2, 2);


}

void Solv::TemUpdate(double x_tmp_[], double p_tmp_[], double L_tmp_[], double R_tmp_[], double y0_[], int i_)
{
	x_tmp_[0] = y0_[0 + i_ * STATE_SIZE];
	x_tmp_[1] = y0_[1 + i_ * STATE_SIZE];
	x_tmp_[2] = y0_[2 + i_ * STATE_SIZE];

	p_tmp_[0] = y0_[12 + i_ * STATE_SIZE];
	p_tmp_[1] = y0_[13 + i_ * STATE_SIZE];
	p_tmp_[2] = y0_[14 + i_ * STATE_SIZE];

	L_tmp_[0] = y0_[15 + i_ * STATE_SIZE];
	L_tmp_[1] = y0_[16 + i_ * STATE_SIZE];
	L_tmp_[2] = y0_[17 + i_ * STATE_SIZE];

	R_tmp_[0] = y0_[3 + i_ * STATE_SIZE];
	R_tmp_[1] = y0_[4 + i_ * STATE_SIZE];
	R_tmp_[2] = y0_[5 + i_ * STATE_SIZE];
	R_tmp_[3] = y0_[6 + i_ * STATE_SIZE];
	R_tmp_[4] = y0_[7 + i_ * STATE_SIZE];
	R_tmp_[5] = y0_[8 + i_ * STATE_SIZE];
	R_tmp_[6] = y0_[9 + i_ * STATE_SIZE];
	R_tmp_[7] = y0_[10 + i_ * STATE_SIZE];
	R_tmp_[8] = y0_[11 + i_ * STATE_SIZE];



}