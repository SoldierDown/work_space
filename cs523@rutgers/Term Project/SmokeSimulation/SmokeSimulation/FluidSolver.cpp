
#include "FluidSolver.h"

//----------------------------- Check Input type of ScalarIndex-----------------------------------------------------------------------
FluidSolver::FluidSolver()
{
	Init();

}

FluidSolver::~FluidSolver()
{

}

void FluidSolver::Init()
{
	for (size_t i = 0; i < N + 2; i++)
	{

		for (size_t j = 0; j < N + 2; j++)
		{

			if (Distance(i*cell_width, j*cell_width) < NON_ZERO_REGION*cell_width)
			{
				density_latest[ScalarIndex(i, j)] = 1.0;
			}
			else
			{
				density_latest[ScalarIndex(i, j)] = 0.0;
			}


			//if (i <= NON_ZERO_REGION && abs(j - (N + 2) / 2.0) < NON_ZERO_REGION / 2.0)
			//{
			//	density_latest[ScalarIndex(i, j)] = 1;
			//}
			//else
			//{
			//	density_latest[ScalarIndex(i, j)] = 0;
			//}

			divV[ScalarIndex(i, j)] = 0.0;
			pressure_bar[ScalarIndex(i, j)] = 0.0;
		}
	}

	for (size_t i = 0; i < N + 1; i++)
	{
		for (size_t j = 0; j < N + 2; j++)
		{

			if (Distance((i + 0.5)*cell_width, j*cell_width) < NON_ZERO_REGION*cell_width)
			{
				vx_latest[VxIndex(i, j)] = -Distance(j*cell_width, (i + 0.5)*cell_width) * Sin((i + 0.5)*cell_width, j*cell_width);
			}
			else
			{
				vx_latest[VxIndex(i, j)] = 0.0;
			}


			if (Distance(j*cell_width, (i + 0.5)*cell_width) < NON_ZERO_REGION*cell_width)
			{
				vy_latest[VyIndex(j, i)] = Distance(j*cell_width, (i + 0.5)*cell_width) * Cos(j*cell_width, (i + 0.5)*cell_width);
				
			}
			else
			{
				vy_latest[VyIndex(j, i)] = 0.0;
			}




		}

	}

	SyncV();
	SyncDensity();

}

// for property stored in the center of a cell
// so far so good
int FluidSolver::ScalarIndex(int i_, int j_)
{
	return j_ * (N + 2) + i_;
}
// so far so good
int FluidSolver::VxIndex(int i_, int j_)
{
	return j_ * (N + 1) + i_;
}
// so far so good
int FluidSolver::VyIndex(int i_, int j_)
{
	return j_ * (N + 2) + i_;
}
// so far so good
void FluidSolver::DensityAdvect()
{
	double x, y;
	double vx, vy;
	double dx, dy;
	int i0, i1, j0, j1;
	int vx_i, vx_j;
	int vy_i, vy_j;
	for (size_t i = 1; i <= N; i++)
	{
		for (size_t j = 1; j <= N; j++)
		{
			// vx= ...;
			vx_i = i - 1;
			vx_j = j;
			vx = 0.5*(vx_latest[VxIndex(vx_i, vx_j)] + vx_latest[VxIndex(vx_i + 1, vx_j)]);
			// vy= ...;
			vy_i = i;
			vy_j = j - 1;
			vy = 0.5*(vy_latest[VyIndex(vy_i, vy_j)] + vy_latest[VyIndex(vy_i, vy_j + 1)]);

			x = i * cell_width - dt * vx;
			y = j * cell_width - dt * vy;

			if (x < 0) x = 0;
			if (x > (N + 1)*cell_width) x = (N + 1)*cell_width;
			i0 = floor(x / cell_width);	i1 = i0 + 1;
			dx = x - i0 * cell_width;
			dx = dx / cell_width;

			if (y < 0) y = 0;
			if (y > (N + 1)*cell_width) y = (N + 1)*cell_width;
			j0 = floor(y / cell_width);	j1 = j0 + 1;
			dy = y - j0 * cell_width;
			dy = dy / cell_width;


			density_latest[ScalarIndex(i, j)] = dy * (dx * density_backup[ScalarIndex(i1, j1)] + (1 - dx)*density_backup[ScalarIndex(i0, j1)])
				+ (1 - dy)*(dx * density_backup[ScalarIndex(i1, j0)] + (1 - dx)*density_backup[ScalarIndex(i0, j0)]);
		
		}

	}
	SetDensityBoundary();

}
// so far so good
void FluidSolver::VelocityAdvect()
{

	int interpolate_i;
	int interpolate_j;
	double vx_interpolate;
	double vy_interpolate;

	double x, y;
	double vx = 0.0, vy = 0.0;
	double dx, dy;
	int i0, i1, j0, j1;


	// advect vx
	for (size_t i = 1; i <= N - 1; i++)
	{
		for (size_t j = 1; j <= N; j++)
		{
			interpolate_i = i;
			interpolate_j = j - 1;
			vy_interpolate = 0.25*(vy_backup[VyIndex(interpolate_i, interpolate_j)] + vy_backup[VyIndex(interpolate_i + 1, interpolate_j)]
				+ vy_backup[VyIndex(interpolate_i, interpolate_j + 1)] + vy_backup[VyIndex(interpolate_i + 1, interpolate_j + 1)]);


			//vy_backup[VyIndex(interpolate_i, interpolate_j)];
			//std::cout << "vy_interpolate: " << vy_interpolate << std::endl;

			x = (i + 0.5) * cell_width - dt * vx_backup[VxIndex(i, j)];
			y = j * cell_width - dt * vy_interpolate;

			if (x < 0.5*cell_width) x = 0.5*cell_width;
			if (x > (0.5 + N)*cell_width) x = (0.5 + N)*cell_width;
			i0 = floor((x - 0.5*cell_width) / cell_width);	i1 = i0 + 1;
			dx = x - (i0 + 0.5) * cell_width;
			dx = dx / cell_width;

			if (y < 0) y = 0;
			if (y > (N + 1)*cell_width) y = (N + 1)*cell_width;
			j0 = floor(y / cell_width);	j1 = j0 + 1;
			dy = y - j0 * cell_width;
			dy = dy / cell_width;
			//std::cout << "dx: " << dx << " dy: " << dy << std::endl;

			//std::cout << "advect vx dx: " << dx << ", dy: " << dy << std::endl;
			vx_latest[VxIndex(i, j)] = dy * (dx * vx_backup[VxIndex(i1, j1)] + (1 - dx)*vx_backup[VxIndex(i0, j1)])
				+ (1 - dy)*(dx * vx_backup[VxIndex(i1, j0)] + (1 - dx)*vx_backup[VxIndex(i0, j0)]);


		}

	}




	// advect vy
	for (size_t i = 1; i <= N; i++)
	{
		for (size_t j = 1; j <= N - 1; j++)
		{
			interpolate_i = i - 1;
			interpolate_j = j;
			vx_interpolate = 0.25*(vx_backup[VxIndex(interpolate_i, interpolate_j)] + vx_backup[VxIndex(interpolate_i + 1, interpolate_j)]
				+ vx_backup[VxIndex(interpolate_i, interpolate_j + 1)] + vx_backup[VxIndex(interpolate_i + 1, interpolate_j + 1)]);



			x = i * cell_width - dt * vx_interpolate;
			y = (j + 0.5) * cell_width - dt * vy_backup[VyIndex(i, j)];
			if (x < 0) x = 0;
			if (x > (N + 1)*cell_width) x = (N + 1)*cell_width;
			i0 = floor(x / cell_width);	i1 = i0 + 1;
			dx = x - i0 * cell_width;
			dx = dx / cell_width;

			if (y < 0.5*cell_width) y = 0;
			if (y > (N + 0.5)*cell_width) y = (N + 0.5)*cell_width;
			j0 = floor((y - 0.5*cell_width) / cell_width);	j1 = j0 + 1;
			dy = y - (j0 + 0.5) * cell_width;
			dy = dy / cell_width;

			//std::cout << "advect vy dx: " << dx << ", dy: " << dy << std::endl;
			vy_latest[VyIndex(i, j)] = dy * (dx * vy_backup[VyIndex(i1, j1)] + (1 - dx)*vy_backup[VyIndex(i0, j1)])
				+ (1 - dy)*(dx * vy_backup[VyIndex(i1, j0)] + (1 - dx)*vy_backup[VyIndex(i0, j0)]);


		}

	}


}
// so far so good
void FluidSolver::SetDensityBoundary()
{

	for (size_t i = 1; i <= N; i++)
	{
		density_latest[ScalarIndex(i, 0)] = density_latest[ScalarIndex(i, 1)];
		density_latest[ScalarIndex(i, N + 1)] = density_latest[ScalarIndex(i, N)];

		density_latest[ScalarIndex(0, i)] = density_latest[ScalarIndex(1, i)];
		density_latest[ScalarIndex(N + 1, i)] = density_latest[ScalarIndex(N, i)];
	}

	// at corner
	density_latest[ScalarIndex(0, 0)] = 0.5*(density_latest[ScalarIndex(1, 0)] + density_latest[ScalarIndex(0, 1)]);
	density_latest[ScalarIndex(0, N + 1)] = 0.5*(density_latest[ScalarIndex(1, N + 1)] + density_latest[ScalarIndex(0, N)]);
	density_latest[ScalarIndex(N + 1, 0)] = 0.5*(density_latest[ScalarIndex(N, 0)] + density_latest[ScalarIndex(N + 1, 1)]);
	density_latest[ScalarIndex(N + 1, N + 1)] = 0.5*(density_latest[ScalarIndex(N, N + 1)] + density_latest[ScalarIndex(N + 1, N)]);



}
// so far so good
void FluidSolver::SetPressureBoundary()
{

	for (size_t i = 1; i <= N; i++)
	{
		pressure_bar[ScalarIndex(i, 0)] = pressure_bar[ScalarIndex(i, 1)];
		pressure_bar[ScalarIndex(i, N + 1)] = pressure_bar[ScalarIndex(i, N)];

		pressure_bar[ScalarIndex(0, i)] = pressure_bar[ScalarIndex(1, i)];
		pressure_bar[ScalarIndex(N + 1, i)] = pressure_bar[ScalarIndex(N, i)];
	}

	// at corner
	pressure_bar[ScalarIndex(0, 0)] = 0.5*(pressure_bar[ScalarIndex(1, 0)] + pressure_bar[ScalarIndex(0, 1)]);
	pressure_bar[ScalarIndex(0, N + 1)] = 0.5*(pressure_bar[ScalarIndex(1, N + 1)] + pressure_bar[ScalarIndex(0, N)]);
	pressure_bar[ScalarIndex(N + 1, 0)] = 0.5*(pressure_bar[ScalarIndex(N, 0)] + pressure_bar[ScalarIndex(N + 1, 1)]);
	pressure_bar[ScalarIndex(N + 1, N + 1)] = 0.5*(pressure_bar[ScalarIndex(N, N + 1)] + pressure_bar[ScalarIndex(N + 1, N)]);

}
// so far so good
void FluidSolver::SetDivVBoundary()
{

	for (size_t i = 1; i <= N; i++)
	{
		divV[ScalarIndex(i, 0)] = divV[ScalarIndex(i, 1)];
		divV[ScalarIndex(i, N + 1)] = divV[ScalarIndex(i, N)];

		divV[ScalarIndex(0, i)] = divV[ScalarIndex(1, i)];
		divV[ScalarIndex(N + 1, i)] = divV[ScalarIndex(N, i)];
	}

	// at corner
	divV[ScalarIndex(0, 0)] = 0.5*(divV[ScalarIndex(1, 0)] + divV[ScalarIndex(0, 1)]);
	divV[ScalarIndex(0, N + 1)] = 0.5*(divV[ScalarIndex(1, N + 1)] + divV[ScalarIndex(0, N)]);
	divV[ScalarIndex(N + 1, 0)] = 0.5*(divV[ScalarIndex(N, 0)] + divV[ScalarIndex(N + 1, 1)]);
	divV[ScalarIndex(N + 1, N + 1)] = 0.5*(divV[ScalarIndex(N, N + 1)] + divV[ScalarIndex(N + 1, N)]);

}

void FluidSolver::SetVelocityBoundary()
{
	for (size_t i = 1; i <= N - 1; i++)
	{
		vx_latest[VxIndex(i, 0)] = vx_latest[VxIndex(i, 1)];
		vx_latest[VxIndex(i, N + 1)] = vx_latest[VxIndex(i, N)];

		vy_latest[VyIndex(0, i)] = vy_latest[VyIndex(1, i)];
		vy_latest[VyIndex(N + 1, i)] = vy_latest[VyIndex(N, i)];
	}

	for (size_t i = 1; i <= N; i++)
	{
		vx_latest[VxIndex(0, i)] = 0.0;
		vx_latest[VxIndex(N, i)] = 0.0;

		vy_latest[VyIndex(i, 0)] = 0.0;
		vy_latest[VyIndex(i, N)] = 0.0;

	}

	// at corners
	vx_latest[VxIndex(0, 0)] = 0.5*(vx_latest[VxIndex(1, 0)] + vx_latest[VxIndex(0, 1)]);
	vx_latest[VxIndex(N, 0)] = 0.5*(vx_latest[VxIndex(N - 1, 0)] + vx_latest[VxIndex(N, 1)]);
	vx_latest[VxIndex(0, N + 1)] = 0.5*(vx_latest[VxIndex(0, N)] + vx_latest[VxIndex(1, N + 1)]);
	vx_latest[VxIndex(N, N + 1)] = 0.5*(vx_latest[VxIndex(N - 1, N + 1)] + vx_latest[VxIndex(N, N)]);


	vy_latest[VyIndex(0, 0)] = 0.5*(vy_latest[VyIndex(1, 0)] + vy_latest[VyIndex(0, 1)]);
	vy_latest[VyIndex(N + 1, 0)] = 0.5*(vy_latest[VyIndex(N, 0)] + vy_latest[VyIndex(N + 1, 1)]);
	vy_latest[VyIndex(0, N)] = 0.5*(vy_latest[VyIndex(0, N - 1)] + vy_latest[VyIndex(1, N)]);
	vy_latest[VyIndex(N + 1, N)] = 0.5*(vy_latest[VyIndex(N, N)] + vy_latest[VyIndex(N + 1, N - 1)]);




}

void FluidSolver::Project()
{


	for (size_t i = 1; i <= N; i++)
	{
		for (size_t j = 1; j <= N; j++)
		{
			divV[ScalarIndex(i, j)] = -(vx_latest[VxIndex(i, j)] - vx_latest[VxIndex(i - 1, j)]
				+ vy_latest[VyIndex(i, j)] - vy_latest[VyIndex(i, j - 1)]) / cell_width;
			pressure_bar[ScalarIndex(i, j)] = 0.0;
		}
	}
	SetDivVBoundary();
	SetPressureBoundary();


	for (size_t k = 0; k < 20; k++)
	{
		for (size_t i = 1; i <= N; i++)
		{
			for (size_t j = 1; j <= N; j++)
			{
				pressure_bar[ScalarIndex(i, j)] = (cell_width*cell_width*divV[ScalarIndex(i, j)]
					+ pressure_bar[ScalarIndex(i - 1, j)] + pressure_bar[ScalarIndex(i + 1, j)]
					+ pressure_bar[ScalarIndex(i, j - 1)] + pressure_bar[ScalarIndex(i, j + 1)]) *0.25;
			}
		}

		SetPressureBoundary();

	}


	for (size_t i = 1; i <= N - 1; i++)
	{
		for (size_t j = 1; j <= N; j++)
		{
			// interpolate pressure to vx, vy here 
			vx_latest[VxIndex(i, j)] -= (pressure_bar[ScalarIndex(i + 1, j)] - pressure_bar[ScalarIndex(i, j)]) / cell_width;

			vy_latest[VyIndex(j, i)] -= (pressure_bar[ScalarIndex(j, i + 1)] - pressure_bar[ScalarIndex(j, i)]) / cell_width;

		}
	}
	SetVelocityBoundary();





}

void FluidSolver::SyncV()
{
	for (size_t i = 0; i < (N + 1)*(N + 2); i++)
	{
		vx_backup[i] = vx_latest[i];
		vy_backup[i] = vy_latest[i];
	}

}

void FluidSolver::SyncDensity()
{

	for (size_t i = 0; i < (N + 2)*(N + 2); i++)
	{
		density_backup[i] = density_latest[i];
		//std::cout << "thread_num: " << omp_get_thread_num() << std::endl;
	}
}

void FluidSolver::Save(double density_data_[])
{

	for (size_t i = 0; i < (N + 2)*(N + 2); i++)
	{

		density_data_[i] = density_latest[i];
	}
}

double FluidSolver::Distance(double x_, double y_)
{
	return sqrt((x_ - (N + 2) / 2 * cell_width)*(x_ - (N + 2) / 2 * cell_width) + (y_ - (N + 2) / 2 * cell_width)*(y_ - (N + 2) / 2 * cell_width));
}

double FluidSolver::Sin(double x_, double y_)
{
	if (Distance(x_, y_) == 0)
	{
		return 0.0;
	}
	else
	{
		return (y_ - (N + 2) / 2 * cell_width) / Distance(x_, y_);
	}

}
double FluidSolver::Cos(double x_, double y_)
{

	if (Distance(x_, y_) == 0)
	{
		return 0.0;
	}
	else
	{
		return (x_ - (N + 2) / 2 * cell_width) / Distance(x_, y_);
	}

}

void FluidSolver::Update()
{
	SetPermantSource();
	SyncV();
	SyncDensity();
	DensityDiffuse();

	VelocityDiffuse();

	Project();

	SyncV();

	SyncDensity();

	DensityAdvect();

	VelocityAdvect();

	SyncV();

	SyncDensity();

	//total = 0.0;
	//for (size_t i = 0; i < (N+2)*(N+2); i++)
	//{
	//	total += density_latest[i];
	//}
	//std::cout << "total num: " << total << std::endl;

}

void FluidSolver::DensityDiffuse()
{
	double diff = 0.1;
	double a = dt * diff / (cell_width*cell_width);
	for (size_t k = 0; k < 20; k++) {
		for (size_t i = 1; i <= N; i++) {
			for (size_t j = 1; j <= N; j++) {
				density_latest[ScalarIndex(i, j)] = (density_backup[ScalarIndex(i, j)] + a * (density_latest[ScalarIndex(i - 1, j)] + density_latest[ScalarIndex(i + 1, j)] +
					density_latest[ScalarIndex(i, j - 1)] + density_latest[ScalarIndex(i, j + 1)])) / (1 + 4 * a);
			}
		}
		SetDensityBoundary();
	}
}

void FluidSolver::VelocityDiffuse()
{
	double diff = 0.1;
	double a = dt * diff / (cell_width*cell_width);



	for (size_t k = 0; k < 20; k++) {
		for (size_t i = 1; i <= N - 1; i++)
		{
			for (size_t j = 1; j <= N; j++)
			{
				vx_latest[VxIndex(i, j)] = (vx_backup[VxIndex(i, j)] + a * (vx_latest[VxIndex(i - 1, j)] + vx_latest[VxIndex(i + 1, j)] +
					vx_latest[VxIndex(i, j - 1)] + vx_latest[VxIndex(i, j + 1)])) / (1 + 4 * a);

				vy_latest[VyIndex(j, i)] = (vy_backup[VyIndex(j, i)] + a * (vy_latest[VyIndex(j, i - 1)] + vy_latest[VyIndex(j, i + 1)] +
					vy_latest[VyIndex(j - 1, i)] + vy_latest[VyIndex(j + 1, i)])) / (1 + 4 * a);
			}
		}
		SetVelocityBoundary();
	}

}

void FluidSolver::Draw(QPainter& painter)
{
	double x, y, width_of_cell = 800.0 / (N + 2);
	int grey;
	//painter.setPen(QPen(QColor(0, 0, 0, 0), 1));
	for (size_t j = 0; j < N + 2; j++)
	{
		y = 0 + width_of_cell * j;
		for (size_t i = 0; i < N + 2; i++)
		{
			x = 0 + width_of_cell *i;
		 QPointF points[4] = {
				QPointF(x, y),
				QPointF(x+width_of_cell, y),
				QPointF(x+width_of_cell, y+width_of_cell),
				QPointF(x, y+width_of_cell)
			};
			grey = 255 * density_latest[ScalarIndex(i, j)];
			painter.setBrush(QColor(grey, grey, grey));
			painter.setPen(QPen(Qt::NoPen));
			painter.drawPolygon(points, 4);
		
		}

	}

}



void FluidSolver::SetInstantSource(double x1_, double y1_, double x2_, double y2_)
{
	int i1, j1;
	int i2, j2;
	double dx, dy, dl;

	dx = x2_ - x1_;
	dy = y2_ - y1_;
	dl = sqrt(dx*dx + dy * dy);

	i1 = floor((x1_ - 0.5*visual_width) / visual_width);
	j1 = floor((y1_ - 0.5*visual_width) / visual_width);
	i2 = floor((x2_ - 0.5*visual_width) / visual_width);
	j2 = floor((y2_ - 0.5*visual_width) / visual_width);
	int imax = i2 > i1 ? i2 : i1;
	int imin = i2 < i1 ? i2 : i1;
	int jmax = j2 > j1 ? j2 : j1;
	int jmin = j2 < j1 ? j2 : j1;
	for (size_t i = imin; i <= imax; i++)
	{
		for (size_t j = jmin; j <= jmax; j++)
		{
			density_latest[ScalarIndex(i, j)] = 1.0;
		}
	}

	i1 = floor((x1_ - visual_width) / visual_width);
	j1 = floor((y1_ - 0.5*visual_width) / visual_width);
	i2 = floor((x2_ - visual_width) / visual_width);
	j2 = floor((y2_ - 0.5*visual_width) / visual_width);
	if (i1 < 0) i1 = 0;
	if (i2 < 0) i2 = 0;
	if (j1 < 0) j1 = 0;
	if (j2 < 0) j2 = 0;
	imax = i2 > i1 ? i2 : i1;
	imin = i2 < i1 ? i2 : i1;
	jmax = j2 > j1 ? j2 : j1;
	jmin = j2 < j1 ? j2 : j1;
	for (size_t i = imin; i <= imax; i++)
	{
		for (size_t j = jmin; j <= jmax; j++)
		{
			vx_latest[VxIndex(i, j)] = 10.0 * dx / dl;
		}
	}

	j1 = floor((y1_ - visual_width) / visual_width);
	i1 = floor((x1_ - 0.5*visual_width) / visual_width);
	j2 = floor((y2_ - visual_width) / visual_width);
	i2 = floor((x2_ - 0.5*visual_width) / visual_width);
	if (i1 < 0) i1 = 0;
	if (i2 < 0) i2 = 0;
	if (j1 < 0) j1 = 0;
	if (j2 < 0) j2 = 0;
	imax = i2 > i1 ? i2 : i1;
	imin = i2 < i1 ? i2 : i1;
	jmax = j2 > j1 ? j2 : j1;
	jmin = j2 < j1 ? j2 : j1;
	for (size_t i = imin; i <= imax; i++)
	{
		for (size_t j = jmin; j <= jmax; j++)
		{
			vy_latest[VyIndex(i, j)] = 10.0 * dy / dl;
		}
	}


	SyncDensity();
	SyncV();


}



void FluidSolver::AddPermantSource(double x1_, double y1_, double x2_, double y2_)
{
	std::vector<double> permanent_source;
	permanent_source.push_back(x1_);
	permanent_source.push_back(y1_);
	permanent_source.push_back(x2_);
	permanent_source.push_back(y2_);
	permanent_sources.push_back(permanent_source);
}


void FluidSolver::SetPermantSource()
{
	if (permanent_sources.size()==0)
	{
		return;
	}
	int i1, j1;
	int i2, j2;
	double dx, dy, dl;
	for (size_t i = 0; i < permanent_sources.size(); i++)
	{
		dx = permanent_sources[i][2] - permanent_sources[i][0];
		dy = permanent_sources[i][3] - permanent_sources[i][1];
		dl = sqrt(dx*dx + dy * dy);
		i1 = floor((permanent_sources[i][0] - 0.5*visual_width) / visual_width);
		j1 = floor((permanent_sources[i][1] - 0.5*visual_width) / visual_width);
		i2 = floor((permanent_sources[i][2] - 0.5*visual_width) / visual_width);
		j2 = floor((permanent_sources[i][3] - 0.5*visual_width) / visual_width);

		int imax = i2 > i1 ? i2 : i1;
		int imin = i2 < i1 ? i2 : i1;
		int jmax = j2 > j1 ? j2 : j1;
		int jmin = j2 < j1 ? j2 : j1;
		for (size_t i = imin; i <= imax; i++)
		{
			for (size_t j = jmin; j <= jmax; j++)
			{
				density_latest[ScalarIndex(i, j)] = 1.0;
			}
		}

		i1 = floor((permanent_sources[i][0]-visual_width) / visual_width);
		j1 = floor((permanent_sources[i][1]-0.5*visual_width) / visual_width);
		i2 = floor((permanent_sources[i][2]-visual_width) / visual_width);
		j2 = floor((permanent_sources[i][3]-0.5*visual_width) / visual_width);
		if (i1 < 0) i1 = 0;
		if (i2 < 0) i2 = 0;
		if (j1 < 0) j1 = 0;
		if (j2 < 0) j2 = 0;
		
		imax = i2 > i1 ? i2 : i1;
		imin = i2 < i1 ? i2 : i1;
		jmax = j2 > j1 ? j2 : j1;
		jmin = j2 < j1 ? j2 : j1;
		for (size_t i = imin; i <= imax; i++)
		{
			for (size_t j = jmin; j <= jmax; j++)
			{
				vx_latest[VxIndex(i, j)] = 10 * dx / dl;
			}
		}



		j1 = floor((permanent_sources[i][0] - visual_width) / visual_width);
		i1 = floor((permanent_sources[i][1] - 0.5*visual_width) / visual_width);
		j2 = floor((permanent_sources[i][2] - visual_width) / visual_width);
		i2 = floor((permanent_sources[i][3] - 0.5*visual_width) / visual_width);
		if (i1 < 0) i1 = 0;
		if (i2 < 0) i2 = 0;
		if (j1 < 0) j1 = 0;
		if (j2 < 0) j2 = 0;

		imax = i2 > i1 ? i2 : i1;
		imin = i2 < i1 ? i2 : i1;
		jmax = j2 > j1 ? j2 : j1;
		jmin = j2 < j1 ? j2 : j1;
		for (size_t i = imin; i <= imax; i++)
		{
			for (size_t j = jmin; j <= jmax; j++)
			{
				vy_latest[VyIndex(i, j)] = 10 * dy / dl;
			}
		}
				




	}

}

void FluidSolver::ClearAll()
{
	for (size_t i = 0; i < N + 2; i++)
	{

		for (size_t j = 0; j < N + 2; j++)
		{


			{
				density_latest[ScalarIndex(i, j)] = 0.0;
			}


			//if (i <= NON_ZERO_REGION && abs(j - (N + 2) / 2.0) < NON_ZERO_REGION / 2.0)
			//{
			//	density_latest[ScalarIndex(i, j)] = 1;
			//}
			//else
			//{
			//	density_latest[ScalarIndex(i, j)] = 0;
			//}

			divV[ScalarIndex(i, j)] = 0.0;
			pressure_bar[ScalarIndex(i, j)] = 0.0;
		}
	}

	for (size_t i = 0; i < N + 1; i++)
	{
		for (size_t j = 0; j < N + 2; j++)
		{

			{
				vx_latest[VxIndex(i, j)] = 0.0;
			}



			{
				vy_latest[VyIndex(j, i)] = 0.0;
			}




		}

	}

	SyncV();
	SyncDensity();
	permanent_sources.clear();
}

void FluidSolver::ClearPermanetSource()
{
	permanent_sources.clear();

}