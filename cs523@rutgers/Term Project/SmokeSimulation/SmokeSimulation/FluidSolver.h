#pragma once
#include<qpainter.h>
#include<vector>
#include <QPainter>
#include <QPoint>
#define N 300
#define NON_ZERO_REGION 20
#include<math.h>
#include<vector>
#include<iostream>
#include<omp.h>
class FluidSolver
{
public:
	FluidSolver();
	~FluidSolver();


	void Init();
	int ScalarIndex(int i_, int j_);
	int VxIndex(int i_, int j_);
	int VyIndex(int i_, int j_);


	void DensityAdvect();
	void VelocityAdvect();

	void DensityDiffuse();

	void VelocityDiffuse();

	void Project();

	void SyncV();

	void SyncDensity();

	void SetDensityBoundary();
	void SetPressureBoundary();
	void SetDivVBoundary();

	void SetInstantSource(double x1_,double y1_,double x2_,double y2_);
	void AddPermantSource(double x1_, double y1_, double x2_, double y2_);
	void SetPermantSource();

	void ClearAll();
	void ClearPermanetSource();

	void Draw(QPainter& painter);

	void Update();
	void SetVelocityBoundary();

	void Save(double density_date_[]);

	double Distance(double x_, double y_);
	double Cos(double x_, double y_);
	double Sin(double x_, double y_);


	


	double cell_width = 1;
	double visual_width = 800.0 / (N + 2);
	double dt = 0.05;

	std::vector<std::vector<double>> permanent_sources;

	double density_latest[(N + 2)*(N + 2)];
	double density_backup[(N + 2)*(N + 2)];

	double pressure_bar[(N + 2)*(N + 2)];
	double divV[(N + 2)*(N + 2)];


	double vx_latest[(N + 1)*(N + 2)];
	double vx_backup[(N + 1)*(N + 2)];

	double vy_latest[(N + 2)*(N + 1)];
	double vy_backup[(N + 2)*(N + 1)];


	double total;
};

