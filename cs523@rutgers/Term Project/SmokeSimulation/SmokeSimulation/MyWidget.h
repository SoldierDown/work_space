#pragma once

#include <QWidget>
#include "ui_MyWidget.h"
#include<vector>
#include"FluidSolver.h"
#define N 200
#define RUNNING_TIME 100

using namespace std;


class MyWidget : public QWidget
{
	Q_OBJECT

public:
	MyWidget(QWidget *parent = Q_NULLPTR);
	~MyWidget(void);
	void	mousePressEvent(QMouseEvent *event);
	void	mouseMoveEvent(QMouseEvent *event);
	void	mouseReleaseEvent(QMouseEvent *event);
	void	paintEvent(QPaintEvent *);

	void	Start(void);
	void	Suspend(void);
	void	Next(void);
	void	AddPermanentSource(void);
	void	AddInstantSource(void);
	void	ClearAll(void);
	void	ClearPermanentSource(void);
	void	Forbid(void);


private:
	Ui::MyWidget ui;
	bool			rightClicked;
	bool			go_on = false;
	bool			first_time = true;
	bool			add_permanent_source = false;
	bool			add_instant_source = false;
	double			x1, y1, x2, y2;
	FluidSolver		fs;
};



