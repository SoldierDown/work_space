#include "MyWidget.h"
#include<qevent.h>
#include<qpainter.h>
#include<iostream>


//#include<rectypes.h>
using namespace std;

MyWidget::MyWidget(QWidget *parent)
	: QWidget(parent)
{ 
	ui.setupUi(this);
	rightClicked = false;
}

// avoid memory leak
MyWidget::~MyWidget(void)
{


}

void MyWidget::mousePressEvent(QMouseEvent *event)
{
	if (Qt::LeftButton == event->button())
	{
		x1 = event->pos().x();
		y1 = event->pos().y();
	}

}

void MyWidget::mouseMoveEvent(QMouseEvent *event)
{

}

//鼠标释放事件
void MyWidget::mouseReleaseEvent(QMouseEvent *event)
{
	x2 = event->pos().x();
	y2 = event->pos().y();
	if (add_permanent_source)
	{
		fs.AddPermantSource(x1, y1, x2, y2);
		fs.SetPermantSource();
	}
	if (add_instant_source)
	{
		fs.SetInstantSource(x1,y1,x2,y2);
	}
}


void MyWidget::paintEvent(QPaintEvent *event)
{

	QPainter painter(this);
	if (go_on)
	{
		fs.Draw(painter);
		first_time = false;
		update();
		Start();
	}
	else if (!go_on && ! first_time)
	{
		fs.Draw(painter);
		update();
	}
	else
	{
		update();
	}


}



void MyWidget::Start(void)
{

		go_on = true;
		//fs.Init();

		fs.Update();
}


void MyWidget::Next(void)
{
	fs.Update();
}


void MyWidget::Suspend()
{
	go_on = false;
}


void MyWidget::AddPermanentSource()
{
	add_permanent_source = true;
	add_instant_source = false;
}

void MyWidget::AddInstantSource()
{
	add_instant_source = true;
	add_permanent_source = false;
}


void MyWidget::Forbid()
{
	add_instant_source = false;
	add_permanent_source = false;
}



void MyWidget::ClearAll()
{
	fs.ClearAll();
}

void MyWidget::ClearPermanentSource()
{
	fs.ClearPermanetSource();
}