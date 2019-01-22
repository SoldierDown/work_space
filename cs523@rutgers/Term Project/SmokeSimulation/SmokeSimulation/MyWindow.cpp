#include "MyWindow.h"
#include "qmessagebox.h"
#include<qwidget.h>


MyWindow::MyWindow(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);
	this->setFixedSize(800,800);
	Init();
}


MyWindow::~MyWindow(void)
{

}


void MyWindow::Init(void)
{
	myWidget = new MyWidget();
	CreateActions();
	CreateToolBar();
	CreateMenu();
	setCentralWidget(myWidget);  
}


void MyWindow::CreateActions(void)
{

	aboutMeAction = new QAction(tr("About Me..."), this);
	connect(aboutMeAction, &QAction::triggered, this, &MyWindow::AboutMe);

	aboutThisAction = new QAction(tr("About This..."), this);
	connect(aboutThisAction, &QAction::triggered, this, &MyWindow::AboutThis);



	startAction = new QAction(tr("&Start"), this);
	connect(startAction, &QAction::triggered, myWidget, &MyWidget::Start);

	suspendAction = new QAction(tr("&Suspend"), this);
	connect(suspendAction, &QAction::triggered, myWidget, &MyWidget::Suspend);

	nextAction = new QAction(tr("&Next"), this);
	connect(nextAction, &QAction::triggered, myWidget, &MyWidget::Next);	
	
	clearAllAction = new QAction(tr("&Clear All"), this);
	connect(clearAllAction, &QAction::triggered, myWidget, &MyWidget::ClearAll);	
	
	clearPermanentSourceAction = new QAction(tr("&Clear Permanent"), this);
	connect(clearPermanentSourceAction, &QAction::triggered, myWidget, &MyWidget::ClearPermanentSource);

	addPermanentSourceAction = new QAction(tr("&Permanet"), this);
	connect(addPermanentSourceAction, &QAction::triggered, myWidget, &MyWidget::AddPermanentSource);

	addInstantSourceAction = new QAction(tr("&Instant"), this);
	connect(addInstantSourceAction, &QAction::triggered, myWidget, &MyWidget::AddInstantSource);

	forbidSourceAction = new QAction(tr("&Forbid"), this);
	connect(forbidSourceAction, &QAction::triggered, myWidget, &MyWidget::Forbid);

}


void MyWindow::CreateToolBar(void)
{
	mainToolbar = addToolBar(tr("Tool"));
	mainToolbar->addAction(startAction);
	mainToolbar->addAction(suspendAction);
	mainToolbar->addAction(nextAction);
	mainToolbar->addAction(clearPermanentSourceAction);
	mainToolbar->addAction(clearAllAction);
	

}


void MyWindow::CreateMenu(void)
{
	addMenu = menuBar()->addMenu(tr("&Add"));
	addMenu->addAction(addPermanentSourceAction);
	addMenu->addAction(addInstantSourceAction);
	addMenu->addAction(forbidSourceAction);


	aboutMenu = menuBar()->addMenu(tr("&About"));
	aboutMenu->addAction(aboutMeAction);
	aboutMenu->addAction(aboutThisAction);
	

}


void MyWindow::AboutMe(void)
{
	QMessageBox::about(this, tr("About Me..."), tr("Welcome! \n" "Here is Hz Su. "));
}


void MyWindow::AboutThis(void)
{
	QMessageBox::about(this, tr("About This..."), tr("This is an interactive version of 2-D smoke simulation developed by Hz Su."));
}