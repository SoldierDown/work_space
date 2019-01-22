#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_MyWindow.h"
#include"MyWidget.h"


class MyWindow : public QMainWindow
{
	Q_OBJECT

public:
	MyWindow(QWidget *parent = 0);
	~MyWindow(void);

private:
	Ui::MyWindowClass ui;
	MyWidget*	myWidget;
	QAction*	aboutMeAction;
	QAction*	aboutThisAction;
	QAction*	startAction;
	QAction*	suspendAction;
	QAction*	nextAction;
	QAction*	clearAllAction;
	QAction*	clearPermanentSourceAction;
	QAction*	addInstantSourceAction;
	QAction*	addPermanentSourceAction;
	QAction*	forbidSourceAction;
	QMenu*		addMenu;
	QMenu*		aboutMenu;
	QToolBar*	mainToolbar;








	void	Init(void);
	void	CreateActions(void);
	void	CreateToolBar(void);
	void	CreateMenu(void);
	void	AboutMe(void);
	void	AboutThis(void);
};
