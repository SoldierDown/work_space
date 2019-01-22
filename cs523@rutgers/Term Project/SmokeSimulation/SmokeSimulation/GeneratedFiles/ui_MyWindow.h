/********************************************************************************
** Form generated from reading UI file 'MyWindow.ui'
**
** Created by: Qt User Interface Compiler version 5.9.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MYWINDOW_H
#define UI_MYWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_MyWindowClass
{
public:
    QMenuBar *menuBar;
    QToolBar *mainToolBar;
    QWidget *centralWidget;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *MyWindowClass)
    {
        if (MyWindowClass->objectName().isEmpty())
            MyWindowClass->setObjectName(QStringLiteral("MyWindowClass"));
        MyWindowClass->resize(600, 400);
        menuBar = new QMenuBar(MyWindowClass);
        menuBar->setObjectName(QStringLiteral("menuBar"));
        MyWindowClass->setMenuBar(menuBar);
        mainToolBar = new QToolBar(MyWindowClass);
        mainToolBar->setObjectName(QStringLiteral("mainToolBar"));
        MyWindowClass->addToolBar(mainToolBar);
        centralWidget = new QWidget(MyWindowClass);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        MyWindowClass->setCentralWidget(centralWidget);
        statusBar = new QStatusBar(MyWindowClass);
        statusBar->setObjectName(QStringLiteral("statusBar"));
        MyWindowClass->setStatusBar(statusBar);

        retranslateUi(MyWindowClass);

        QMetaObject::connectSlotsByName(MyWindowClass);
    } // setupUi

    void retranslateUi(QMainWindow *MyWindowClass)
    {
        MyWindowClass->setWindowTitle(QApplication::translate("MyWindowClass", "MyWindow", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class MyWindowClass: public Ui_MyWindowClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MYWINDOW_H
