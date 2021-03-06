// Ball_Movement.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#define _CRTDBG_MAP_ALLOC  
#include <stdlib.h>  
#include <crtdbg.h>
#include <iostream>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>


#define GL_GLEXT_PROTOTYPES
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include<GL/glut.h>
#endif



#include <iostream>
#include<Eigen/Dense>
#include"RigidBody.h"
#include"Solv.h"





#ifdef _DEBUG
#define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
// Replace _NORMAL_BLOCK with _CLIENT_BLOCK if you want the
// allocations to be of _CLIENT_BLOCK type
#else
#define DBG_NEW new
#endif
void display();
void timer(int i_);
double frame;
void specialKeys(int key, int x, int y);
double rotate_y = 0;
double rotate_x = 0;
double vertex_data[TOTAL_RUNNING_TIME*NBODIES * 3 * NVertex] =
{

};



int time = 0;
int main(int argc, char* argv[])
{
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	//std::cout << sqrt(0.08) << std::endl;
	Solv solver;
	RigidBody bodies[NBODIES];
	std::cout << "Enter Frame:" << std::endl;
	std::cin >> frame;

	double ori_vertex_data[NBODIES * NVertex * 3];
	
	
		/*std::cout << "Enter Px, Py, Pz for " << i+1 << " separately:" << std::endl;
		std::cin >> px[i] >> py[i] >> pz[i];

		std::cout << "Enter Lx, Ly, Lz for " << i+1 << " separately:" << std::endl;
		std::cin >> lx[i] >> ly[i] >> lz[i];*/
		bodies[0] = RigidBody(0, -700, 0, 1, 1, 1);
		bodies[1] = RigidBody();
		for (size_t i = 0; i < NBODIES; i++)
		{
			for (size_t j = 0; j < NVertex * 3; j++)
			{
				ori_vertex_data[i*NVertex * 3 + j] = bodies[i].GetVerteices()[j];
			}
		}



	
	



	solver.RunSimulation(bodies,ori_vertex_data, vertex_data);
	
	for (size_t i = 0; i < TOTAL_RUNNING_TIME * NBODIES*NVertex * 3; i++)
	{
		//std::cout << vertex_data[i] << " ";
		if ((i + 1) % 3 == 0)
		{
			if ((i + 1) % (NVertex * 3) == 0)
			{/*
				std::cout << "Average: "
					<< " " << (vertex_data[i - 23] + vertex_data[i - 20] + vertex_data[i - 17] + vertex_data[i - 14] + vertex_data[i - 11] + vertex_data[i - 8] + vertex_data[i - 5] + vertex_data[i - 2]) / 8.0
					<< " " << (vertex_data[i - 22] + vertex_data[i - 19] + vertex_data[i - 16] + vertex_data[i - 13] + vertex_data[i - 10] + vertex_data[i - 7] + vertex_data[i - 4] + vertex_data[i - 1]) / 8.0
					<< " " << (+vertex_data[i - 21] + vertex_data[i - 18] + vertex_data[i - 15] + vertex_data[i - 12] + vertex_data[i - 9] + vertex_data[i - 6] + vertex_data[i - 3] + vertex_data[i]) / 8.0;
			*/
			}
			
			//std::cout << std::endl;
		}

	}
	//  Initialize GLUT and process user parameters
	glutInit(&argc, argv);

	//  Request double buffered true color window with Z-buffer
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	
	// Create window
	glutInitWindowSize(1024, 1024);
	glutCreateWindow("Awesome Cube");
	
	//  Enable Z-buffer depth test
	glEnable(GL_DEPTH_TEST);
	


	// Callback functions
	
	glutDisplayFunc(display);
	glutSpecialFunc(specialKeys);
	glutTimerFunc(1000.0 / frame, timer, 0);
	//  Pass control to GLUT for events
	glutMainLoop();

	//  Return to OS
	
	return 0;

	
	_CrtDumpMemoryLeaks(); 
}

void display() {

	//  Clear screen and Z-buffer
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Reset transformations

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();


	glMatrixMode(GL_MODELVIEW);

	glLoadIdentity();
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glRotatef(0, 1.0f, 0.0f, 0.0f);
	glRotatef(0, 0.0f, 1.0f, 0.0f);
	glRotatef(0, 0.0f, 0.0f, 1.0f);
	//glTranslated(1, 1, 1);




	
	

	// Rotate when user changes rotate_x and rotate_y
	glRotatef(rotate_x, 1.0, 0.0, 0.0);
	glRotatef(rotate_y, 0.0, 1.0, 0.0);

	// Other Transformations
	// glScalef( 2.0, 2.0, 0.0 );          // Not included
	for (size_t i = 0; i < NBODIES; i++)
	{
		if (i == 0)
		{
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		}
		
		else
		{
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		}

		//R side - FRONT 0473
		glBegin(GL_POLYGON);
		//glColor3f(1.0, 0.0, 0.0);
		glColor3f(1.0, 1.0, 0.0);
		glVertex3d(vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 0 * 3], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 0 * 3 + 1], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 0 * 3 + 2]);
		glVertex3d(vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 4 * 3], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 4 * 3 + 1], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 4 * 3 + 2]);
		glVertex3d(vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 7 * 3], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 7 * 3 + 1], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 7 * 3 + 2]);
		glVertex3d(vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 3 * 3], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 3 * 3 + 1], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 3 * 3 + 2]);
		glEnd();

		// G side - BACK 2651

		glBegin(GL_POLYGON);
		//glColor3f(0.0, 0.5, 0.5);
		glVertex3d(vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 2 * 3], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 2 * 3 + 1], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 2 * 3 + 2]);
		glVertex3d(vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 6 * 3], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 6 * 3 + 1], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 6 * 3 + 2]);
		glVertex3d(vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 5 * 3], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 5 * 3 + 1], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 5 * 3 + 2]);
		glVertex3d(vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 1 * 3], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 1 * 3 + 1], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 1 * 3 + 2]);
		glEnd();

		// B side - RIGHT 1540
		glBegin(GL_POLYGON);
		//glColor3f(0.0, 0.0, 1.0);
		glVertex3d(vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 1 * 3], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 1 * 3 + 1], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 1 * 3 + 2]);
		glVertex3d(vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 5 * 3], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 5 * 3 + 1], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 5 * 3 + 2]);
		glVertex3d(vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 4 * 3], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 4 * 3 + 1], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 4 * 3 + 2]);
		glVertex3d(vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 0 * 3], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 0 * 3 + 1], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 0 * 3 + 2]);
		glEnd();

		// Green side - LEFT 3762
		glBegin(GL_POLYGON);
		//glColor3f(1.0, 1.0, 0.0);
		glVertex3d(vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 3 * 3], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 3 * 3 + 1], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 3 * 3 + 2]);
		glVertex3d(vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 7 * 3], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 7 * 3 + 1], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 7 * 3 + 2]);
		glVertex3d(vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 6 * 3], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 6 * 3 + 1], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 6 * 3 + 2]);
		glVertex3d(vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 2 * 3], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 2 * 3 + 1], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 2 * 3 + 2]);
		glEnd();

		// Blue side - TOP 1032
		glBegin(GL_POLYGON);
		//glColor3f(1.0, 0.0, 1.0);
		glVertex3d(vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 1 * 3], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 1 * 3 + 1], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 1 * 3 + 2]);
		glVertex3d(vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 0 * 3], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 0 * 3 + 1], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 0 * 3 + 2]);
		glVertex3d(vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 3 * 3], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 3 * 3 + 1], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 3 * 3 + 2]);
		glVertex3d(vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 2 * 3], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 2 * 3 + 1], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 2 * 3 + 2]);
		glEnd();

		// Red side - BOTTOM 4567
		glBegin(GL_POLYGON);
		//glColor3f(0.0, 0.5, 0.5);
		glVertex3d(vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 4 * 3], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 4 * 3 + 1], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 4 * 3 + 2]);
		glVertex3d(vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 5 * 3], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 5 * 3 + 1], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 5 * 3 + 2]);
		glVertex3d(vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 6 * 3], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 6 * 3 + 1], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 6 * 3 + 2]);
		glVertex3d(vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 7 * 3], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 7 * 3 + 1], vertex_data[time*NBODIES*NVertex * 3 + i * NVertex * 3 + 7 * 3 + 2]);
		glEnd();
	}
	
	
	glFlush();
	glutSwapBuffers();

}

void specialKeys(int key, int x, int y) {

	//  Right arrow - increase rotation by 5 degree
	if (key == GLUT_KEY_RIGHT)
		rotate_y += 5;

	//  Left arrow - decrease rotation by 5 degree
	else if (key == GLUT_KEY_LEFT)
		rotate_y -= 5;

	else if (key == GLUT_KEY_UP)
		rotate_x += 5;

	else if (key == GLUT_KEY_DOWN)
		rotate_x -= 5;

	//  Request display update
	glutPostRedisplay();

}

void timer(int i_)
{
	time += 1;
	glutPostRedisplay();
	glutTimerFunc(1000.0 / frame, timer, 0);
}
