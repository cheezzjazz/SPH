#include <iostream>
#include <vector>
#include <GL/glut.h>
#include "Simulator.h"

using namespace std;

unsigned int scr_width = 600;
unsigned int scr_height = 600;

double view_width = 20;
double view_height = 20;

bool bStart = false;

// for Mouse callback
int Mouse_PrevCoord[2] = {0,};
unsigned char Mouse_LeftEvent = 0;

Simulator mySimulator;

void InitGL()
{
	glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
	glutInitWindowPosition(400, 100);
	glutInitWindowSize(scr_width, scr_height);
	glutCreateWindow("SPH 2D");

	glEnable(GL_DEPTH_TEST);
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
}

void Display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(0.0f, 0.0f, 20.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);

	mySimulator.Render();

	glutSwapBuffers();
}

void Reshape(int _w, int _h)
{
	glViewport(0, 0, _w, _h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-view_width, view_width, -view_height, view_height, 0, 100);

	glutPostRedisplay();
}

void Update()
{
	if (bStart)
	{
		mySimulator.Update();
	}

	glutPostRedisplay();
}

void DoKeyboard(unsigned char key, int x, int y)
{
	switch (key)
	{
	case'r':
	case 'R':
		mySimulator.Init_Simulation();
		break;
	case 27://esc
	case 'q':
	case 'Q':
		exit(0);
		break;
	case 's':
	case 'S':
	case ' ': //space
		bStart = !bStart;
		break;
	default:
		break;
	}
}

void DoMouse(int mouse_event, int state, int x, int y)
{
	Mouse_PrevCoord[0] = x;
	Mouse_PrevCoord[1] = y;

	switch (mouse_event)
	{
	case GLUT_LEFT_BUTTON:
		Mouse_LeftEvent = (GLUT_DOWN == state) ? 1 : 0;
		break;
	default:
		break;
	}

	glutPostRedisplay();
}

void DoMotion(int x, int y)
{
	int Xdiff = x - Mouse_PrevCoord[0];
	int Ydiff = y - Mouse_PrevCoord[1];
	Mouse_PrevCoord[0] = x;
	Mouse_PrevCoord[1] = y;

	if (Mouse_LeftEvent)
	{
		//cloth -> add_force
	}

	glutPostRedisplay();
}

int main(int argc, char ** argv)
{
	glutInit(&argc, argv);
	InitGL();

	glutDisplayFunc(Display);
	glutReshapeFunc(Reshape);
	glutIdleFunc(Update);
	glutKeyboardFunc(DoKeyboard);
	glutMouseFunc(DoMouse);	//for mouse click
	glutMotionFunc(DoMotion);	//for mouse motion

	mySimulator.Init_Simulation();

	glutMainLoop();

	return 0;
}