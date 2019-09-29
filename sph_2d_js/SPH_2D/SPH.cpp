#include <iostream>
#include <GL/glut.h>
#include <vector>

using namespace std;

#define MAX_PARTICLES 100

unsigned int scr_width = 640;
unsigned int scr_height = 640;

int view_width = 100;
int view_height = 100;

bool bStart = false;

struct Vector2D {
	Vector2D() {}
	Vector2D(float _x, float _y) : x(_x), y(_y) {}
	float x;
	float y;
};

struct Particle {
	Particle(float _x, float _y) : position(_x, _y), velocity(0.0f, 0.0f), density(0.0f), pressure(0.0f), Mass(1) {}
	Vector2D position;
	Vector2D velocity;
	float density;
	float pressure;
	int Mass;
};

vector<Particle> particles;

void InitParticle()
{
	for (float y = 20; y < 30; y++)
	{
		for (float x = -5; x < 5; x++)
		{
			if (particles.size() < MAX_PARTICLES)
				particles.push_back(Particle(x, y));
		}
	}
	cout << particles.size()<< " Paricles" << endl;
}

void Reshape(int w, int h){
	glViewport(0, 0, (GLsizei)w, (GLsizei)h);
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-view_width, view_width, -view_height, view_height, 0, 100);
	glutPostRedisplay();
}

void Display(void){
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(0.0f, 0.0f, 20.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);

	//render simulation
	glColor3f(1.0f, 1.0f, 1.0f);
	glPointSize(2.0f);
	glBegin(GL_POINTS);
	for (auto &p : particles)
	{
		glVertex2f(p.position.x, p.position.y);
		//cout << p.position.x << ", " << p.position.y << endl;
	}
	//glVertex2f(0.0f, 0.0f);
	glEnd();

	glutSwapBuffers();	
}

void InitGL()
{
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowSize(scr_width, scr_height);
	glutCreateWindow("SPH 2D");

	glEnable(GL_DEPTH_TEST);
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
}

void Dokeyboard(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 'r':
	case 'R':
		//Init simulation
		InitParticle();
		break;
	case 's':
	case 'S':
		//Start simulation
		bStart = true;
		break;
	default:
		break;
	}
	glutPostRedisplay();
}

void Update()
{
	if (bStart)
	{
		//simulation Step

	}

	glutPostRedisplay();
}

int main(int argc, char** argv){
	glutInit(&argc, argv);

	InitGL();

	glutDisplayFunc(Display);
	glutReshapeFunc(Reshape);
	glutIdleFunc(Update);
	glutKeyboardFunc(Dokeyboard);

	//Init simulation
	InitParticle();

	glutMainLoop();

	return 0;
}