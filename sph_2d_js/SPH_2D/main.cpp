#include <iostream>
#include <GL/glut.h>
#include <vector>

using namespace std;

#define MAX_PARTICLES 400
#define GAS_CONSTANT 100
//#define POLY6 4.0f/pi

unsigned int scr_width = 400;
unsigned int scr_height = 400;

float pi = 3.1415926;
float pointR = 1.0f;
float dt = 0.00315;

int view_width = 20;
int view_height = 20;

bool bStart = true;

struct Vector2D {
	Vector2D() {}
	Vector2D(float _x, float _y) : x(_x), y(_y) {}
	float normsqr()
	{
		return x*x + y*y;
	}
	float norm() {
		return sqrt(x*x + y*y);
	}

	float x;
	float y;
};

Vector2D operator+(Vector2D v1, Vector2D v2)
{
	Vector2D result;
	result.x = v2.x + v1.x;
	result.y = v2.y + v1.y;

	return result;
}

Vector2D operator-(Vector2D v1, Vector2D v2)
{
	Vector2D result;
	result.x = v1.x - v2.x;
	result.y = v1.y - v2.y;

	return result;
}

Vector2D operator-=(Vector2D v1, Vector2D v2)
{
	Vector2D result;
	result.x = v1.x - v2.x;
	result.y = v1.y - v2.y;

	return result;
}

Vector2D operator+=(Vector2D v1, Vector2D v2)
{
	Vector2D result;
	result.x = v1.x + v2.x;
	result.y = v1.y + v2.y;

	return result;
}

bool operator==(Vector2D v1, Vector2D v2)
{
	bool result = false;

	if (v1.x == v2.x && v1.y == v2.y)
		result = true;

	return result;
}

Vector2D operator*(float f, Vector2D v)
{
	Vector2D result;
	result.x = f*v.x;
	result.y = f*v.y;

	return result;
}

Vector2D operator*(Vector2D v, float f)	//교환법칙 왜안되지
{
	Vector2D result;
	result.x = f*v.x;
	result.y = f*v.y;

	return result;
}

Vector2D operator*=(Vector2D v, float f)
{
	Vector2D result;
	result.x = f*v.x;
	result.y = f*v.y;

	return result;
}

Vector2D operator/(Vector2D v, float f)
{
	Vector2D result;
	result.x = v.x/f;
	result.y = v.y/f;

	return result;
}

struct Particle {
	Particle(float _x, float _y) : position(_x, _y), velocity(0.0f, 0.0f), fpressure(0.0f, 0.0f), 
		fviscosity(0.0f, 0.0f), accel(0.0f, 0.0f), density(0.0f), pressure(0.0f), Mass(1), idx(0) {}
	Particle(float _x, float _y, int _idx) : position(_x, _y), velocity(0.0f, 0.0f), fpressure(0.0f, 0.0f),
		fviscosity(0.0f, 0.0f), accel(0.0f, 0.0f), density(0.0f), pressure(0.0f), Mass(1), idx(_idx) {}
	Vector2D position;
	Vector2D velocity;
	Vector2D fpressure;
	Vector2D fviscosity;
	Vector2D accel;

	int idx;
	float density;
	float pressure;
	int Mass;
};

bool operator==(Particle p1, Particle p2)
{
	bool result = false;

	if (p1.position.x == p2.position.x && p1.position.y == p2.position.y)
		result = true;

	return result;
}

vector<Particle> particles;

void InitParticle()
{
	while(!particles.empty())
	{
		particles.pop_back();
	}

	int idx = 0;
	for (float y = -(view_height/2); y < (view_height / 2); y+=0.5)
	{
		//for (float x = -(view_width/2); x < (view_width / 2); x+=0.5)
		for (float x = -(view_width-pointR ); x < (view_width / 6); x += 0.5)
		{
			if (particles.size() < MAX_PARTICLES)
				particles.push_back(Particle(x, y, idx++));
				//particles.push_back(Particle(x, y));
		}
	}
	cout << particles.size() << " Paricles" << endl;
}

//poly6 kernel
float poly6 = 4.0f / pi ;
float densKernel(Vector2D x, float h)
{
	/*float result = 0.0f;
	float poly6 = 4.0f / pi / (h*h*h*h*h*h*h*h);
	result = poly6* (h*h-x.normsqr())*(h*h - x.normsqr())*(h*h - x.normsqr());*/
	float result = 0.0f;
	float q = x.norm() / h;
	
	//result = poly6 * (1 - q*q) * (1 - q * q) * (1 - q * q);
	result = (1 - q * q) * (1 - q * q) * (1 - q * q);
	return result;
}

float Kernel(Vector2D x, float h)
{
	float q;
	float result;
	float constant;
	constant = 40.0f / (h * h * 7.0f * pi);

	q = x.norm() / h;


	if (q >= 0.0f && q <= 0.5)
	{
		result = constant * (6.0f * pow(q, 3) - 6.0f * (q * q) + 1.0f);

		return result;
	}
	else if (q >= 0.5 && q <= 1.0f)
	{
		result = constant * (2.0f * pow(1 - q, 3));
		return result;
	}
	else if (q > 1.0f)
	{
		return 0;
	}


}

Vector2D GradientKernel(Vector2D x, float h)
{
	float q;
	Vector2D gradientq;
	float constant;
	constant = 6.0f * 40.0f / (h * h * 7.0f * pi);

	q = x.norm() / h;
	gradientq = x / (x.norm() * h);
	
	if (q > 0 && q <= 0.5)
	{
		gradientq = gradientq * constant * (3.0f * q * q - 2.0f * q);

		return gradientq;
	}
	else if (q >= 0.5 && q <= 1.0f)
	{
		gradientq = gradientq * constant * (-(1.0f - q) * (1.0f - q));

		return gradientq;
	}
	else if (q > 1.0f || q == 0)
	{
		gradientq = Vector2D(0.0, 0.0);

		return gradientq;
	}
}

float LapKernel(Vector2D x, float h)
{
	float q;
	float laplacianq;
	float result;
	float constant;
	constant = 6.0f * 40.0f / (h * h * 7.0f * pi);
	q = x.norm() / h;
	laplacianq = 2.0f / (x.norm() * h + 0.0000000001);

	if (q >= 0.0 && q <= 0.5)
	{

		result = constant * ((6.0f * q - 2.0f) / (h * h) + (3.0f * q * q - 2.0f * q) * laplacianq);

		return result;
	}
	else if (q >= 0.5f && q <= 1.0f)
	{

		result = constant * ((2.0f * (1.0f - q)) / (h * h) + (-(1.0f - q) * (1.0f - q)) * laplacianq);
		return result;
	}

	else if (q > 1)
	{
		return 0;
	}

}

//spiky gradient kernel
//float spiky = 15 / pi;
float spiky = -30 / pi ;
Vector2D gradientKernel(Vector2D x, float h)
{
	float q = x.norm() / h;
	//float spiky = -30 / pi / (h*h*h*h);
	//float temp = spiky*(1 - q)*(1 - q) / q;
	float temp = (1 - q) * (1 - q) / q;
	Vector2D result(temp*x.x, temp*x.y);
	return result;
	//return ((1.0f - x.norm())*(1.0f - x.norm()) / x.norm())*x;
}

//viscosity laplacian kernel
float viscosity = 40 / pi;
float laplacianKernel(Vector2D x, float h)
{
	//float viscosity = 40 / pi / (h*h*h*h);

	//return viscosity * (1 - x.norm());
	return (1 - x.norm());
}

void computeDensity()
{
	for (auto &pi : particles)
	{
		pi.density = 0.0f;
		for (auto &pj : particles)
		{
			if (pi.idx == pj.idx)
				continue;

			Vector2D rij = pi.position - pj.position;
			float h = 1.0f;
			float q = rij.norm() / h;
			if (0<= q && q <= 1.0f)								
			{
				pi.density = pi.density + pj.Mass*densKernel(rij, h);//JS
				//pi.density = pi.density + pj.Mass * Kernel(rij, h);//QM
			}
		}
		pi.density = pi.density * poly6;//JS
		//cout << "dens " << pi.density << endl;
	}
}

float calcSubpressure(float gas_constant, Particle pi, Particle pj, float rest_density)
{
	float pres_i = gas_constant*(pi.density - rest_density);
	float pres_j = gas_constant*(pj.density - rest_density);
	return (pres_i + pres_j) / (2.0f*pj.density);
	//return gas_constant*(pi.density + pj.density - 2.0*rest_density)/pj.density;
}

Vector2D calcSubviscosity(Particle pi, Particle pj)
{
	if (pj.density == 0.0f)
		cout << "den = 0!" << endl;
	return (pj.velocity - pi.velocity) / pj.density;
}

void computePressure()
{
	for (auto &pi : particles)
	{
		pi.fpressure = Vector2D(0.0f, 0.0f);
		for (auto &pj : particles)
		{
			if (pi.idx == pj.idx)
				continue;

			Vector2D rij = pi.position - pj.position;
			float h = 1.0f;
			float q = rij.norm() / h;
			if (0 < q && q <= 1.0f)
			{
				float rest_density = 1.0f;
				
				pi.fpressure = pi.fpressure + pj.Mass*calcSubpressure(GAS_CONSTANT, pi, pj, rest_density)*gradientKernel(rij, h);//JS
				
				//pi.fpressure = pi.fpressure + pj.Mass * calcSubpressure(GAS_CONSTANT, pi, pj, rest_density) * GradientKernel(rij, h);//QM
			}
		}
		pi.fpressure = -1.0f * pi.fpressure * spiky;//JS
		//pi.fpressure = -1.0f * pi.fpressure;//QM
		//cout <<"pressure"<< pi.fpressure.x<<"," <<pi.fpressure.y<< endl;
	}
}

void computeViscosity()
{
	float u = 10.0f;	// viscosity constant
	for (auto &pi : particles)
	{
		pi.fviscosity = Vector2D(0.0f, 0.0f);
		for (auto &pj : particles)
		{
			if (pi.idx == pj.idx)
				continue;

			Vector2D rij = pi.position - pj.position;
			float h = 1.0f;
			float q = rij.norm() / h;
			if (0 <= q && q <= 1.0f)
			{
				pi.fviscosity = pi.fviscosity + pj.Mass*laplacianKernel(rij, h)*calcSubviscosity(pi, pj);//JS
				//pi.fviscosity = pi.fviscosity + u*pj.Mass * LapKernel(rij, h) * calcSubviscosity(pi, pj);//QM
			}
		}
		pi.fviscosity = pi.fviscosity * u * viscosity;//JS
		//cout << "fviscosity" << pi.fviscosity.x << "," << pi.fviscosity.y << endl;
	}


}

void Integrate()
{
	for (auto &p : particles)
	{
		Vector2D fgrav = p.Mass*Vector2D(0.0f, -9.8f);
		//Vector2D fgrav = p.Mass * Vector2D(0.0f, -1.9f);
		p.accel = (p.fpressure + p.fviscosity)/p.density + fgrav;///p.density;
		//cout << "accel" << p.accel.x << "," << p.accel.y << endl;
		p.velocity = p.velocity + p.accel*dt;
		//cout << "velocity" << p.velocity.x << "," << p.velocity.y << endl;
		p.position = p.position + p.velocity*dt;
		//cout << "position" << p.position.x << "," << p.position.y << endl;

		float restitution = 0.5f;//0.3f;
		float damping = 0.5f;
		if (p.position.x - pointR < -view_width)
		{
			Vector2D n = Vector2D(1.0f, 0.0f);
			float length = abs((p.velocity.x)*n.x + (p.velocity.y)*n.y);//JS
			//float length = (-p.velocity.x) * n.x + (-p.velocity.y) * n.y;//QM
			Vector2D n1 = length*n;
			p.velocity = p.velocity + (1.0f + restitution)*n1;//JS
			//p.velocity.x = 0.7 * n1.x + p.velocity.x;//QM

			//p.velocity.x = p.velocity.x*damping;
			p.position.x = -view_width + pointR;
			////cout << "!" << endl;
		}
		if (p.position.x + pointR > view_width)
		{
			Vector2D n = Vector2D(-1.0f, 0.0f);
			float length = abs(p.velocity.x*n.x + p.velocity.y*n.y);//JS
			//float length = (-p.velocity.x) * n.x + (-p.velocity.y) * n.y;//QM
			Vector2D n1 = length*n;
			p.velocity = p.velocity + (1.0f + restitution)*n1;//JS
			/*p.velocity.x = p.velocity.x*damping;*/
			//p.velocity.x = n1.x + p.velocity.x;//QM
			p.position.x = view_width - pointR;
			//cout << "!" << endl;
		}
		if (p.position.y - pointR < -view_height)
		{
			Vector2D n = Vector2D(0.0f, 1.0f);
			//cout << p.velocity.norm();
			float length = abs(p.velocity.x*n.x + p.velocity.y*n.y);//JS
			//float length = (-p.velocity.x) * n.x + (-p.velocity.y) * n.y;//QM
			//cout << " l:" << length << endl;
			Vector2D n1 = length*n;
			p.velocity = p.velocity + (1.0f + restitution)*n1;//JS
			//p.velocity.y = p.velocity.y*damping;
			//p.velocity.y = n1.y + p.velocity.y;//QM
			p.position.y = -view_height + pointR;
			//cout << "!" << endl;
		}
		if (p.position.y + pointR > view_height)
		{
			Vector2D n = Vector2D(0.0f, -1.0f);
			float length = abs(p.velocity.x*n.x + p.velocity.y*n.y);//JS
			//float length = (-p.velocity.x) * n.x + (-p.velocity.y) * n.y;//QM
			Vector2D n1 = length*n;
			p.velocity = p.velocity + (1.0f + restitution)*n1;//JS
			//p.velocity.y = p.velocity.y*damping;
			//p.velocity.y = n1.y + p.velocity.y;//QM
			p.position.y = view_height - pointR;
			//cout << "!" << endl;
		}
	}
}

void Reshape(int w, int h) {
	glViewport(0, 0, (GLsizei)w, (GLsizei)h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-view_width, view_width, -view_height, view_height, 0, 100);
	glutPostRedisplay();
}

void Display(void) {
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
		computeDensity();
		computePressure();
		computeViscosity();
		Integrate();
	}

	glutPostRedisplay();
}

int main(int argc, char** argv) {
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