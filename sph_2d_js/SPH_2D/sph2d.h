#ifndef  SPH_H
#define  SPH_H



#include <iostream>
#include <GL\glut.h>

using namespace std;

class sph2d
{
public:
	sph2d() {}
	//sph2d
	void display();
	void step();
	void Init();
private:
	void calculateDensity();
};

#endif // ! SPH_H