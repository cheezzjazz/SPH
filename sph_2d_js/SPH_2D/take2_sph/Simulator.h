#pragma once
#include <iostream>
#include <vector>
#include <GL/glut.h>
#include "SPH.h"

using namespace std;

class Simulator
{
public:
	Simulator();
	~Simulator();

	void Init_Simulation();
	void Update();
	void Render();
private:
	SPH *mySPH;
	float timeStep;

	int nframe;
};