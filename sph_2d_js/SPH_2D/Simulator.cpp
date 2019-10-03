#include "Simulator.h"

Simulator::Simulator() : timeStep(0.001), nframe(0)
{
	mySPH = new SPH();
}

Simulator::~Simulator()
{
	delete mySPH;
}

void Simulator::Init_Simulation()
{
	//SPH Initilize
	//mySPH->Init();
	mySPH->InitPouring();
}

void Simulator::Update()
{
	Vector2D gravity(0.0f, -9.8f);
	//timeStep = 0.00115;
	timeStep = 0.00315;
	nframe++;
	if(nframe % 10 == 0)
		mySPH->InitPouring();

	for(int i = 0; i < 4; i++)
	{
		//SPH simulation Steps
		mySPH->Update(timeStep, gravity);
	}
}

void Simulator::Render()
{
	//SPH display vertices
	mySPH->Draw();
}

