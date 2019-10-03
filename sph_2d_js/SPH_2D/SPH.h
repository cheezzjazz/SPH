#pragma once
#include "Particle.h"
#include <vector>

using namespace std;
#define MAX_PARTICLES 1600
#define GAS_CONSTANT 100
#define GRIDSIZE 40
#define HASHSIZE GRIDSIZE+2//103

class SPH
{
public:
	vector<Particle *> particles;
	float view_height = 20;	//main 이랑 연결해야함
	float view_width = 20;
	float pointR = 0.0f;// 1.0f;
	float PI = 3.1415926;
	int index = 0;
private:
	float rest_density = 4.5f;
	float u = 1.0f;	// viscosity constant
	float h = 1.0f;
	float surfaceCoef = 1.0f;
	enum kernels {
		poly6 = 0,
		spiky,
		viscosity
	};
public:
	SPH();
	~SPH();

	void ResetParticle();
	void Init();
	void InitPouring();
	void add_force(Vector2D ext_force); //확인필요
	void Update(float dt, Vector2D gravity);
	void Draw();
private:
	//poly6 kernel
	float cpoly6 = 4.0f / PI;
	float cpoly6G = -24.0f / PI;
	float cpoly6L = -48.0f / PI;
	float cspikyG = -30.0f / PI;
	float cviscosityL = 40 / PI;
	float Kernel(Vector2D rij, float norm, kernels kernel);
	Vector2D gradientKernel(Vector2D rij, float norm, kernels kernel);
	float laplacianKernel(Vector2D rij, float norm, kernels kernel);

	//float densKernel(Vector2D x, float h_);
	float calcSubpressure(float gas_constant, Particle *pi, Particle *pj, float rest_density);
	//spiky gradient kernel
	//float spiky = 15 / PI;
	
	Vector2D gradientKernel(Vector2D x, float q_);
	
	Vector2D gradientKernelPoly(Vector2D x, float norm_);
	
	float laplacianKernelPoly(Vector2D x, float norm_);
	//viscosity laplacian kernel
	
	float laplacianKernel(Vector2D x, float q_);
	Vector2D calcSubviscosity(Particle* pi, Particle* pj);

	void computeDensity();
	void computePressure();
	void computeViscosity();
	void computePresnVisco();
	void computeSurfacetension();
	void Integrate(float dt, Vector2D gravity);
private:
	vector<int>		hashGrid[HASHSIZE][HASHSIZE];
	//vecntor<Particle *> hashGrid[HASHSIZE][HASHSIZE];

	void HashTable();
	vector<int>GetNeighbor(float x, float y, float radius);
};