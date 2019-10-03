#pragma once
#include "Particle.h"
#include <vector>

using namespace std;
#define MAX_PARTICLES 2601
#define GAS_CONSTANT 80
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
	float surfaceCoef = 20.0f;
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

	float calcSubpressure(float gas_constant, Particle *pi, Particle *pj, float rest_density);
	Vector2D calcSubviscosity(Particle* pi, Particle* pj);

	void computeDensity();
	void computePressure();
	void computeViscosity();
	void computePresnVisco();
	void computeSurfacetension();
	void Integrate(float dt, Vector2D gravity);
private:
	//vector<int>		hashGrid[HASHSIZE][HASHSIZE];
	vector<Particle *> hashGrid[HASHSIZE][HASHSIZE];

	void HashTable();
	//vector<int>GetNeighbor(float x, float y, float radius);
	vector<Particle*> GetNeighbor(float x, float y, float radius);
	vector<Particle*> GetNeighbor(int gridx, int gridy, float radius, vector<Particle*>& mine);
	vector<Particle*> GetNeighbor(int gridx, int gridy, float radius);
	vector<Particle*> GetMine(int gridx, int gridy, float radius);
};