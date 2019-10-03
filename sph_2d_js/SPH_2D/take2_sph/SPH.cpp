#include "SPH.h"
#include <time.h>
SPH::SPH()
{

}
SPH::~SPH()
{
	while (!particles.empty())
	{
		particles.pop_back();
	}
	particles.clear();
}

void SPH::ResetParticle()
{
	if (index >= MAX_PARTICLES)
	{
		index = 0;
		while (!particles.empty())
		{
			particles.pop_back();
		}
	}
}

void SPH::Init()
{
	ResetParticle();

	int idx = 0;
	//for (float y = -(view_height / 4); y < (view_height / 2); y += 0.4)
	for (float y = view_height - pointR; y > (view_height / 2); y -= 0.4)
	{
		for (float x = -(view_width / 4); x < (view_width / 2); x += 0.4)
			//for (float x = -(view_width - pointR); x < (view_width / 6); x += 0.5)
		{
			if (particles.size() < MAX_PARTICLES)
			{
				Particle *p = new Particle(x, y, idx++, Particle::fluid);
				particles.push_back(p);
			}
		}
	}
	cout << "SPH" << particles.size() << " Paricles" << endl;
}

void SPH::InitPouring()
{
	//ResetParticle();

	/*for (float y = view_height / 2; y > view_height / 2 - 1.0f; y -= 0.4)
	{
		for (float x = -1.0f; x < 1.0f; x += 0.4)
		{
			if (particles.size() < MAX_PARTICLES)
			{
				Particle *p = new Particle(x, y, index++, Particle::fluid);
				p->velocity.y = -10.0;
				particles.push_back(p);
			}
		}
	}*/
	//cout << "SPH" << particles.size() << " Paricles" << endl;
	for (float y = 10.0f; y < 13.0f; y += 0.4)
	{
		for (float x = -view_width; x < -view_width+3.0f; x += 0.4)
		{
			if (particles.size() < MAX_PARTICLES)
			{
				Particle *p = new Particle(x, y, index++, Particle::fluid);
				p->velocity.x = 20.0f;
				particles.push_back(p);
			}
		}
	}
}

//확인필요
void SPH::add_force(Vector2D ext_force)
{
	for (auto &p : particles)
	{
		p->add_force(ext_force);
	}
}

void SPH::Update(float dt, Vector2D gravity)
{
	HashTable();
	clock_t start = clock();
	computeDensity();
	computePresnVisco();
	//computePressure();
	//computeViscosity();
	computeSurfacetension();
	clock_t end = clock();
	float time = (float)(end - start) ;
	//if (time != 0.0f)
	printf("%0.5f\n", time);	//cout << "Time : " << time<<endl;
	//system("cls");
	Integrate(dt, gravity);
}

void SPH::Draw()
{
	for (auto &p : particles)
		p->draw();
}

void SPH::computeDensity()
{
	////float tmpDens = 0.0f;
	//for (auto &pi : particles)
	//{
	//	pi->density = 0.0f;
	//	vector<int> neighbors = GetNeighbor(pi->position.x, pi->position.y, h);//using hashgrid
	//	for (auto &pj : particles)
	//		//for(auto &j : neighbors)//using hashgrid
	//	{
	//		//Particle &pj = particles.at(j);//using hashgrid
	//		//cout << pi.idx << "," << pj.idx << endl;
	//		//if (pi.idx == pj.idx)
	//		//	continue;

	//		Vector2D rij = pi->position - pj->position;

	//		float norm = rij.norm();
	//		float q = norm / h;
	//		if (0 <= q && q < 1.0f)//q <= 1.0f)
	//		{
	//			pi->density = pi->density + pj->mass*Kernel(rij, norm, poly6);//JS
	//		}
	//	}
	//	pi->density = pi->density * cpoly6;//JS
	//	//tmpDens += pi->density;

	//}
	////tmpDens /= particles.size();
	////printf("%f\n", tmpDens);

	//float tmpDens = 0.0f;

	for (int i = 1; i < HASHSIZE - 1; i++)
	{
		for (int j = 1; j < HASHSIZE - 1; j++)
		{
			vector<Particle*> pis;
			vector<Particle*> pjs = GetNeighbor(i, j, h, pis);
			for (auto &pi : pis)
			{
				pi->density = 0.0f;
				for (auto &pj : pjs)
				{
					Vector2D rij = pi->position - pj->position;

					float norm = rij.norm();
					float q = norm / h;
					if (0 <= q && q < 1.0f)//q <= 1.0f)
					{
						pi->density = pi->density + pj->mass*Kernel(rij, norm, poly6);//JS
					}

				}
				pi->density = pi->density * cpoly6;//JS
				//tmpDens += pi->density;
			}
		}
	}

	//tmpDens /= particles.size();
	//printf("%f\n", tmpDens);
}

float SPH::Kernel(Vector2D rij, float norm, kernels kernel)
{
	float result = 0.0f;
	float temp = 0.f;
	switch (kernel)
	{
	case kernels::poly6:
		temp = 1 - norm*norm;
		result = temp*temp*temp;
		break;
	case kernels::spiky:
		break;
	case kernels::viscosity:
		break;
	default:
		break;
	}
	return result;
}

Vector2D SPH::gradientKernel(Vector2D rij, float norm, kernels kernel)
{
	Vector2D result = Vector2D(0.0f, 0.0f);
	float temp = 0.f;
	switch (kernel)
	{
	case kernels::poly6:
		temp = 1.0f - norm*norm;
		temp *= temp;
		result = Vector2D(temp*rij.x, temp*rij.y);
		break;
	case kernels::spiky://norm = q
		temp = (1 - norm) * (1 - norm) / norm;
		result = Vector2D(temp*rij.x, temp*rij.y);
		break;
	case kernels::viscosity:
		break;
	default:
		break;
	}
	return result;
}

float SPH::laplacianKernel(Vector2D rij, float norm, kernels kernel)
{
	float result = 0.0f;
	float temp = 0.f;
	switch (kernel)
	{
	case kernels::poly6:
		temp = norm*norm;
		result = (1 - temp)*(1 - 3 * temp);
		break;
	case kernels::spiky:
		break;
	case kernels::viscosity: // norm = q
		result = (1 - norm);
		break;
	default:
		break;
	}
	return result;
}

float SPH::calcSubpressure(float gas_constant, Particle *pi, Particle *pj, float rest_density)
{
	float pres_i = gas_constant*(pi->density - rest_density);
	float pres_j = gas_constant*(pj->density - rest_density);
	return (pres_i + pres_j) / (2.0f*pj->density);
	//return gas_constant*(pi.density + pj.density - 2.0*rest_density)/pj.density;
}

void SPH::computePressure()
{
	//for (auto &pi : particles)
	//{
	//	pi->fpressure = Vector2D(0.0f, 0.0f);
	//	vector<int> neighbors = GetNeighbor(pi->position.x, pi->position.y, h);//using hashgrid
	//																		   //for (auto &pj : particles)
	//	for (auto &j : neighbors)//using hashgrid
	//	{
	//		Particle *&pj = particles.at(j);
	//		if (pi->idx == pj->idx)
	//			continue;

	//		Vector2D rij = pi->position - pj->position;
	//		float q = rij.norm() / h;
	//		if (0 < q && q <= 1.0f)
	//		{

	//			pi->fpressure = pi->fpressure + gradientKernel(rij, q, spiky)*pj->mass*calcSubpressure(GAS_CONSTANT, pi, pj, rest_density);//JS
	//			//pi.fpressure = pi.fpressure + pj.Mass * calcSubpressure(GAS_CONSTANT, pi, pj, rest_density) * GradientKernel(rij, h);//QM
	//		}
	//	}
	//	pi->fpressure = pi->fpressure * -1.0f * cspikyG;//JS
	//}
	for (int i = 1; i < HASHSIZE - 1; i++)
	{
		for (int j = 1; j < HASHSIZE - 1; j++)
		{
			vector<Particle*> pis;
			vector<Particle*> pjs = GetNeighbor(i, j, h, pis);
			for (auto &pi : pis)
			{
				pi->fpressure = Vector2D(0.0f, 0.0f);
				for (auto &pj : pjs)
				{
					if (pi->idx == pj->idx)
						continue;

					Vector2D rij = pi->position - pj->position;
					float q = rij.norm() / h;
					if (0 < q && q <= 1.0f)
					{
						pi->fpressure = pi->fpressure + gradientKernel(rij, q, spiky)*pj->mass*calcSubpressure(GAS_CONSTANT, pi, pj, rest_density);//JS
					}
				}
				pi->fpressure = pi->fpressure * -1.0f * cspikyG;//JS
			}
		}
	}
}

Vector2D SPH::calcSubviscosity(Particle* pi, Particle* pj)
{
	if (pj->density == 0.0f)
	{
		cout << "den = 0!" << endl;
		pj->density = 0.0001f;
	}

	return (pj->velocity - pi->velocity) / pj->density;
}

void SPH::computeViscosity()
{
	//for (auto &pi : particles)
	//{
	//	pi->fviscosity = Vector2D(0.0f, 0.0f);
	//	vector<int> neighbors = GetNeighbor(pi->position.x, pi->position.y, h);//using hashgrid
	//																		   //for (auto &pj : particles)
	//	for (auto &j : neighbors)//using hashgrid
	//	{
	//		Particle *&pj = particles.at(j);
	//		if (pi->idx == pj->idx)
	//			continue;

	//		Vector2D rij = pi->position - pj->position;
	//		float q = rij.norm() / h;
	//		if (0 <= q && q <= 1.0f)
	//		{
	//			pi->fviscosity = pi->fviscosity + calcSubviscosity(pi, pj)*pj->mass*laplacianKernel(rij, q, viscosity);//JS
	//		}
	//	}
	//	pi->fviscosity = pi->fviscosity * u * cviscosityL;//JS
	//}
	for (int i = 1; i < HASHSIZE - 1; i++)
	{
		for (int j = 1; j < HASHSIZE - 1; j++)
		{
			vector<Particle*> pis;
			vector<Particle*> pjs = GetNeighbor(i, j, h, pis);
			for (auto &pi : pis)
			{
				pi->fviscosity = Vector2D(0.0f, 0.0f);
				for (auto &pj : pjs)
				{
					if (pi->idx == pj->idx)
						continue;

					Vector2D rij = pi->position - pj->position;
					float q = rij.norm() / h;
					if (0 < q && q <= 1.0f)
					{
						pi->fviscosity = pi->fviscosity + calcSubviscosity(pi, pj)*pj->mass*laplacianKernel(rij, q, viscosity);//JS
					}
				}
				pi->fviscosity = pi->fviscosity * u * cviscosityL;//JS
			}
		}
	}
}

void SPH::computePresnVisco()
{
	//for (auto &pi : particles)
	//{
	//	pi->fpressure = Vector2D(0.0f, 0.0f);
	//	pi->fviscosity = Vector2D(0.0f, 0.0f);
	//	vector<int> neighbors = GetNeighbor(pi->position.x, pi->position.y, h);//using hashgrid
	//	for (auto &j : neighbors)//using hashgrid
	//	{
	//		Particle *&pj = particles.at(j);
	//		if (pi->idx == pj->idx)
	//			continue;

	//		Vector2D rij = pi->position - pj->position;
	//		float q = rij.norm() / h;
	//		if (0 < q && q <= 1.0f)
	//		{
	//			pi->fpressure = pi->fpressure + gradientKernel(rij, q, spiky)*pj->mass*calcSubpressure(GAS_CONSTANT, pi, pj, rest_density);//JS
	//			pi->fviscosity = pi->fviscosity + calcSubviscosity(pi, pj)*pj->mass*laplacianKernel(rij, q, viscosity);//JS
	//		}
	//	}
	//	pi->fpressure = pi->fpressure * -1.0f * cspikyG;//JS
	//	pi->fviscosity = pi->fviscosity * u * cviscosityL;//JS
	//}
	for (int i = 1; i < HASHSIZE - 1; i++)
	{
		for (int j = 1; j < HASHSIZE - 1; j++)
		{
			vector<Particle*> pis;
			vector<Particle*> pjs = GetNeighbor(i, j, h, pis);
			for (auto &pi : pis)
			{
				pi->fpressure = Vector2D(0.0f, 0.0f);
				pi->fviscosity = Vector2D(0.0f, 0.0f);
				for (auto &pj : pjs)
				{
					if (pi->idx == pj->idx)
						continue;

					Vector2D rij = pi->position - pj->position;
					float q = rij.norm() / h;
					if (0 < q && q <= 1.0f)
					{
						pi->fpressure = pi->fpressure + gradientKernel(rij, q, spiky)*pj->mass*calcSubpressure(GAS_CONSTANT, pi, pj, rest_density);//JS
						pi->fviscosity = pi->fviscosity + calcSubviscosity(pi, pj)*pj->mass*laplacianKernel(rij, q, viscosity);//JS
					}
				}
				pi->fpressure = pi->fpressure * -1.0f * cspikyG;//JS
				pi->fviscosity = pi->fviscosity * u * cviscosityL;//JS
			}
		}
	}
	
}

void SPH::computeSurfacetension()
{
	////color field
	//for (auto &pi : particles)
	//{
	//	pi->colorFnormal = Vector2D(0.0f, 0.0f);
	//	float tempSurface = 0.0f;
	//	vector<int> neighbors = GetNeighbor(pi->getPosX(), pi->getPosY(), h);
	//	for (auto &j : neighbors)
	//	{
	//		Particle *&pj = particles.at(j);
	//		Vector2D rij = pi->position - pj->position;

	//		float norm = rij.norm();
	//		float q = norm / h;
	//		if (q <= 1)
	//		{
	//			float tmp = pj->mass / pj->density;
	//			pi->colorFnormal = pi->colorFnormal + gradientKernel(rij, norm, poly6) * tmp;
	//			tempSurface = tempSurface + laplacianKernel(rij, norm, poly6) * tmp;
	//		}
	//	}
	//	tempSurface = tempSurface * cpoly6L;
	//	pi->colorFnormal = pi->colorFnormal * cpoly6G;
	//	if (pi->colorFnormal.norm() > 0.0f)
	//		pi->fsurface = pi->colorFnormal / pi->colorFnormal.norm() * tempSurface * surfaceCoef * (-1.0f);
	//}
	//color field
	for (int i = 1; i < HASHSIZE - 1; i++)
	{
		for (int j = 1; j < HASHSIZE - 1; j++)
		{
			vector<Particle*> pis;
			vector<Particle*> pjs = GetNeighbor(i, j, h, pis);
			for (auto &pi : pis)
			{
				pi->colorFnormal = Vector2D(0.0f, 0.0f);
				float tempSurface = 0.0f;
				for (auto &pj : pjs)
				{
					Vector2D rij = pi->position - pj->position;

					float norm = rij.norm();
					float q = norm / h;
					if (q <= 1)
					{
						float tmp = pj->mass / pj->density;
						pi->colorFnormal = pi->colorFnormal + gradientKernel(rij, norm, poly6) * tmp;
						tempSurface = tempSurface + laplacianKernel(rij, norm, poly6) * tmp;
					}
				}
				tempSurface = tempSurface * cpoly6L;
				pi->colorFnormal = pi->colorFnormal * cpoly6G;
				if (pi->colorFnormal.norm() > 0.0f)
					pi->fsurface = pi->colorFnormal / pi->colorFnormal.norm() * tempSurface * surfaceCoef * (-1.0f);
			}
		}
	}
}

void SPH::Integrate(float dt, Vector2D gravity)
{
	for (auto &p : particles)
	{
		p->integrate(dt, gravity);
	}
}

/*void SPH::HashTable()
{
	for (int i = 0; i < HASHSIZE; i++)
	{
		for (int j = 0; j < HASHSIZE; j++)
		{
			hashGrid[i][j].clear();
		}
	}
	int i = 0;
	for (auto &p : particles)
	{
		float x = (p->getPosX() + (float)view_width) / ((float)view_width * 2.0f);
		float y = (p->getPosY() + (float)view_height) / ((float)view_height * 2.0f);
		int gridx = (int)(GRIDSIZE * x) + 1;
		int gridy = (int)(GRIDSIZE * y) + 1;

		if (gridx < 1) gridx = 1;
		if (gridx > HASHSIZE - 2) gridx = HASHSIZE - 2;
		if (gridy < 1) y = 1;
		if (gridy > HASHSIZE - 2) gridy = HASHSIZE - 2;

		hashGrid[gridx][gridy].push_back(p->idx);
	}
}*/

/*vector<int> SPH::GetNeighbor(float x, float y, float radius)
{
	float tempx = (x + (float)view_width) / ((float)view_width * 2.0f);
	float tempy = (y + (float)view_height) / ((float)view_height * 2.0f);
	int gridx = (int)(GRIDSIZE * tempx) + 1;
	int gridy = (int)(GRIDSIZE * tempy) + 1;

	vector<int>res;
	for (int i = gridx - radius; i <= gridx + radius; i++)
	{
		for (int j = gridy - radius; j <= gridy + radius; j++)
		{
			if (i<0 || i>HASHSIZE - 1 || j<0 || j>HASHSIZE - 1)continue;
			for (int k = 0; k < hashGrid[i][j].size(); k++)
			{
				res.push_back(hashGrid[i][j].at(k));
			}
		}
	}
	return res;
}*/
void SPH::HashTable()
{
	for (int i = 0; i < HASHSIZE; i++)
	{
		for (int j = 0; j < HASHSIZE; j++)
		{
			hashGrid[i][j].clear();
		}
	}
	int i = 0;
	for (auto &p : particles)
	{
		float x = (p->getPosX() + (float)view_width) / ((float)view_width * 2.0f);
		float y = (p->getPosY() + (float)view_height) / ((float)view_height * 2.0f);
		int gridx = (int)(GRIDSIZE * x) + 1;
		int gridy = (int)(GRIDSIZE * y) + 1;

		if (gridx < 1) gridx = 1;
		if (gridx > HASHSIZE - 2) gridx = HASHSIZE - 2;
		if (gridy < 1) y = 1;
		if (gridy > HASHSIZE - 2) gridy = HASHSIZE - 2;

		hashGrid[gridx][gridy].push_back(p);
	}
}

vector<Particle *> SPH::GetNeighbor(float x, float y, float radius)
{
	float tempx = (x + (float)view_width) / ((float)view_width * 2.0f);
	float tempy = (y + (float)view_height) / ((float)view_height * 2.0f);
	int gridx = (int)(GRIDSIZE * tempx) + 1;
	int gridy = (int)(GRIDSIZE * tempy) + 1;

	return GetNeighbor(gridx, gridy, radius);
}

vector<Particle *> SPH::GetNeighbor(int gridx, int gridy, float radius)
{
	vector<Particle *>res;

	for (int i = gridx - radius; i <= gridx + radius; i++)
	{
		for (int j = gridy - radius; j <= gridy + radius; j++)
		{
			if (i<0 || i>HASHSIZE - 1 || j<0 || j>HASHSIZE - 1)continue;
			for (int k = 0; k < hashGrid[i][j].size(); k++)
			{
				res.push_back(hashGrid[i][j].at(k));
			}
		}
	}
	return res;
}

vector<Particle *> SPH::GetNeighbor(int gridx, int gridy, float radius, vector<Particle*> &mine)
{
	vector<Particle *>res;
	mine.clear();
	for (int i = gridx - radius; i <= gridx + radius; i++)
	{
		for (int j = gridy - radius; j <= gridy + radius; j++)
		{
			if (i<0 || i>HASHSIZE - 1 || j<0 || j>HASHSIZE - 1)continue;
			for (int k = 0; k < hashGrid[i][j].size(); k++)
			{
				res.push_back(hashGrid[i][j].at(k));
				if (i == gridx && j == gridy)
				{
					mine.push_back(hashGrid[i][j].at(k));
				}
			}
		}
	}
	return res;
}

vector<Particle *> SPH::GetMine(int gridx, int gridy, float radius)
{
	vector<Particle *>res;
	if (gridx >= 0 && gridx <= HASHSIZE - 1 && gridy >= 0 && gridy <= HASHSIZE - 1)
	{
		for (int k = 0; k < hashGrid[gridx][gridy].size(); k++)
		{
			res.push_back(hashGrid[gridx][gridy].at(k));
		}
	}
	
	return res;
}