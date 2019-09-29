//
//#include<Windows.h>
//#include<math.h>
//
//#include<vector>
//#include<GL/GL.h>
////#include<GL/glu_.h>
//#include<GL/glut.h>
//#include"vector.h"
//#include<iostream>
// 
//
//using namespace std;
//
//float	pi = 3.1415926;
//float	dt=0.015;
//int		num_of_particles = 1600;
//bool	StartSimulation=false;
//
//float h=1;
//float k_1=100; 
//float restDensity=1;
//float visc=10;
//float surfacetension = 10.0;
//#define	HashSize 100
//vector<int>		hashGrid[HashSize][HashSize];
//
//struct index2D{
//	int x, y;
//};
//
//
//struct particle{
// 
//	vector2D	Position;
//	vector2D	Velocity;
//	vector2D    f_pressure;//3
//	vector2D    f_viscosity;
//	vector2D	f_surface;
//	vector2D    Acceleration;
//
//	vector2D		colorF;
//	float		Density;//1
//	float       Mass;
//	float       Pressure;//2
//
//	vector2D		n;
//
// 
//};
//
//particle *waterParticle;
// 
//void HashTable()
//{
//	for (int x = 0; x < HashSize; x++)
//	{
//		for (int y = 0; y < HashSize; y++)
//		{ 
//				hashGrid[x][y].clear(); 
//		}
//	}
//
//	for (int i = 0; i <num_of_particles; i++)
//	{
//		int x = waterParticle[i].Position.x;
//		int y = waterParticle[i].Position.y;
// 
//
//		if (x < 0) x = 0;
//		if (x > HashSize - 1)x = HashSize - 1;
//
//		if (y < 0) y = 0;
//		if (y > HashSize - 1)y = HashSize - 1;
//
//		hashGrid[x][y].push_back(i);
//	}
//}
//
//
//vector<int>GetNeighbor(float x, float y, float radius)
//{
//	vector<int>res;
//	for (int x0 = x - radius; x0 <= x + radius; x0++)
//	{
//		for (int y0 = y - radius; y0 <=  y + radius; y0++)
//		{
//
//			if (x0 <  0 || x0>HashSize - 1 || y0 <  0 || y0>HashSize - 1) continue;
//
//			for (int i = 0; i < hashGrid[x0][y0].size(); i++)
//			{
//				int j = hashGrid[x0][y0][i];
//				 
//				res.push_back(j);
//			}
//		}
//	}
//	return res;
//}
//
// 
//void init()
//{
//	waterParticle = (particle*)malloc(num_of_particles*sizeof(particle));
//	memset(waterParticle, 0, num_of_particles*sizeof(particle));
//
//	for (int i = 0; i < 40; i++)
//	{
//		for (int j = 0; j < 40; j++)
//		{
//			waterParticle[i + 40 * j].Position.x = i / 1.02  +1;
//			waterParticle[i + 40 * j].Position.y = j / 1.02  +3;
//		}
//	}
//
//	for (int i = 0; i < num_of_particles; i++)
//	{
//		waterParticle[i].f_surface = vector2D(0, 0);
//	}
//}
//float Kernel(vector2D x, float h)
//{
//	float q;
//	float result;
//	float constant;
//	constant= 40.0f / (h*h * 7.0f * pi);
//
//	q = x.Magnitude() / h;
//
//
//	if (q >= 0.0f && q <= 0.5)
//	{
//		result = constant * (6.0f *pow(q,3) - 6.0f * (q*q) + 1.0f);
//	 
//		return result;
//	}
//	else if (q >= 0.5 && q <= 1.0f)
//	{
//		result = constant * (2.0f * pow(1 - q, 3));
//		return result;
//	}
//	else if (q > 1.0f)
//	{
//		return 0;
//	}
// 
//
//}
//vector2D GradientKernel(vector2D x, float h)
//{
//	float q;
//	vector2D gradientq;
//	float constant;
//	constant= 6.0f*40.0f / (h*h * 7.0f* pi);
//
//	q = x.Magnitude() / h ;
//	gradientq = x / (x.Magnitude()*h );
//	 
//
//	if (q > 0 && q <= 0.5	)
//	{
//
//		gradientq = gradientq *constant* (3.0f* q*q - 2.0f* q) ;
//	 
//		return gradientq;
//	}
//	else if (q >= 0.5 && q <= 1.0f)
//	{
//		gradientq = gradientq  *constant* (-(1.0f - q)*(1.0f- q));
//	 
//		return gradientq;
//	}
//	else if (q > 1.0f||q==0   )
//	{
//		gradientq  = vector2D(0.0,0.0);
//	 
//		return gradientq;
//	}
//	
//	
//}
//
//float LapKernel(vector2D x,float h)
//{
//	float q;
//	float laplacianq;
//	float result;
//	float constant;
//	constant= 6.0f*40.0f / (h*h * 7.0f* pi);
//	q = x.Magnitude() / h;
//	laplacianq = 2.0f / (x.Magnitude() *h + 0.0000000001);
//
//	if (q >= 0.0 && q <=0.5)
//	{
//	
//		result = constant*((6.0f* q - 2.0f) / (h*h) + (3.0f* q*q - 2.0f* q) *laplacianq);
//	 
//		return result;
//	}
//	else if (q >= 0.5f && q <= 1.0f)
//	{
//	 
//		result = constant*((2.0f * (1.0f - q)) / (h*h) + (-(1.0f - q)*(1.0f - q))*laplacianq);
//		return result;
//	}
// 
//	else if (q>1)
//	{
//		return 0;
//	}
// 
//}
//
///*
//float Kernel(vector2D r, float h)
//{
//	float q;
//	float result;
//	q = sqrt((r.x)*(r.x) + (r.y)*(r.y));
//	float constant;
//	constant = 4.0 / pi*pow(h, 8);
//	if (q >= 0.0f && q <= h){
//		result = constant*pow(h*h - q*q, 3);
//		return result;
//	}
//	else{ return 0; }
//
//}
//vector2D GradientKernel(vector2D r, float h)
//{
//	float q;
//	float result;
//	vector2D gradientq;
//	q = sqrt((r.x)*(r.x) + (r.y)*(r.y));
//	float constant;
//	constant = -30.0f / pi*pow(h, 5);
//	
//	if (q >= 0.0f && q <= h){
//		gradientq.x = r.x / (sqrt((r.x)*(r.x) + (r.y)*(r.y)));
//		gradientq.y = r.y / (sqrt((r.x)*(r.x) + (r.y)*(r.y)));
//		gradientq.x = constant*pow(h - q, 2)*gradientq.x;
//		gradientq.y = constant*pow(h - q, 2)*gradientq.y ;
//		return gradientq;
//	}
//	else{
//		gradientq.x = 0.0;
//		gradientq.y = 0.0;
//		return gradientq;
//	}
//}
//float LapKernel(vector2D r, float h)
//{
//	float q;
//	float result;
//	float constant;
//	q = sqrt((r.x)*(r.x) + (r.y)*(r.y));
//	constant = 20.0f / 3.0f * pi*pow(h, 5);
//	if (q >= 0.0f && q <= h){
//	result = constant*(h - q);
//	}
//	else{
//	return 0.0;
//	}
//	}*/
//void computeDensity(float rest_density, float l, float h)
//{
//	for (int i = 0; i < num_of_particles; i++)
//	{
//		vector2D r_i = waterParticle[i].Position;
//
//		waterParticle[i].Density = 0.0f;
//		float dens = 0.0;
//		vector<int> res = GetNeighbor(waterParticle[i].Position.x, waterParticle[i].Position.y, h);
//
//		for (int k = 0; k < res.size(); k++)
//		{
//			int j = res[k];
//
//			vector2D r_j = waterParticle[j].Position;
//
//			vector2D   x1 = r_i - r_j;
//			waterParticle[j].Mass = 1.0f;
//			dens = dens + waterParticle[j].Mass*Kernel(x1, h);
//		}
//		waterParticle[i].Density = dens;
//		waterParticle[i].Pressure = l*(waterParticle[i].Density - rest_density);
//
//		// if (waterParticle[i].Density - rest_density <0.0) waterParticle[i].Pressure = 0.0;
//	}
//}
//
//void computePressureForce(float h){
//
//	for (int i = 0; i < num_of_particles; i++)
//	{
//		vector2D r_i = waterParticle[i].Position;
//		waterParticle[i].f_pressure = vector2D(0, 0);
//
//		vector<int> res = GetNeighbor(waterParticle[i].Position.x, waterParticle[i].Position.y, h);
//		float tmp = 0.0;
//		vector2D tmpGK;
//		for (int k = 0; k < res.size(); k++)
//		{
//			int j = res[k];
//
//			waterParticle[j].Mass = 1.0f;
//			vector2D r_j = waterParticle[j].Position;
//			vector2D   x1 = r_i - r_j;
//
//			tmp = tmp + waterParticle[j].Mass*((waterParticle[i].Pressure + waterParticle[j].Pressure) / (2.0f* waterParticle[j].Density));
//
//			tmpGK = tmpGK + GradientKernel(x1, h);
//		}
//
//		waterParticle[i].f_pressure = -1 * tmp*tmpGK;
//
//	}
//}
//
//void computeViscosityForce(float h, float visc)
//{
//	for (int i = 0; i < num_of_particles; i++)
//	{
//		vector2D r_i = waterParticle[i].Position;
//
//		waterParticle[i].f_viscosity = vector2D(0, 0);
//		float tmpLK;
//		vector2D tmp;
//		vector<int> res = GetNeighbor(waterParticle[i].Position.x, waterParticle[i].Position.y, h);
//
//
//		for (int k = 0; k < res.size(); k++)
//		{
//			int j = res[k];
//
//			waterParticle[j].f_viscosity = vector2D(0, 0);
//			vector2D r_j = waterParticle[j].Position;
//
//			vector2D   x1 = r_i - r_j;
//
//			tmpLK = visc*waterParticle[j].Mass* LapKernel(x1, h) / waterParticle[j].Density;
//			tmp = tmp + (waterParticle[j].Velocity - waterParticle[i].Velocity)*tmpLK;
//
//
//		}
//		waterParticle[i].f_viscosity = tmp;
//	}
//}
//
//
////void computeSurfaceForce(float h)
////{
////	/*for (int i = 0; i < num_of_particles; i++)
////	{
//// 
////		vector2D surf=vector2D(0,0);
////		vector2D r_i = waterParticle[i].Position;
////		waterParticle[i].f_surface = vector2D(0.0, 0.0);
////		vector<int> res = GetNeighbor(waterParticle[i].Position.x, waterParticle[i].Position.y, h);
////		for (int k = 0; k < res.size(); k++)
////		{
////			int j = res[k];
////			vector2D r_j = waterParticle[j].Position;
////			vector2D   x1 = r_i - r_j;
////
////			x1.Normalize();
////	 
////			surf = surf + x1* (-1 * surfacetension)*waterParticle[j].Mass*LapKernel(x1, h) / (waterParticle[j].Density + 0.01);
////
////		}
//// 	waterParticle[i].f_surface = surf;
////
////	}*/
////}
//
//
//void computeSurfaceForce(float h)
//{
//	
//	for (int i = 0; i < num_of_particles; i++)
//	{
//		vector2D r_i = waterParticle[i].Position;
//
//		waterParticle[i].f_surface = vector2D(0.0, 0.0);
//		vector2D colorF=vector2D(0,0);	 
//		vector<int> res = GetNeighbor(waterParticle[i].Position.x, waterParticle[i].Position.y, h);
//		for (int k = 0; k < res.size(); k++)
//		{
//			int j = res[k];
//			vector2D r_j = waterParticle[j].Position;
//			vector2D   x1 = r_i - r_j;
//			float dist = x1.Magnitude();
//			if (dist < h)
//			{
//				colorF = colorF + waterParticle[j].Mass*GradientKernel(x1,h) / waterParticle[j].Density;
//			}
//		}
//		waterParticle[i].colorF = colorF;
//	}
//
//	/*for (int i = 0; i < num_of_particles; i++)
//	{
//		vector2D r_i = waterParticle[i].Position;
//
//		waterParticle[i].n = vector2D(0.0, 0.0);
//		vector2D nor = vector2D(0, 0);
//		float tmp = 0.0;
//		vector2D tmpGK  ;
//		vector<int> res = GetNeighbor(waterParticle[i].Position.x, waterParticle[i].Position.y, h);
//
//		for (int k = 0; k < res.size(); k++)
//		{
//			int j = res[k];
//			vector2D r_j = waterParticle[j].Position;
//
//			vector2D   x1 = r_i - r_j;
//			float dist = x1.Magnitude();
//			if (x1.Magnitude < h)
//			{
//			
//				tmp = tmp + waterParticle[j].Mass*((waterParticle[i].colorF + waterParticle[j].colorF) / (2.0f* waterParticle[j].Density));
//
//				tmpGK = tmpGK + GradientKernel(x1, h);
//			}
//			 
//		}
//		waterParticle[i].n = tmpGK;
//	}*/
//
//	for (int i = 0; i < num_of_particles; i++)
//	{
//		vector2D r_i = waterParticle[i].Position;
//
//		waterParticle[i].f_surface = vector2D(0.0, 0.0);
//		vector2D tmpSur = vector2D(0, 0);
//		vector2D cohesionfor = vector2D(0.0, 0.0);
//		vector2D curvaturefor = vector2D(0.0, 0.0);
//		waterParticle[i].Mass=1.0;
//		vector<int> res = GetNeighbor(waterParticle[i].Position.x, waterParticle[i].Position.y, h);
//
//		for (int k = 0; k < res.size(); k++)
//		{
//			int j = res[k];
//			vector2D r_j = waterParticle[j].Position;
//
//			vector2D   x1 = r_i - r_j;
//			float dist = x1.Magnitude();
//			if (dist < h)
//			{
//				float kkk = 2.0*restDensity / (waterParticle[i].Density + waterParticle[j].Density);
// 
//				curvaturefor = -200*waterParticle[i].Mass*(waterParticle[i].colorF  -waterParticle[j].colorF);
//				/*float tmpLK = waterParticle[j].Mass* LapKernel(x1, h) / waterParticle[j].Density;
//				tmpSur = tmpSur + (waterParticle[j].colorF - waterParticle[i].colorF)*tmpLK*(waterParticle[i].n / waterParticle[i].n.Magnitude());*/
//				tmpSur = tmpSur + kkk*curvaturefor;
//			 
//			} 
//		}
//		waterParticle[i].f_surface =  1*tmpSur;
//		//cout << waterParticle[i].f_surface.x << "   " << waterParticle[i].f_surface.y << endl;
//	}
//}
//
//float radius =30.0;
//
//float getDistance(vector2D pos){
//	vector2D center=vector2D(0.0,0.0);
//	return pow(pos.x - center.x, 2.0) + pow(pos.y - center.y, 2.0) - (radius*radius);
//}
//
//vector2D getGradient(vector2D pos)
//{
//	vector2D gradient;
//	float dx = 0.001;
//	gradient.x = getDistance(vector2D(pos.x + dx, pos.y) - vector2D(pos.x - dx, pos.y));
//	gradient.y = getDistance(vector2D(pos.x, pos.y + dx) - vector2D(pos.x, pos.y - dx));
//	return gradient;
//}
//
//
//
//bool IsInside(vector2D pos)
//{
//	return getDistance(pos) <= 0.0;
//}
//
//void position()
//{
//	//float g = -1.0f;
//	vector2D g = vector2D(0, -1.9);
//	for (int i = 0; i < num_of_particles; i++)
//	{
//		
//
//		waterParticle[i].Acceleration = vector2D(0, 0);
//		waterParticle[i].Acceleration = g + (waterParticle[i].f_pressure + waterParticle[i].f_viscosity + waterParticle[i].f_surface) / waterParticle[i].Density;
//		waterParticle[i].Velocity = waterParticle[i].Velocity + waterParticle[i].Acceleration*dt;
//		waterParticle[i].Position = waterParticle[i].Position + waterParticle[i].Velocity*dt;
//
//		/*if (!IsInside(waterParticle[i].Position)){
//			vector2D normal = getGradient(waterParticle[i].Position);
//			float update_vel = normal.x*(-waterParticle[i].Velocity.x) + normal.y*(-waterParticle[i].Velocity.y);
//			vector2D n1;
//			n1 = update_vel*normal;
//			waterParticle[i].Velocity.x = 2.0 * n1.x + waterParticle[i].Velocity.x;
//			waterParticle[i].Velocity.y = 2.0 * n1.y + waterParticle[i].Velocity.y;
//			waterParticle[i].Position.x = waterParticle[i].Position.x*( radius/sqrt(pow(waterParticle[i].Position.x, 2) + pow(waterParticle[i].Position.y, 2)));
//			waterParticle[i].Position.y = waterParticle[i].Position.y* (radius/sqrt(pow(waterParticle[i].Position.x, 2) + pow(waterParticle[i].Position.y, 2)));
//			}
//
//			/*double getDistance(vec3 pos)
//			{
//			return pow(pos.x - sphere.cp.x, 2.0) + pow(pos.y - sphere.cp.y, 2.0) + pow(pos.z - sphere.cp.z, 2.0) - (sphere.r*sphere.r);
//			}
//
//			vec3 getGradient(vec3 pos)
//			{
//			vec3 grad;
//			double dx = 0.001;
//			grad.x = getDistance(vec3(pos.x + dx, pos.y, pos.z)) - getDistance(vec3(pos.x - dx, pos.y, pos.z));
//			grad.y = getDistance(vec3(pos.x, pos.y + dx, pos.z)) - getDistance(vec3(pos.x - dx, pos.y, pos.z));
//			grad.z = getDistance(vec3(pos.x + dx, pos.y, pos.z)) - getDistance(vec3(pos.x - dx, pos.y, pos.z));
//			return grad;
//			}
//
//			bool IsInside(vec3 pos)
//			{
//			return getDistance(pos) < 0.0;
//			}
//
//			if (!sphere.IsInside(particle.pos)) {
//			// derivative
//			vec3 n = getGradient(pos);
//
//			float a = n.x*(-waterParticle[i].Velocity.x) + n.y*(-waterParticle[i].Velocity.y);
//			vector2D n1;
//			n1 = a*n;
//			waterParticle[i].Velocity.x = 2.0 * n1.x + waterParticle[i].Velocity.x;
//			waterParticle[i].Velocity.y = 2.0 * n1.y + waterParticle[i].Velocity.y;
//			//waterParticle[i].Position.y = 0.01f;
//			}
//			*/
//
//		float para = 1.0;
//		if (waterParticle[i].Position.y < 0.0f)
//		{
//			vector2D n;
//			n = vector2D(0.0, 1.0);
//			float a = n.x*(-waterParticle[i].Velocity.x) + n.y*(-waterParticle[i].Velocity.y);
//			vector2D n1;
//			n1 = a*n;
//			waterParticle[i].Position.y = 0.0f;
//			 waterParticle[i].Velocity.y = para * n1.y + waterParticle[i].Velocity.y;
//		 
//		}
//
//		if (waterParticle[i].Position.y > 99)
//		{
//			vector2D n;
//			n = vector2D(0.0, -1.0);
//			float a = n.x*(-waterParticle[i].Velocity.x) + n.y*(-waterParticle[i].Velocity.y);
//			vector2D n1;
//			n1 = a*n;
//			waterParticle[i].Position.y = 99;
//			//	waterParticle[i].Velocity.x = 2.0 * n1.x + waterParticle[i].Velocity.x;
//			waterParticle[i].Velocity.y = para * n1.y + waterParticle[i].Velocity.y;
// 
//		}
//
//		if (waterParticle[i].Position.x < 0.1f)
//		{
//			vector2D n;
//			n = vector2D(1.0, 0.0);
//			float a = n.x*(-waterParticle[i].Velocity.x) + n.y*(-waterParticle[i].Velocity.y);
//			vector2D n1;
//			n1 = a*n;
//			 	waterParticle[i].Velocity.x = 0.7*para * n1.x + waterParticle[i].Velocity.x;
//			waterParticle[i].Position.x = 0.01f;
//		 
//		}
//
//
//		if (waterParticle[i].Position.x > 99.0f)
//		{
//
//			vector2D n;
//			n = vector2D(-1.0, 0.0);
//
//			float a = n.x*(-waterParticle[i].Velocity.x) + n.y*(-waterParticle[i].Velocity.y);
//			vector2D n1;
//			n1 = a*n;
//
//			 	waterParticle[i].Velocity.x = para * n1.x + waterParticle[i].Velocity.x;
// 
//			waterParticle[i].Position.x = 99.1f;
//	 
//		}
//	}
//}
//
//void SPH(float h, float k, float restDensity, float visc)
//{
//	HashTable();
//	if (StartSimulation == true)
//	{
//	
//		computeDensity(restDensity, k, h);
//		computePressureForce(h);
//		computeViscosityForce(h, visc);
//		computeSurfaceForce(h);
//		position();
//	}
//}
//
//void reshape(int width, int height)
//{
//	glViewport(0.0, 0.0, width, height);
//	glMatrixMode(GL_PROJECTION);
//	glLoadIdentity();
//	glOrtho(0.0, 200, 0.0, 200,-200,200);
//}
//void Keyboard(unsigned char key, int x, int y)
//{
//	switch (key)
//	{
//	case 'q':
//	case VK_ESCAPE:
//		exit(0);
//
//	case' ':
//		StartSimulation =! StartSimulation;
//		break;
//	}
//}
//void display()
//{
//	glClear(GL_COLOR_BUFFER_BIT);
//	glMatrixMode(GL_MODELVIEW);
//	glLoadIdentity();	
//
//
//	glTranslatef(50,50,-10);
//	glPointSize(2.0);
//	
//	glBegin(GL_POINTS);
//	for (int i = 0; i < num_of_particles; i++){
//    	glColor3f(1.0,1.0,1.0);
//		glVertex2f(waterParticle[i].Position.x, waterParticle[i].Position.y);
//	}
//	glEnd();
//   /*glBegin(GL_POINTS);
//	for (int i = 0; i < num_of_particles; i++){
//		glColor3f(0.0, 0.0, 1.0);
//		glVertex2f(waterParticle[i].HashIndex.x, waterParticle[i].HashIndex.y);
//		
//	}
//	glEnd();*/
//	glColor3f(1.0, 0.0, 0.0);
//	glBegin(GL_LINE_LOOP);
//	glVertex2f(0, 0);
//	glVertex2f(100, 0);
//	glVertex2f(100, 100);
//	glVertex2f(0, 100);
//	glEnd();
//
//	SPH(h,k_1,restDensity,visc);
//	glutSwapBuffers();
//	glutPostRedisplay();
//}
//
//int main(int argc, char **argv)
//{
//	glutInit(&argc, argv);
//	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
//	glutInitWindowPosition(2000, 0);
//	glutInitWindowSize(800, 800);
//	glutCreateWindow("Fluid");
//	glutReshapeFunc(reshape);
//	glutDisplayFunc(display);
//	glutKeyboardFunc(Keyboard);
//	init();
//	glutMainLoop();
//}