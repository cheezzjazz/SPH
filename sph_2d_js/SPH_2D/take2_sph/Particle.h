#pragma once
#include <GL/glut.h>
#include "Vector.h"

class Particle {
public:	//common
	enum mode
	{
		node = 0, 
		fluid
	};
	int mass;
	Vector2D position;
	Vector2D velocity;
	Vector2D accel;
	mode type;
public://spring
	Vector2D force;
//	float air_friction;
//	Vector2D normal;
//	bool isFixed;
public://SPH
	float density;
	Vector2D colorFnormal;
	int idx;
	Vector2D fpressure;
	Vector2D fviscosity;
	Vector2D fsurface;
	float pointR = 0.0f;//1.0f;
	float view_height = 20;	//main 이랑 연결해야함
	float view_width = 20;

	Particle();
	Particle(float _x, float _y, mode m);
	Particle(float _x, float _y, int _idx, mode m);
	~Particle();
	float getPosX() { return position.getX(); }
	float getPosY() { return position.getY(); }

	void add_force(Vector2D ext_force);

	void integrate(float dt, Vector2D gravity);

	void draw();
};