#include "Particle.h"

Particle::Particle()
{
}

Particle::Particle(float _x, float _y, mode m) : position(_x, _y), velocity(0.0f, 0.0f), accel(0.0f, 0.0f), mass(1), type(m)
{
	force = Vector2D(0.0f, 0.0f);
	fpressure = Vector2D(0.0f, 0.0f);
	fviscosity = Vector2D(0.0f, 0.0f);
	fsurface = Vector2D(0.0f, 0.0f);
	colorFnormal = Vector2D(0.0f, 0.0f);
	density = 0.0f;
	idx = -1;
}
Particle::Particle(float _x, float _y, int _idx, mode m) : position(_x, _y), velocity(0.0f, 0.0f), accel(0.0f, 0.0f), mass(1), type(m)
{
	force = Vector2D(0.0f, 0.0f);
	fpressure = Vector2D(0.0f, 0.0f);
	fviscosity = Vector2D(0.0f, 0.0f);
	colorFnormal = Vector2D(0.0f, 0.0f);
	fsurface = Vector2D(0.0f, 0.0f);
	density = 0.0f;
	idx = _idx;
}
Particle::~Particle()
{

}

void Particle::add_force(Vector2D ext_force)
{
	force += ext_force;
}

void Particle::integrate(float dt, Vector2D gravity)
{
	Vector2D fgrav = gravity * mass;
	//Vector2D fgrav = p.Mass * Vector2D(0.0f, -1.9f);
	if (density == 0.0f)
	{
		density = 0.0001f;
		std::cout << "integrate(000!)" << std::endl;
	}
	accel = (fpressure + fviscosity + fsurface) / density + fgrav;
	velocity = velocity + accel*dt;
	position = position + velocity*dt;

	float restitution = 0.6f;//0.3f;
	float damping = 0.5f;
	if (position.x - pointR < -view_width)
	{
		Vector2D n = Vector2D(1.0f, 0.0f);
		float length = abs((velocity.x)*n.x + (velocity.y)*n.y);//JS
																//float length = (-p.velocity.x) * n.x + (-p.velocity.y) * n.y;//QM
		Vector2D n1 = n*length;
		velocity = velocity + n1*(1.0f + restitution);//JS
													  //p.velocity.x = 0.7 * n1.x + p.velocity.x;//QM
													  //p.velocity.x = p.velocity.x*damping;
		position.x = -view_width + pointR;
	}
	if (position.x + pointR > view_width)
	{
		Vector2D n = Vector2D(-1.0f, 0.0f);
		float length = abs(velocity.x*n.x + velocity.y*n.y);//JS
															//float length = (-p.velocity.x) * n.x + (-p.velocity.y) * n.y;//QM
		Vector2D n1 = n*length;
		velocity = velocity + n1*(1.0f + restitution);//JS
													  /*p.velocity.x = p.velocity.x*damping;*/
													  //p.velocity.x = n1.x + p.velocity.x;//QM
		position.x = view_width - pointR;
	}
	if (position.y - pointR < -view_height)
	{
		Vector2D n = Vector2D(0.0f, 1.0f);
		float length = abs(velocity.x*n.x + velocity.y*n.y);//JS
															//float length = (-p.velocity.x) * n.x + (-p.velocity.y) * n.y;//QM

		Vector2D n1 = n*length;
		velocity = velocity + n1*(1.0f + restitution);//JS
													  //p.velocity.y = p.velocity.y*damping;
													  //p.velocity.y = n1.y + p.velocity.y;//QM
		position.y = -view_height + pointR;
	}
	if (position.y + pointR > view_height)
	{
		Vector2D n = Vector2D(0.0f, -1.0f);
		float length = abs(velocity.x*n.x + velocity.y*n.y);//JS
															//float length = (-p.velocity.x) * n.x + (-p.velocity.y) * n.y;//QM
		Vector2D n1 = n*length;
		velocity = velocity + n1*(1.0f + restitution);//JS
													  //p.velocity.y = p.velocity.y*damping;
													  //p.velocity.y = n1.y + p.velocity.y;//QM
		position.y = view_height - pointR;
	}

}

void Particle::draw()
{
	glColor3f(1.0f, 1.0f, 1.0f);
	glPointSize(5.0f);
	glEnable(GL_POINT_SMOOTH);
	glBegin(GL_POINTS);
	glVertex2f(getPosX(), getPosY());
	glEnd();

}
