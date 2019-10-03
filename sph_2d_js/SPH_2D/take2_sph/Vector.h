#pragma once
#include <iostream>

class Vector2D {
public:
	float x;
	float y;
public:
	Vector2D():x(0.0f), y(0.0f) { }
	Vector2D(float _x, float _y) : x(_x), y(_y) { }
	~Vector2D() {}

	float getX() { return x; }
	float getY() { return y; }
	float setX(float _x) { this->x = _x; }
	float setY(float _y) { this->y = _y; }

	float norm() 
	{
		return sqrt(x*x + y*y);
	}
	float normsqr()
	{
		return x*x + y*y;
	}
	float dist(Vector2D v)
	{
		return sqrt((x - v.x)*(x - v.x) + (y - v.y)*(y - v.y));
	}
	void operator=(Vector2D &v)
	{
		x = v.x;
		y = v.y;
	}

	Vector2D operator+(Vector2D &v)
	{
		Vector2D result;
		result.x = x + v.x;
		result.y = y + v.y;
		return result;
	}

	Vector2D operator-(Vector2D &v)
	{
		Vector2D result;
		result.x = x - v.x;
		result.y = y - v.y;
		return result;
	}

	Vector2D operator*(float f)
	{
		Vector2D result;
		result.x = x * f;
		result.y = y * f;
		return result;
	}

	float operator*(Vector2D &v)
	{
		return this->x*v.x + this->y*v.y;
	}

	Vector2D operator/(float f)
	{
		Vector2D result;
		result.x = x / f;
		result.y = y / f;
		return result;
	}

	Vector2D &operator+=(Vector2D &v)
	{
		x += v.x;
		y += v.y;
		return *this;
	}

	Vector2D &operator-=(Vector2D &v)
	{
		x -= v.x;
		y -= v.y;
		return *this;
	}

	bool operator==(Vector2D &v)
	{
		return (x == v.x && y == v.y)? true : false;
	}

	Vector2D &operator*=(float f)
	{
		x *= f;
		y *= f;
		return *this;
	}
};