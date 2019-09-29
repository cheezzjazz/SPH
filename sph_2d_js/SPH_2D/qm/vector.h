class vector2D
{

public:
	vector2D(){
		x = y = 0.0f;
	};
	vector2D(float x, float y){
		this->x = x;
		this->y = y;
	};
	~vector2D(){};



	float x;
	float y;


	vector2D operator+(vector2D v)
	{
		vector2D result = (*this);
		result.x += v.x;
		result.y += v.y;

		return result;
	}

	 
	vector2D operator-(vector2D v)
	{
		vector2D result = (*this);
		result.x -= v.x;
		result.y -=v.y;

		return result;
	}

	vector2D operator*(float val)
	{
		vector2D result = (*this);
		result.x *= val;
		result.y *= val;
		return result;
	}
 
	vector2D operator/(float val)
	{
		vector2D result = (*this);
		result.x /= val;
		result.y /= val;

		return result;
	}
	 
	vector2D operator=(float val)
	{
		vector2D result = (*this);
		result.x = val;
		result.y = val;

		return result;
	}
 
	friend vector2D operator*(float val, vector2D v)
	{
		v.x *= val;
		v.y *= val;

		return v;
	}
 
	friend vector2D operator/(float val, vector2D v)
	{
		v.x /= val;
		v.y /= val;

		return v;
	}

	float Magnitude()
	{
		return sqrt(x*x + y*y);
	}

	void Normalize()
	{
		float w = Magnitude();
		if (w < 0.00001) return;
		x /= w;
		y /= w;

	}
};
