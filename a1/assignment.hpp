#pragma once

#include <atlas/math/Math.hpp>
#include <atlas/math/Ray.hpp>

#include <fmt/printf.h>
#include <stb_image.h>
#include <stb_image_write.h>

#include <vector>
#include <stdlib.h>

using namespace atlas;
using Colour = math::Vector;


void saveToBMP(std::string const& filename,
               std::size_t width,
               std::size_t height,
               std::vector<Colour> const& image);

// Your code here.
struct ShadeRec {
	ShadeRec() : c{ 0, 0, 0 }, t{ 0 }
	{}

	ShadeRec(Colour m, float n) : c{ m }, t{ n }
	{}

	~ShadeRec() {};

	void setC(Colour m) { c = m; }
	void setT(float n) { t = n; }
	void closer(ShadeRec temp) {
		if (temp.dist() < t) {
			t = temp.dist();
			c = temp.colour();
		}
		return;
	}

	Colour colour() { return c; }
	float dist() const { return t; }

	Colour c;
	float t;
};

struct Ray
{
	Ray() : o{ 0,0,0 }, d{ 0,0,0 }
	{}

	Ray(math::Point m, math::Vector n) : o { m }, d { n }
	{}
	
	void setO(math::Point m) { o = m; }
	void setD(math::Vector n) { d = n; }
	
	math::Point getO() { return o; }
	math::Vector getD() { return d; }

	math::Point eval(float t)
	{
		return o + t * d;
	}

	math::Point o;
	math::Vector d;
};

class Sphere {
public:
	Sphere() : r{ 0 }, r2{ 0 }, centre{ 0,0,0 }
	{}

	Sphere(float m, math::Point n) : r{ m }, r2{ m * m }, centre{ n }
	{}

	bool hit(Ray ray, ShadeRec& sr)
	{
		math::Vector temp = ray.getO() - centre;
		float a = glm::dot(ray.getD(), ray.getD());
		float b = 2.0f * glm::dot(temp, ray.getD());
		float c = glm::dot(temp, temp) - (r2);
		float disc = b * b - 4.0f * a * c;

		if (disc < 0.0f) { return false; }
		else if (disc == 0.0f) {
			float t = (0.0f - b) / (2 * a);
			sr.setT(t);
		}
		else {
			float t1 = (0.0f - b + disc) / (2 * a);
			float t2 = (0.0f - b - disc) / (2 * a);

			if (t1 < t2) { sr.setT(t1); }
			else { sr.setT(t2); }
		}
		sr.setC(colour);
		return true;
	}

	void setColour(Colour m) { colour = m; }

private:
	float r;
	float r2;
	math::Point centre;
	Colour colour = { 0,0,0 };
};

class Plane {
public:
	Plane() : p{ 0,0,0 }, v{ 0,0,0 }
	{}

	Plane(math::Point m, math::Vector n) : p{ m }, v{ n }
	{}

	bool hit(Ray ray, ShadeRec& sr) {
		float a = v.x * (p.x - ray.getO().x);
		float b = v.y * (p.y - ray.getO().y);
		float c = v.z * (p.z - ray.getO().z);
		float d = glm::dot(v, ray.getD());

		if (d == 0) {
			return false;
		}
		else {
			sr.setT((a + b + c) / d);
			sr.setC(colour);
			return true;
		}
	}

	void setColour(Colour m) { colour = m; }

private:
	math::Point p;
	math::Vector v;
	Colour colour = { 0,0,0 };
};
