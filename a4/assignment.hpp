#pragma once

#include <atlas/math/Math.hpp>
#include <atlas/math/Ray.hpp>

#include <fmt/printf.h>
#include <stb_image.h>
#include <stb_image_write.h>

#include <vector>
#include <atlas/core/Float.hpp>
#include <atlas/math/Random.hpp>
#include <limits>

//using namespace atlas;
using Colour = atlas::math::Vector;


void saveToBMP(std::string const& filename,
               std::size_t width,
               std::size_t height,
               std::vector<Colour> const& image);

// Your code here.
class Tracer;
class Camera;
class Shape;
class Sampler;
class BRDF;
class Material;
class Light;

struct World {
	std::size_t width, height;
	Colour background;
	std::shared_ptr<Sampler> sampler;
	std::vector<std::shared_ptr<Shape>> scene;
	std::vector<Colour> image;
	std::vector<std::shared_ptr<Light>> lights;
	std::shared_ptr<Light> ambient;
	int max_depth;
};

struct ShadeRec {
	bool hit;
	Colour color;
	float t;
	atlas::math::Normal normal;
	atlas::math::Ray<atlas::math::Vector> ray;
	atlas::math::Point hitPoint;
	std::shared_ptr<Material> material;
	std::shared_ptr<World> world;
	int depth;
};

class Tracer{
public:
	Tracer();
	Tracer(std::shared_ptr<World> w);
	virtual Colour trace_ray([[maybe_unused]] atlas::math::Ray<atlas::math::Vector> const& ray, [[maybe_unused]] float & tMin, [[maybe_unused]] int const depth) const = 0;
protected:
	std::shared_ptr<World> mWorld;
};

class RayCast : public Tracer {
public:
	RayCast();
	RayCast(std::shared_ptr<World> w);
	Colour trace_ray([[maybe_unused]] atlas::math::Ray<atlas::math::Vector> const& ray, [[maybe_unused]] float& tMin, [[maybe_unused]] int const depth) const;
};

// Abstract classes defining the interfaces for concrete entities

class Camera {
public:
	Camera();
	virtual ~Camera() = default;
	virtual void renderScene(std::shared_ptr<World>& world) const = 0;
	void checkGamut(std::vector<Colour>& px) const;
	void setEye(atlas::math::Point const& eye);
	void setLookAt(atlas::math::Point const& lookAt);
	void setUpVector(atlas::math::Vector const& up);
	void computeUVW();

protected:
	atlas::math::Point mEye;
	atlas::math::Point mLookAt;
	atlas::math::Point mUp;
	atlas::math::Vector mU, mV, mW;
};

class Pinhole : public Camera {
public:
	Pinhole();
	void renderScene(std::shared_ptr<World>& world) const;
	void setDistance(float dist);
	void setZoom(float zoom);
	atlas::math::Vector rayDirection(atlas::math::Point const& p) const;

protected:
	float mDistance;
	float mZoom;
};

class Sampler {
public:
	Sampler(int numSamples, int numSets);
	virtual ~Sampler() = default;
	int getNumSamples() const;
	void setupShuffledIndeces();
	virtual void generateSamples() = 0;
	atlas::math::Point sampleUnitSquare();
	void mapToHemisphere(float const e);
	atlas::math::Point sampleHemisphere();

protected:
	std::vector<atlas::math::Point> mSamples;
	std::vector<atlas::math::Point> mHSamples;
	std::vector<int> mShuffledIndeces;
	int mNumSamples;
	int mNumSets;
	unsigned long mCount;
	int mJump;
};

class BRDF {
public:
	virtual ~BRDF() = default;
	virtual Colour fn(ShadeRec const& sr,
		atlas::math::Vector const& reflected,
		atlas::math::Vector const& incoming) const = 0;
	virtual Colour rho(ShadeRec const& sr,
		atlas::math::Vector const& reflected) const = 0;
};

class Lambertian : public BRDF {
public:
	Lambertian();
	Lambertian(float kd, Colour cd);
	Colour fn([[maybe_unused]] ShadeRec const& sr, [[maybe_unused]] atlas::math::Vector const& reflected, [[maybe_unused]] atlas::math::Vector const& incoming) const;
	Colour rho([[maybe_unused]] ShadeRec const& sr, [[maybe_unused]] atlas::math::Vector const& reflected) const;
	void setDiffuseReflection(float kd);
	void setDiffuseColour(Colour cd);

protected:
	float kD;
	Colour cD;
};

class FresnelReflector : public BRDF {
public:
	FresnelReflector();
	Colour fn([[maybe_unused]] ShadeRec const& sr, [[maybe_unused]] atlas::math::Vector const& reflected, [[maybe_unused]] atlas::math::Vector const& incoming) const;
	Colour rho([[maybe_unused]] ShadeRec const& sr, [[maybe_unused]] atlas::math::Vector const& reflected) const;
	Colour sample_f([[maybe_unused]] ShadeRec const& sr, [[maybe_unused]] atlas::math::Vector const& reflected, [[maybe_unused]] atlas::math::Vector & incoming) const;
	float fresnel(ShadeRec const& sr) const;

protected:
	Colour cR;
	float eta_out;
	float eta_in;
};

class FresnelTransmitter : public BRDF {
public:
	FresnelTransmitter();
	Colour fn([[maybe_unused]] ShadeRec const& sr, [[maybe_unused]] atlas::math::Vector const& reflected, [[maybe_unused]] atlas::math::Vector const& incoming) const;
	Colour rho([[maybe_unused]] ShadeRec const& sr, [[maybe_unused]] atlas::math::Vector const& reflected) const;
};

class GlossySpecular : public BRDF {
public:
	GlossySpecular();
	GlossySpecular(float ks, Colour cs, float exp);
	Colour fn(ShadeRec const& sr, atlas::math::Vector const& reflected, atlas::math::Vector const& incoming)  const;
	Colour rho([[maybe_unused]] ShadeRec const& sr, [[maybe_unused]] atlas::math::Vector const& reflected) const;
	void setSpecularReflection(float ks);
	void setSpecularColour(Colour cs);
	void setExp(float exp);

protected:
	float kS;
	Colour cS;
	float eP;
};

class Material {
public:
	virtual ~Material() = default;
	virtual Colour shade(ShadeRec& sr) = 0;
};

class Matte : public Material {
public:
	Matte();
	Matte(float kA, float kD, Colour c);
	void setDiffuseReflection(float dR);
	void setAmbientReflection(float aR);
	void setDiffuseColour(Colour dC);
	Colour shade(ShadeRec& sr);

protected:
	Lambertian amb;
	Lambertian dif;
};

class Phong : public Material {
public:
	Phong();
	Phong(float kA, float kD, Colour c, float ks, Colour cs, float exp);
	void setDiffuseReflection(float dR);
	void setAmbientReflection(float aR);
	void setDiffuseColour(Colour c);
	void setSpecularReflection(float ks);
	void setSpecularColour(Colour cs);
	void setExp(float exp);
	Colour shade(ShadeRec& sr);

protected:
	Lambertian amb;
	Lambertian dif;
	GlossySpecular spec;
};

class Dielectric : public Phong {
public:
	Dielectric();
	Colour shade([[maybe_unused]]ShadeRec& sr);

protected:
	Colour cf_in;
	Colour cf_out;
	FresnelReflector fresnel_brdf;
	FresnelTransmitter fresnel_btdf;
};

class Light {
public:
	virtual atlas::math::Vector getDirection(ShadeRec& sr) = 0;
	virtual Colour L(ShadeRec& sr);
	void scaleRadiance(float b);
	void setColour(Colour const& c);
	void setShadow(bool s);
	bool castShadows();
	virtual bool inShadow(atlas::math::Ray<atlas::math::Vector>& shadowRay, ShadeRec& sr) const = 0;

protected:
	Colour mColour;
	float mRadiance;
	bool shadow;
};

class Ambient : public Light {
public:
	Ambient();
	atlas::math::Vector getDirection([[maybe_unused]] ShadeRec& sr);
	bool inShadow([[maybe_unused]] atlas::math::Ray<atlas::math::Vector>& shadowRay, [[maybe_unused]] ShadeRec& sr) const;

private:
	atlas::math::Vector mDirection;
};

class PointL : public Light {
public:
	PointL();
	PointL(atlas::math::Point const& p);
	void setPoint(atlas::math::Point const& p);
	atlas::math::Vector getDirection(ShadeRec& sr);
	bool inShadow(atlas::math::Ray<atlas::math::Vector>& shadowRay, ShadeRec& sr) const;

private:
	atlas::math::Point mPoint;
	atlas::math::Vector mDirection;
};

class Directional : public Light {
public:
	Directional();
	Directional(atlas::math::Vector const& d);
	void setDirection(atlas::math::Vector const& d);
	atlas::math::Vector getDirection(ShadeRec& sr);

private:
	atlas::math::Vector mDirection;
};

class AmbientOccluder : public Light {
public:
	AmbientOccluder();
	atlas::math::Vector getDirection([[maybe_unused]] ShadeRec& sr);
	bool inShadow([[maybe_unused]] atlas::math::Ray<atlas::math::Vector>& shadowRay, [[maybe_unused]] ShadeRec& sr) const;
	void setSampler(Sampler* s_ptr);
	Colour L(ShadeRec& sr);
	void setMin(float m);

protected:
	atlas::math::Vector mU;
	atlas::math::Vector mV;
	atlas::math::Vector mW;
	Sampler* mSampler;
	Colour mMin;
};

class Shape {
public:
	Shape();
	virtual ~Shape() = default;

	// if t computed is less than the t in sr, it and the color should be updated in sr
	virtual bool hit(atlas::math::Ray<atlas::math::Vector> const& ray, ShadeRec& sr) const = 0;
	void setColour(Colour const& col);
	Colour getColour() const;
	void setMaterial(std::shared_ptr<Material> const& material);
	std::shared_ptr<Material> getMaterial() const;
	virtual bool shadowHit(atlas::math::Ray<atlas::math::Vector>& shadowRay, float& tmin) const = 0;
	void setEpsilon(float kE);

protected:
	virtual bool intersectRay(atlas::math::Ray<atlas::math::Vector> const& ray, float& tMin) const = 0;
	Colour mColour;
	std::shared_ptr<Material> mMaterial;
	float kEpsilon;
};

// Concrete classes which we can construct and use in our ray tracer

class Sphere : public Shape {
public:
	Sphere(atlas::math::Point center, float radius);
	bool hit(atlas::math::Ray<atlas::math::Vector> const& ray, ShadeRec& sr) const;
	bool shadowHit(atlas::math::Ray<atlas::math::Vector>& shadowRay, float& tmin) const;

private:
	bool intersectRay(atlas::math::Ray<atlas::math::Vector> const& ray, float& tMin) const;
	atlas::math::Point mCentre;
	float mRadius;
	float mRadiusSqr;
};

class Triangle : public Shape {
public:
	Triangle();
	Triangle(atlas::math::Point i, atlas::math::Point j, atlas::math::Point k);
	bool hit(atlas::math::Ray<atlas::math::Vector> const& ray, ShadeRec& sr) const;
	bool shadowHit(atlas::math::Ray<atlas::math::Vector>& shadowRay, float& tmin) const;

private:
	bool intersectRay(atlas::math::Ray<atlas::math::Vector> const& ray, float& tMin) const;
	atlas::math::Point p1;
	atlas::math::Point p2;
	atlas::math::Point p3;
	atlas::math::Vector mNormal;
};

class Plane : public Shape {
public:
	Plane(atlas::math::Point mPoint, atlas::math::Vector mNormal);
	bool hit(atlas::math::Ray<atlas::math::Vector> const& ray, ShadeRec& sr) const;
	bool shadowHit(atlas::math::Ray<atlas::math::Vector>& shadowRay, float& tmin) const;

private:
	bool intersectRay(atlas::math::Ray<atlas::math::Vector> const& ray, float& tMin) const;
	atlas::math::Point mPoint;
	atlas::math::Vector mNormal;
};

class Random : public Sampler
{
public:
	Random(int numSamples, int numSets);

	void generateSamples();
};