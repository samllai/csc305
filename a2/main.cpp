#include "assignment.hpp"
Shape::Shape() : mColour{0, 0, 0}
{}

void Shape::setColour(Colour const& col) {
    mColour = col;
}

Colour Shape::getColour() const {
    return mColour;
}

void Shape::setMaterial(std::shared_ptr<Material> const& material) {
    mMaterial = material;
}

std::shared_ptr<Material> Shape::getMaterial() const {
    return mMaterial;
}

// ***** Camera function members *****
Camera::Camera() :
    mEye{0.0f, 0.0f, 500.0f},
    mLookAt{0.0f},
    mUp{0.0f, 1.0f, 0.0f},
    mU{1.0f, 0.0f, 0.0f},
    mV{0.0f, 1.0f, 0.0f},
    mW{0.0f, 0.0f, 1.0f}
{}

void Camera::setEye(atlas::math::Point const& eye) {
    mEye = eye;
}

void Camera::setLookAt(atlas::math::Point const& lookAt) {
    mLookAt = lookAt;
}

void Camera::setUpVector(atlas::math::Vector const& up) {
    mUp = up;
}

void Camera::computeUVW() {
	mW = glm::normalize(mEye - mLookAt);
	mU = glm::normalize(glm::cross(mUp, mW));
	mV = glm::cross(mW, mU);
}

Pinhole::Pinhole() : Camera{}, mDistance{500.0f}, mZoom{1.0f}
{}

void Pinhole::setDistance(float dist) {
	mDistance = dist;
}

void Pinhole::setZoom(float zoom){
	mZoom = zoom;
}

atlas::math::Vector Pinhole::rayDirection(atlas::math::Point const& p) const {
	const auto dir = p.x*mU + p.y * mV - mDistance * mW;
	return glm::normalize(dir);
}

void Pinhole::renderScene(std::shared_ptr<World>& world) const {
	atlas::math::Ray<atlas::math::Vector> ray;
	ray.o = mEye;

	atlas::math::Point sample_point;
	atlas::math::Point pixel_point;

	float avg{ 1.0f / world->sampler->getNumSamples() };
	for (int i = 0; i < world->height; i++) {
		for (int j = 0; j < world->width; j++) {
			Colour avg_colour{ 0,0,0 };
			ShadeRec sr;
			sr.color = { 0,0,0 };
			sr.t = std::numeric_limits<float>::max();
			sr.world = world;

			for (int s = 0; s < world->sampler->getNumSamples(); s++) {
				sample_point = world->sampler->sampleUnitSquare();
				pixel_point.x = (j - 0.5f * (world->width - 1.0f));
				pixel_point.y = (i - 0.5f * (world->height - 1.0f));
				sample_point += pixel_point;
				ray.d = rayDirection(sample_point);
				bool hit{ false };

				for (auto const& object : world->scene) {
					hit |= object->hit(ray, sr);
				}
				if (hit) {
					avg_colour += sr.material->shade(sr);
				}
				else {
					avg_colour = world->background;
				}	
			}
			world->image.push_back({avg_colour.r * avg,
									avg_colour.g * avg,
									avg_colour.b * avg});
		}
	}
}

// ***** Sampler function members *****
Sampler::Sampler(int numSamples, int numSets) : mNumSamples{numSamples}, mNumSets{numSets}, mCount{0}, mJump{0}
{
    mSamples.reserve(mNumSets * mNumSamples);
    setupShuffledIndeces();
}

int Sampler::getNumSamples() const {
    return mNumSamples;
}

void Sampler::setupShuffledIndeces() {
    mShuffledIndeces.reserve(mNumSamples * mNumSets);
    std::vector<int> indices;

    std::random_device d;
    std::mt19937 generator(d());

    for (int j = 0; j < mNumSamples; ++j) {
        indices.push_back(j);
    }

    for (int p = 0; p < mNumSets; ++p) {
        std::shuffle(indices.begin(), indices.end(), generator);

        for (int j = 0; j < mNumSamples; ++j) {
            mShuffledIndeces.push_back(indices[j]);
        }
    }
}

atlas::math::Point Sampler::sampleUnitSquare() {
    if (mCount % mNumSamples == 0) {
        atlas::math::Random<int> engine;
        mJump = (engine.getRandomMax() % mNumSets) * mNumSamples;
    }

    return mSamples[mJump + mShuffledIndeces[mJump + mCount++ % mNumSamples]];
}

// ***** Light function members *****
Colour Light::L([[maybe_unused]] ShadeRec& sr) {
    return mRadiance * mColour;
}

void Light::scaleRadiance([[maybe_unused]] float b) {
	mRadiance = b;
}

void Light::setColour([[maybe_unused]] Colour const& c) {
	mColour = c;
}

Ambient::Ambient() : Light{}
{}

atlas::math::Vector Ambient::getDirection([[maybe_unused]] ShadeRec& sr) {
	return atlas::math::Vector{ 0,0,0 };
}

PointL::PointL() : Light{}
{}

PointL::PointL(atlas::math::Point const& p) : Light{} 
{
	setPoint(p);
}

void PointL::setPoint(atlas::math::Point const& p) {
	mPoint = p;
}

atlas::math::Vector PointL::getDirection(ShadeRec& sr) {
	atlas::math::Point pSurface = sr.ray.o + (sr.ray.d * sr.t);
	return glm::normalize(mPoint - pSurface);
}

Directional::Directional() : Light {}
{}

Directional::Directional(atlas::math::Vector const& d) : Light{}
{
	setDirection(d);
}

void Directional::setDirection(atlas::math::Vector const& d) {
	mDirection = glm::normalize(d);
}

atlas::math::Vector Directional::getDirection([[maybe_unused]] ShadeRec& sr) {
	return mDirection;
}

//Lambertian
Lambertian::Lambertian() : kD{}, cD{}
{}

Lambertian::Lambertian(float kd, Colour cd) : kD{kd}, cD{cd}
{}

Colour Lambertian::fn([[maybe_unused]] ShadeRec const& sr, [[maybe_unused]] atlas::math::Vector const& reflected, [[maybe_unused]] atlas::math::Vector const& incoming) const {
	return (kD * cD) * glm::one_over_pi<float>();
}

Colour Lambertian::rho([[maybe_unused]] ShadeRec const& sr, [[maybe_unused]] atlas::math::Vector const& reflected) const {
	return (kD * cD);
}

void Lambertian::setDiffuseReflection(float kd) {
	kD = kd;
}

void Lambertian::setDiffuseColour(Colour cd) {
	cD = cd;
}

//GlossySpec
GlossySpecular::GlossySpecular() : kS{}, cS{}, eP{}
{}

GlossySpecular::GlossySpecular(float ks, Colour cs, float exp) : kS{ ks }, cS{ cs }, eP{exp}
{}

Colour GlossySpecular::fn(ShadeRec const& sr, atlas::math::Vector const& reflected, atlas::math::Vector const& incoming)  const {
	const float kEpsilon{ 0.01f };
	float nWi{ glm::dot(sr.normal, incoming) };
	atlas::math::Vector r{ 2.0f * sr.normal * nWi -incoming };
	float rWo{ glm::dot(r, reflected) };

	if (rWo > kEpsilon) {
		return cS * kS * glm::pow<float, float>(rWo, eP);
	}
	else {
		return Colour{ 0,0,0 };
	}
}

Colour GlossySpecular::rho([[maybe_unused]] ShadeRec const& sr, [[maybe_unused]] atlas::math::Vector const& reflected) const {
	return Colour{ 0,0,0 };
}

void GlossySpecular::setSpecularReflection(float ks) {
	kS = ks;
}
void GlossySpecular::setSpecularColour(Colour cs) {
	cS = cs;
}
void GlossySpecular::setExp(float exp) {
	eP = exp;
}

//Materials
Matte::Matte() : Material{}, dif{Lambertian()}, amb{Lambertian()}
{}
	
Matte::Matte(float kA, float kD, Colour c) : Matte{}
{
	setDiffuseReflection(kD);
	setAmbientReflection(kA);
	setDiffuseColour(c);
}

void Matte::setDiffuseReflection(float dR) {
	dif.setDiffuseReflection(dR);
}

void Matte::setAmbientReflection(float aR) {
	amb.setDiffuseReflection(aR);
}

void Matte::setDiffuseColour(Colour c) {
	dif.setDiffuseColour(c);
	amb.setDiffuseColour(c);
}

Colour Matte::shade(ShadeRec& sr) {
	atlas::math::Vector wo{ -sr.ray.o };
	Colour ret{ amb.rho(sr, wo) * sr.world->ambient->L(sr) };
	size_t numL{ sr.world->lights.size() };
	
	for(size_t i{0}; i < numL; ++i) {
		atlas::math::Vector wi{ sr.world->lights[i]->getDirection(sr) };
		float nWi{ glm::dot(sr.normal, wi) };
		if (nWi > 0.0f) {
			ret += dif.fn(sr, wo, wi) * sr.world->lights[i]->L(sr) * nWi;
		}
	}
	return ret;
}

Phong::Phong() : Material{}, dif{Lambertian()}, amb{Lambertian()}, spec{GlossySpecular()}
{}

Phong::Phong(float kA, float kD, Colour c, float ks, Colour cs, float exp) : Phong{}
{
	setDiffuseReflection(kD);
	setAmbientReflection(kA);
	setDiffuseColour(c);
	setSpecularReflection(ks);
	setSpecularColour(cs);
	setExp(exp);
}

void Phong::setDiffuseReflection(float dR) {
	dif.setDiffuseReflection(dR);
}

void Phong::setAmbientReflection(float aR) {
	amb.setDiffuseReflection(aR);
}

void Phong::setDiffuseColour(Colour c) {
	dif.setDiffuseColour(c);
	amb.setDiffuseColour(c);
}

void Phong::setSpecularReflection(float ks) {
	spec.setSpecularReflection(ks);
}

void Phong::setSpecularColour(Colour cs) {
	spec.setSpecularColour(cs);
}

void Phong::setExp(float exp) {
	spec.setExp(exp);
}

Colour Phong::shade(ShadeRec& sr) {
	atlas::math::Vector wo{ -sr.ray.d };
	Colour ret{ amb.rho(sr, wo) * sr.world->ambient->L(sr) };
	size_t numL{ sr.world->lights.size() };

	for (size_t i{ 0 }; i < numL; ++i) {
		atlas::math::Vector wi{ sr.world->lights[i]->getDirection(sr) };
		float nWi{ glm::dot(sr.normal, wi) };
		if (nWi > 0.0f) {
			ret += (dif.fn(sr, wo, wi) + spec.fn(sr, wo, wi)) * sr.world->lights[i]->L(sr) * nWi;
		}
	}
	return ret;
}

// ***** Sphere function members *****
Sphere::Sphere(atlas::math::Point center, float radius) : mCentre{center}, mRadius{radius}, mRadiusSqr{radius * radius}
{}

bool Sphere::hit(atlas::math::Ray<atlas::math::Vector> const& ray, ShadeRec& sr) const {
	atlas::math::Vector tmp = ray.o - mCentre;
	float t{ std::numeric_limits<float>::max() };
    bool intersect{intersectRay(ray, t)};

    // update ShadeRec info about new closest hit
    if (intersect && t < sr.t) {
        sr.normal   = (tmp + t * ray.d) / mRadius;
        sr.ray      = ray;
        sr.color    = mColour;
        sr.t        = t;
        sr.material = mMaterial;
    }
    return intersect;
}

bool Sphere::intersectRay(atlas::math::Ray<atlas::math::Vector> const& ray, float& tMin) const {
    const auto tmp{ray.o - mCentre};
    const auto a{glm::dot(ray.d, ray.d)};
    const auto b{2.0f * glm::dot(ray.d, tmp)};
    const auto c{glm::dot(tmp, tmp) - mRadiusSqr};
    const auto disc{(b * b) - (4.0f * a * c)};

    if (atlas::core::geq(disc, 0.0f)) {
        const float kEpsilon{0.01f};
        const float e{std::sqrt(disc)};
        const float denom{2.0f * a};

        // Look at the negative root first
        float t = (-b - e) / denom;
        if (atlas::core::geq(t, kEpsilon)) {
            tMin = t;
            return true;
        }

        // Now the positive root
        t = (-b + e);
        if (atlas::core::geq(t, kEpsilon)) {
            tMin = t;
            return true;
        }
    }
    return false;
}


Triangle::Triangle() : p1{ 0,0,0 }, p2{ 0,0,1 }, p3{ 1,0,0 }, mNormal{0,1,0}
{}

Triangle::Triangle(atlas::math::Point i, atlas::math::Point j, atlas::math::Point k) : p1{ i }, p2{ j }, p3{ k }
{
	atlas::math::Vector n = glm::cross((p2 - p1), (p3 - p1));
	mNormal = glm::normalize(n);
}

bool Triangle::hit(atlas::math::Ray<atlas::math::Vector> const& ray, ShadeRec& sr) const {
	float t{ std::numeric_limits<float>::max() };
	bool intersect{ intersectRay(ray, t) };

	// update ShadeRec info about new closest hit
	if (intersect && t < sr.t) {
		sr.normal = mNormal;
		sr.ray = ray;
		sr.color = mColour;
		sr.t = t;
		sr.material = mMaterial;
	}
	return intersect;
}

bool Triangle::intersectRay(atlas::math::Ray<atlas::math::Vector> const& ray, float& tMin) const {
	float a{ p1.x - p2.x };
	float b{ p1.x - p3.x };
	float c{ ray.d.x };
	float d{ p1.x - ray.o.x };
	float e{ p1.y - p2.y };
	float f{ p1.y - p3.y };
	float g{ ray.d.y };
	float h{ p1.y - ray.o.y };
	float i{ p1.z - p2.z };
	float j{ p1.z - p3.z };
	float k{ ray.d.z };
	float l{ p1.z - ray.o.z };
	float m{ f * k - g * j };
	float n{ h * k - g * l };
	float p{ f * l - h * j };
	float q{ g * i - e * k };
	float s{ e * j - f * i };
	float inv_denom{ 1.0f / (a * m + b * q + c * s) };
	float e1{ d * m - b * n - c * p };
	float beta{ e1 * inv_denom };

	if (beta < 0.0f) { return false; }

	float r{ e * l - h * i };
	float e2{ a * n + d * q + c * r };
	float gamma{ e2 * inv_denom };

	if (gamma < 0.0f) { return false; }

	if (beta + gamma > 1.0f) { return false; }

	float e3{ a * p - b * r + d * s };
	float t{ e3 * inv_denom };

	if (t < 0.0f) { return false; }
	tMin = t;
	return true;
}

Plane::Plane(atlas::math::Point p, atlas::math::Vector n) : mPoint{p}, mNormal{n}
{}

bool Plane::hit(atlas::math::Ray<atlas::math::Vector> const& ray, ShadeRec& sr) const {
	float t{ std::numeric_limits<float>::max() };
	bool intersect{ intersectRay(ray, t) };

	// update ShadeRec info about new closest hit
	if (intersect && t < sr.t) {
		sr.normal = mNormal;
		sr.ray = ray;
		sr.color = mColour;
		sr.t = t;
		sr.material = mMaterial;
	}
	return intersect;
}

bool Plane::intersectRay(atlas::math::Ray<atlas::math::Vector> const& ray, float& tMin) const {
	const float kEpsilon{ 0.01f };
	float t{ glm::dot((mPoint - ray.o), mNormal) / glm::dot(ray.d, mNormal) };
	if (t > kEpsilon) {
		tMin = t;
		return true;
	}
	else {
		return false;
	}
}

// ***** Random function members *****
Random::Random(int numSamples, int numSets) : Sampler{numSamples, numSets}
{
    generateSamples();
}

void Random::generateSamples() {
    atlas::math::Random<float> engine;
    for (int p = 0; p < mNumSets; ++p)
    {
        for (int q = 0; q < mNumSamples; ++q)
        {
            mSamples.push_back(atlas::math::Point{engine.getRandomOne(), engine.getRandomOne(), 0.0f});
        }
    }
}

//Check gamut
void checkGamut(std::vector<Colour>& px) {
	bool over {false};
	float max {1.0f};
	for (int i = 0; i < px.size(); i++){
		if (px[i].r > max) {
			over = true;
			max = px[i].r;
		} else if (px[i].g > max) {
			over = true;
			max = px[i].g;
		} else if (px[i].b > max) {
			over = true;
			max = px[i].b;
		}
	}
	if (over) {
		for (int k = 0; k < px.size(); k++) {
			px[k].r /= max;
			px[k].g /= max;
			px[k].b /= max;
		}
	}
}

int main() {
    // Your code here.
	//initialize camera
	Pinhole phcam;
	phcam.setDistance(1000.0f);
	phcam.setEye(atlas::math::Point{ 0,0, -1000.0f });
	phcam.setLookAt(atlas::math::Point{ 40.0f, -25.0f, 0 });
	phcam.setUpVector(atlas::math::Vector{ -3.0f, -2.0f, -1.0f });
	phcam.setZoom(3.0f);
	phcam.computeUVW();


	//initialize world
	std::shared_ptr<World> w{ std::make_shared<World>() };
	w->width = 600;
	w->height = 600;
	w->background = { 0,0,0 };
	w->sampler = std::make_shared<Random>(16, 16);
	
	//build scene
	w->scene.push_back(std::make_shared<Plane>(atlas::math::Point{ 0,-1500.0f,0 }, atlas::math::Vector{ 0,1.0f,0 }));
	w->scene[0]->setColour({ 0,1,0 });
	w->scene[0]->setMaterial(std::make_shared<Matte>(0.8f, 0.8f, Colour{ 0,1,0 }));
	
	w->scene.push_back(std::make_shared<Triangle>(atlas::math::Point{ 0,10.0f,30.0f }, atlas::math::Point{ -200.0f, -100.0f, 60.0f }, atlas::math::Point{-100.0f, 250.0f, 30.0f}));
	w->scene[1]->setColour({ 0.8f, 0, 0.5f });
	w->scene[1]->setMaterial(std::make_shared<Matte>(0.5f, 0.3f, Colour{ 0.8f, 0, 0.5f }));

	w->scene.push_back(std::make_shared<Sphere>(atlas::math::Point{ 100.0f, -50.0f, 30.0f }, 150.0f));
	w->scene[2]->setColour({ 0.2f, 0.1f, 0.7f });
	w->scene[2]->setMaterial(std::make_shared<Phong>(0.2f, 0.5f, Colour{ 0.2f, 0.1f, 0.7f }, 0.2f, Colour{ 0.2f, 0.1f, 0.7f}, 3.0f));

	//build lights
	w->ambient = std::make_shared<Ambient>();
	w->lights.push_back(std::make_shared<PointL>(atlas::math::Point{ -300.0f, 400.0f, -1000.0f }));
	w->lights[0]->setColour({ 1, 1, 1 });
	w->lights[0]->scaleRadiance(50.0f);
	
	phcam.renderScene(w);
	saveToBMP("noGamut.bmp", 600, 600, w->image);
	checkGamut(w->image);
	saveToBMP("pointScene.bmp", 600, 600, w->image);

	
	w->image = std::vector<Colour>{};
	w->lights[0] = std::make_shared<Directional>(atlas::math::Vector{3.0f, 2.0f, -3.0f});
	w->lights[0]->setColour({ 1, 1, 1 });
	w->lights[0]->scaleRadiance(10.0f);

	phcam.renderScene(w);
	checkGamut(w->image);
	saveToBMP("directionalScene.bmp", 600, 600, w->image);

    return 0;
}

/**
 * Saves a BMP image file based on the given array of pixels. All pixel values
 * have to be in the range [0, 1].
 *
 * @param filename The name of the file to save to.
 * @param width The width of the image.
 * @param height The height of the image.
 * @param image The array of pixels representing the image.
 */
void saveToBMP(std::string const& filename,
               std::size_t width,
               std::size_t height,
               std::vector<Colour> const& image)
{
    std::vector<unsigned char> data(image.size() * 3);

    for (std::size_t i{0}, k{0}; i < image.size(); ++i, k += 3)
    {
        Colour pixel = image[i];
        data[k + 0]  = static_cast<unsigned char>(pixel.r * 255);
        data[k + 1]  = static_cast<unsigned char>(pixel.g * 255);
        data[k + 2]  = static_cast<unsigned char>(pixel.b * 255);
    }

    stbi_write_bmp(filename.c_str(),
                   static_cast<int>(width),
                   static_cast<int>(height),
                   3,
                   data.data());
}
