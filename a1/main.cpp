#include "assignment.hpp"

int main()
{
    // Your code here.

	std::vector <Colour> image{ 600 * 600 };

	const std::size_t imageWidth{ 600 };
	const std::size_t imageHeight{ 600 };

	Ray ray;
	ray.setD({ 0, 0, -1 });
	
	Plane p1({ -2,3,-1 }, { 1,2,3 });
	p1.setColour({ 0,1,0.85f });
	Plane p2({ 50,-20,40 }, { -3,5,2 });
	p2.setColour({ 0.9f,0,0.4f });
	
	Sphere s1(220.0f, { 300,-400,-125 });
	s1.setColour({ 0.8f,0.8f,0.1f });
	Sphere s2(400.0f, { 500,250,3 });
	s2.setColour({ 0.9f,0.5f,0 });
	Sphere s3(250.0f, { -40,500,30 });
	s3.setColour({ 0.7f,0.2f,0.1f });
	Sphere s4(270.0f, { -50,-300,15 });
	s4.setColour({ 0.2f,0.5f,0.6f });
	Sphere s5(120.0f, { 0,0,40 });
	s5.setColour({ 0.5f,0.2f,0.1f });
	Sphere s6(300.0f, { -400,300,-30 });
	s6.setColour({ 1,0,0 });


	ShadeRec closest;
	Colour c = { 0,0,0 };
	float off1 = 0.0f;
	float off2 = 0.0f;
	
	for (std::size_t y{ 0 }; y < imageHeight; ++y)
	{
		for (std::size_t x{ 0 }; x < imageWidth; ++x)
		{			
			ShadeRec srt;
			c = { 0,0,0 };
			closest.setC(c);
			closest.setT(1000.0f);

			for (int i{ 0 }; i < 4; i++) {
				float originX = (x - 0.5f * (imageWidth - 1.0f));
				float originY = (y - 0.5f * (imageHeight - 1.0f));

				off1 = (float)(rand() % 6);
				if (i%2) { off1 *= -1.0f; }
				off1 /= 10.0f;
			
				off2 = (float)(rand() % 6);
				if (i>1) { off2 *= -1.0f; }
				off2 /= 10.0f;
			
				originX += off1;
				originY = (originY + off2) * (-1.0f);
			
				ray.setO({ originX, originY, 100.0f });
			
				if (p1.hit(ray, srt)) {	closest.closer(srt); }

				if (p2.hit(ray, srt)) {	closest.closer(srt); }
			
				if (s1.hit(ray, srt)) {	closest.closer(srt); }

				if (s2.hit(ray, srt)) {	closest.closer(srt); }

				if (s3.hit(ray, srt)) {	closest.closer(srt); }

				if (s4.hit(ray, srt)) {	closest.closer(srt); }

				if (s5.hit(ray, srt)) {	closest.closer(srt); }

				if (s6.hit(ray, srt)) { closest.closer(srt); }
			}
			image[x + y * imageHeight] = closest.colour();
		}
	}

	saveToBMP("a1Jittered.bmp", 600, 600, image);
	
	for (std::size_t y{ 0 }; y < imageHeight; ++y)
	{
		for (std::size_t x{ 0 }; x < imageWidth; ++x)
		{
			ShadeRec srt;
			c = { 0,0,0 };
			closest.setC(c);
			closest.setT(1000.0f);

			for (int i{ 0 }; i < 4; i++) {
				float originX = (x - 0.5f * (imageWidth - 1.0f));
				float originY = (y - 0.5f * (imageHeight - 1.0f));

				off1 = (float)(rand() % 6);
				if (rand() % 2) { off1 *= -1.0f; }
				off1 /= 10.0f;

				off2 = (float)(rand() % 6);
				if (rand() % 2) { off2 *= -1.0f; }
				off2 /= 10.0f;

				originX += off1;
				originY = (originY + off2) * (-1.0f);

				ray.setO({ originX, originY, 100.0f });

				if (p1.hit(ray, srt)) { closest.closer(srt); }

				if (p2.hit(ray, srt)) { closest.closer(srt); }

				if (s1.hit(ray, srt)) { closest.closer(srt); }

				if (s2.hit(ray, srt)) { closest.closer(srt); }

				if (s3.hit(ray, srt)) { closest.closer(srt); }

				if (s4.hit(ray, srt)) { closest.closer(srt); }

				if (s5.hit(ray, srt)) { closest.closer(srt); }

				if (s6.hit(ray, srt)) { closest.closer(srt); }
			}
			image[x + y * imageHeight] = closest.colour();
		}
	}
	saveToBMP("a1Random.bmp", 600, 600, image);
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
