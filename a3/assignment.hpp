#pragma once

#include "paths.hpp"

#include <atlas/glx/Context.hpp>
#include <atlas/glx/ErrorCallback.hpp>
#include <atlas/glx/GLSL.hpp>

#include <fmt/printf.h>

// Your code here.
#include <exception>
#include <iostream>
#include <string>

#include <atlas/glx/Buffer.hpp>
#include <atlas/utils/Cameras.hpp>
#include <atlas/utils/LoadObjFile.hpp>

#include <fmt/printf.h>
#include <magic_enum.hpp>

using namespace atlas;

static constexpr float nearVal{ 1.0f };
static constexpr float farVal{ 10000000000.0f };

static const std::vector<std::string> IncludeDir{ ShaderPath };

struct OpenGLError : std::runtime_error
{
    OpenGLError(const std::string& what_arg) : std::runtime_error(what_arg) {};
    OpenGLError(const char* what_arg) : std::runtime_error(what_arg) {};
};

class Cube;
class PointLight;
class Program;

class Pinhole 
{
public:
    Pinhole();
    Pinhole(glm::vec3 const& eye, glm::vec3 const& look, glm::vec3 const& up);
    void setCamera(glm::vec3 const& eye, glm::vec3 const& look, glm::vec3 const& up);
    void loadShaders();
    void loadCubeDataToGPU(Cube& cube);
    void reloadShaders();
    void render(int width, int height, Cube& cube, PointLight& light);
    void freeGPUData();

private:
    void setupUniformVariables();

    //camera vars
    glm::vec3 cEye;
    glm::vec3 cLook;
    glm::vec3 cUp;

    // Vertex buffers.
    GLuint mVao;
    GLuint mVbo;
    GLuint mEbo;

    // Shader data.
    GLuint mVertHandle;
    GLuint mFragHandle;
    GLuint mProgramHandle;
    glx::ShaderFile vertexSource;
    glx::ShaderFile fragmentSource;

    // Uniform variable data.
    GLuint mUniformModelLoc;
    GLuint mUniformViewLoc;
    GLuint mUniformProjectionLoc;
    GLuint mUniformCubeColour;
    GLuint mUniformLightColour;
    GLuint mUniformLightValue;
    GLuint mUniformLightPosition;

};

class Cube
{
public:
    Cube(std::array<float, 3 * 8> const& vs, std::array<unsigned int, 3 * 2 * 6> const& is, glm::vec3 const& c);
    void setVertices(std::array<float, 3 * 8> const& vs);
    void setIndices(std::array<unsigned int, 3 * 2 * 6> const& is);
    void setColour(glm::vec3 const& c);
    std::array<float, 3 * 8> getVertices();
    std::array<unsigned int, 3 * 2 * 6> getIndices();
    glm::vec3 getColour();

private:
    std::array<float, 3 * 8> mVertices; 
    std::array<unsigned int, 3 * 2 * 6> mIndices;
    glm::vec3 mColour;
};

class PointLight 
{
public:
    PointLight(std::array<float,3> const& point);
    void setPoint(std::array<float, 3> const& point);
    void scaleRadiance(float b);
    void setColour(glm::vec3 const& c);
    glm::vec3 getPoint();
    glm::vec3 getColour();
    glm::vec3 lValue();

private:
    glm::vec3 mColour;
    float mRadiance;
    std::array<float, 3> mPoint;
    glm::vec3 vPoint;
};


class Program
{
public:
    Program(int width, int height, std::string title);

    void run(Pinhole& ph, Cube& cube, PointLight& pLight);

    void freeGPUData();

private:
    static void errorCallback(int code, char const* message)
    {
        fmt::print("error ({}): {}\n", code, message);
    }

    void createGLContext();

    GLFWwindow* mWindow;
    glx::WindowSettings settings;
    glx::WindowCallbacks callbacks;
};
