#include "assignment.hpp"
Pinhole::Pinhole()
{
    // allocate the memory to hold the program and shader data
    mProgramHandle = glCreateProgram();
    mVertHandle = glCreateShader(GL_VERTEX_SHADER);
    mFragHandle = glCreateShader(GL_FRAGMENT_SHADER);

    //default camera settings
    cEye = glm::vec3(1.0f, 1.0f, 1.0f);
    cLook = glm::vec3(0.0f, 0.0f, 0.0f);
    cUp = glm::vec3(0.0f, 1.0f, 0.0f);
}

Pinhole::Pinhole(glm::vec3 const& eye, glm::vec3 const& look, glm::vec3 const& up)
{
    // allocate the memory to hold the program and shader data
    mProgramHandle = glCreateProgram();
    mVertHandle = glCreateShader(GL_VERTEX_SHADER);
    mFragHandle = glCreateShader(GL_FRAGMENT_SHADER);
    cEye = eye;
    cLook = look;
    cUp = up;
}

void Pinhole::setCamera(glm::vec3 const& eye, glm::vec3 const& look, glm::vec3 const& up){
    cEye = eye;
    cLook = look;
    cUp = up;
}

void Pinhole::loadShaders(){
    std::string shaderRoot{ ShaderPath };
    vertexSource =
        glx::readShaderSource(shaderRoot + "cube.vert", IncludeDir);
    fragmentSource =
        glx::readShaderSource(shaderRoot + "cube.frag", IncludeDir);

    if (auto result{ glx::compileShader(vertexSource.sourceString, mVertHandle) };
        result)
    {
        throw OpenGLError(*result);
    }

    if (auto result =
        glx::compileShader(fragmentSource.sourceString, mFragHandle);
        result)
    {
        throw OpenGLError(*result);
    }

    // communicate to OpenGL the shaders used to render the Cube
    glAttachShader(mProgramHandle, mVertHandle);
    glAttachShader(mProgramHandle, mFragHandle);

    if (auto result = glx::linkShaders(mProgramHandle); result)
    {
        throw OpenGLError(*result);
    }

    setupUniformVariables();
}

void Pinhole::loadCubeDataToGPU(Cube& cube){
    std::array<float, 3 * 8> vertices = cube.getVertices();
    std::array<unsigned int, 3 * 2 * 6> indices = cube.getIndices();
    // create buffer to hold triangle vertex data
    glCreateBuffers(1, &mVbo);
    // allocate and initialize buffer to vertex data
    glNamedBufferStorage(mVbo, glx::size<float>(vertices.size()), vertices.data(), 0);

    // create holder for all buffers
    glCreateVertexArrays(1, &mVao);
    // bind vertex buffer to the vertex array
    glVertexArrayVertexBuffer(mVao, 0, mVbo, 0, glx::stride<float>(3));

    // enable attributes for the two components of a vertex
    glEnableVertexArrayAttrib(mVao, 0);

    // specify to OpenGL how the vertices and colors are laid out in the buffer
    glVertexArrayAttribFormat(mVao, 0, 3, GL_FLOAT, GL_FALSE, glx::relativeOffset<float>(0));

    // associate the vertex attributes (coordinates and color) to the vertex
    // attribute
    glVertexArrayAttribBinding(mVao, 0, 0);

    glCreateBuffers(1, &mEbo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mEbo);
    glNamedBufferStorage(mEbo, glx::size<unsigned int>(indices.size()), indices.data(), 0);

}

void Pinhole::reloadShaders(){
    if (glx::shouldShaderBeReloaded(vertexSource))
    {
        glx::reloadShader(
            mProgramHandle, mVertHandle, vertexSource, IncludeDir);
    }

    if (glx::shouldShaderBeReloaded(fragmentSource))
    {
        glx::reloadShader(
            mProgramHandle, mFragHandle, fragmentSource, IncludeDir);
    }
}

void Pinhole::render([[maybe_unused]] int width, [[maybe_unused]] int height, Cube& cube, PointLight& light){
    reloadShaders();

    auto projMat{ glm::perspective(static_cast<float>(glm::radians(45.0f)), static_cast<float>(width) / (height), nearVal, farVal) };
    auto viewMat{ glm::lookAt(cEye, cLook, cUp) };

    //float position = static_cast<float>(glfwGetTime()) * 30.0f;
    auto modelMat{ glm::rotate(math::Matrix4(1.0f), glm::radians(40.0f), glm::vec3(0.0f,1.0f,0.0f)) };

    glm::vec3 cubeColour{ cube.getColour() };
    glm::vec3 lightColour{ light.getColour() };
    glm::vec3 lightValue{ light.lValue() };
    glm::vec3 lightPos{ light.getPoint() };

    // tell OpenGL which program object to use to render the Cube
    glUseProgram(mProgramHandle);

    //bind uniform variables
    glUniformMatrix4fv(mUniformProjectionLoc, 1, GL_FALSE, glm::value_ptr(projMat));
    glUniformMatrix4fv(mUniformViewLoc, 1, GL_FALSE, glm::value_ptr(viewMat));
    glUniformMatrix4fv(mUniformModelLoc, 1, GL_FALSE, glm::value_ptr(modelMat));
    glUniform3fv(mUniformCubeColour, 1, glm::value_ptr(cubeColour));
    glUniform3fv(mUniformLightColour, 1, glm::value_ptr(lightColour));
    glUniform3fv(mUniformLightValue, 1, glm::value_ptr(lightValue));
    glUniform3fv(mUniformLightPosition, 1, glm::value_ptr(lightPos));

    // tell OpenGL which vertex array object to use to render the Triangle
    glBindVertexArray(mVao);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mEbo);
    // actually render the Triangle
    //glDrawArrays(GL_TRIANGLES, 0, 3 * 2 * 6);
    glDrawElements(GL_TRIANGLES, 3 * 2 * 6, GL_UNSIGNED_INT, 0);
}

void Pinhole::freeGPUData(){
    // unwind all the allocations made
    glDeleteVertexArrays(1, &mVao);
    glDeleteBuffers(1, &mVbo);
    glDeleteBuffers(1, &mEbo);
    glDeleteShader(mFragHandle);
    glDeleteShader(mVertHandle);
    glDeleteProgram(mProgramHandle);
}

void Pinhole::setupUniformVariables(){
    mUniformModelLoc = glGetUniformLocation(mProgramHandle, "model");
    mUniformViewLoc = glGetUniformLocation(mProgramHandle, "view");
    mUniformProjectionLoc = glGetUniformLocation(mProgramHandle, "projection");
    mUniformCubeColour = glGetUniformLocation(mProgramHandle, "ccolour");
    mUniformLightColour = glGetUniformLocation(mProgramHandle, "lcolour");
    mUniformLightValue = glGetUniformLocation(mProgramHandle, "lvalue");
    mUniformLightPosition = glGetUniformLocation(mProgramHandle, "lpos");
}


// ===---------------CUBE-----------------===

Cube::Cube(std::array<float, 3 * 8> const& vs, std::array<unsigned int, 3 * 2 * 6> const& is, glm::vec3 const& c)
{
    mVertices = vs;
    mIndices = is;
    mColour = c;
}

void Cube::setVertices(std::array<float, 3 * 8> const& vs){
    mVertices = vs;
}

void Cube::setIndices(std::array<unsigned int, 3 * 2 * 6> const& is) {
    mIndices = is;
}

void Cube::setColour(glm::vec3 const& c) {
    mColour = c;
}

std::array<float, 3 * 8> Cube::getVertices() {
    return mVertices;
}

std::array<unsigned int, 3 * 2 * 6> Cube::getIndices() {
    return mIndices;
}

glm::vec3 Cube::getColour() {
    return mColour;
}

//Point Light

PointLight::PointLight(std::array<float, 3> const& point)
{
    mPoint = point;
    vPoint = glm::vec3(point[0], point[1], point[2]);
    mRadiance = 0.5f;
    mColour = { 1.0f, 1.0f, 1.0f };
}

void PointLight::setPoint(std::array<float, 3> const& point){
    mPoint = point;
    vPoint = glm::vec3(point[0], point[1], point[2]);
}
void PointLight::scaleRadiance(float b){
    mRadiance = b;
}

void PointLight::setColour(glm::vec3 const& c) {
    mColour = c;
}

glm::vec3 PointLight::getPoint(){
    return vPoint;
}

glm::vec3 PointLight::getColour(){
    return mColour;
}

glm::vec3 PointLight::lValue(){
    return mRadiance * mColour;
}

// ===------------IMPLEMENTATIONS-------------===

Program::Program(int width, int height, std::string title) :
    settings{}, callbacks{}, mWindow{ nullptr }
{
    settings.size.width = width;
    settings.size.height = height;
    settings.title = title;

    if (!glx::initializeGLFW(errorCallback))
    {
        throw OpenGLError("Failed to initialize GLFW with error callback");
    }

    mWindow = glx::createGLFWWindow(settings);
    if (mWindow == nullptr)
    {
        throw OpenGLError("Failed to create GLFW Window");
    }

    createGLContext();
}

void Program::run(Pinhole& ph, Cube& cube, PointLight& pLight)
{
    glEnable(GL_DEPTH_TEST);

    while (!glfwWindowShouldClose(mWindow))
    {
        int width;
        int height;

        glfwGetFramebufferSize(mWindow, &width, &height);
        // setup the view to be the window's size
        glViewport(0, 0, width, height);
        // tell OpenGL the what color to clear the screen to
        glClearColor(0, 0, 0, 1);
        // actually clear the screen
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        ph.render(width, height, cube, pLight);

        glfwSwapBuffers(mWindow);
        glfwPollEvents();
    }
}

void Program::freeGPUData()
{
    glx::destroyGLFWWindow(mWindow);
    glx::terminateGLFW();
}

void Program::createGLContext()
{
    using namespace magic_enum::bitwise_operators;

    glx::bindWindowCallbacks(mWindow, callbacks);
    glfwMakeContextCurrent(mWindow);
    glfwSwapInterval(1);

    if (!glx::createGLContext(mWindow, settings.version))
    {
        throw OpenGLError("Failed to create OpenGL context");
    }

    glx::initializeGLCallback(glx::ErrorSource::All,
        glx::ErrorType::All,
        glx::ErrorSeverity::High |
        glx::ErrorSeverity::Medium);
}

// ===-----------------DRIVER-----------------===

int main()
{
    try
    {
        // clang-format off
        std::array<float, 3 * 8> vertices {
            // Vertices          
            -0.4f, -0.4f, 0.4f, //bottom left front
            0.4f, -0.4f, 0.4f, //bottom right front
            0.4f,  0.4f, 0.4f, //top right front
            -0.4f,  0.4f, 0.4f, //top left front
            0.4f, -0.4f, -0.4f, //bottom right back
            0.4f,  0.4f, -0.4f, //top right back
            -0.4f, -0.4f, -0.4f, //bottom left back
            -0.4f,  0.4f, -0.4f, //top left back
        };

        std::array<unsigned int, 3 * 2 * 6> indices{
            0, 1, 2,
            2, 3, 0,
            1, 4, 5,
            5, 2, 1,
            4, 6, 7,
            7, 5, 4,
            6, 0, 3,
            3, 7, 6,
            3, 2, 5,
            5, 7, 3,
            6, 4, 1,
            1, 0, 6
        };

        glm::vec3 cColour{ 1.0f, 0.4f, 0.2f };

        std::array<float, 3> plight{
            6.0f, 9.0f, 3.0f
        };

        // clang-format on

        Program prog{ 1280, 720, "CSC305 Assignment 3" };
        Pinhole pin{};
        Cube cube{vertices, indices, cColour};
        PointLight pLight{ plight };

        pin.setCamera(glm::vec3(2.0f, 1.0f, 3.0f), glm::vec3(0.0f, 0.3f, 0.0f), glm::vec3(3.0f, 1.0f, 0.0f));
        pin.loadShaders();
        pin.loadCubeDataToGPU(cube);

        prog.run(pin, cube, pLight);

        prog.freeGPUData();
        pin.freeGPUData();
    }
    catch (OpenGLError & err)
    {
        fmt::print("OpenGL Error:\n\t{}\n", err.what());
    }

    return 0;
}

