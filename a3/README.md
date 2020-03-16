README: CSC 305 Assignment 3

No additional features implemented
Code renders an orange cube - was not tested on a lab machine

Pinhole:
- manages shaders and data sharing with the GPU
- renders the scene
-- requires a cube and a light
- the pinhole camera values (eye, look at, up) can be set with setCamera
-- these values are used in the view matrix

Cube:
- stores the vertex, index, and colour values for a cube
- data accessed with get/set

PointLight:
- stores the location, radiance, and colour for a point light
-- the position data is initially sent in as an array, then also stored as a vector

Program:
- contains commands to run with OpenGL

Main:
- initializes the data for the camera, cube, and light

cube.vert:
- the vertex shader for the cube
- the MVP matrices are received as uniform data
- outputs the fragment's position and the fragment's normal
-- normal computation was simplified using the knowledge that the cube was specified on the origin
--the mat3(transpose(inverse(model))) multiplication was to convert the "surface" normal into world space

cube.frag:
- the fragment shader for the cube
- cube colour, light colour, light value (radiance * light colour), and light position are received as uniform data
- takes in the fragment position and normal from cube.vert to compute diffuse shading