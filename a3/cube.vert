#version 450 core

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

layout(location = 0) in vec3 position;

out vec3 fragPos;
out vec3 fragNorm;

void main()
{
    gl_Position = projection*view*model*vec4(position, 1.0);
    fragPos = vec3(model* vec4(position, 1.0));
    fragNorm = normalize(mat3(transpose(inverse(model))) * position);
}
