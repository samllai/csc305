#version 450 core

uniform vec3 ccolour;
uniform vec3 lcolour;
uniform vec3 lvalue;
uniform vec3 lpos;

in vec3 fragPos;
in vec3 fragNorm;

out vec4 fragColour;

void main()
{
    vec3 lightDir = normalize(lpos - fragPos);
    float diff = max(dot(fragNorm, lightDir), 0.0);
    vec3 diffuse = diff * lcolour;
    vec3 result = (lvalue + diffuse) * ccolour;
    fragColour = vec4(result, 1.0);
}
