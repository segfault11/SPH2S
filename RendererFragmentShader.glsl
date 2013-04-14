#version 150

uniform float xs;
uniform float ys;
uniform float width;
uniform float height;
uniform vec4 color;

out vec4 fragOutput;

void main ()
{
    fragOutput = color;
}
