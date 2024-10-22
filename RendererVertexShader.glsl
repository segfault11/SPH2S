#version 150

uniform float xs;
uniform float ys;
uniform float width;
uniform float height;

in vec3 position;
in float colors;


out VertexData
{
    float quantity;
}
gVertexData;

void main ()
{
    float x = (position.x - xs)/width*2.0f - 1.0f;
    float y = (position.y - ys)/height*2.0f - 1.0f;
	
	gVertexData.quantity = colors;    


    gl_Position = vec4(x, y, 0.0f, 1.0f);
}
