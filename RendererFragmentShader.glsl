#version 150

uniform float xs;
uniform float ys;
uniform float width;
uniform float height;

out vec4 fragOutput;

in GeometryData
{
    vec2 relCoord;
    float quantity;
}
gGeometryData;

// Returns the appropriate value from the Jet color function.
vec3 getJetColor(float value) 
{
     float fourValue = 4 * value;
     float red   = min(fourValue - 1.5, -fourValue + 4.5);
     float green = min(fourValue - 0.5, -fourValue + 3.5);
     float blue  = min(fourValue + 0.5, -fourValue + 2.5);
 
     return clamp( vec3(red, green, blue), 0.0, 1.0 );
}

void main ()
{
    float x = gGeometryData.relCoord.x;
    float y = gGeometryData.relCoord.y;

    if (x*x + y*y > 1.0f)
    {
        discard;
    }

    fragOutput = vec4(getJetColor(gGeometryData.quantity), 1.0f); //color;
}
