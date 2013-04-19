//------------------------------------------------------------------------------
#ifndef RENDERER_H
#define RENDERER_H
//------------------------------------------------------------------------------
#include "Solver.h"
#include "OpenGL.h"
//------------------------------------------------------------------------------
struct DisplayRectangle
{
	DisplayRectangle (float xs, float ys, float xe, float ye)
	: XMin(xs), YMin(ys), XMax(xe), YMax(ye)
	{

	}
	
	void Translate (float dx, float dy)
	{
		XMin += dx;
		XMax += dx;
		YMin += dy;
		YMax += dy;
	}

	void Scale (float s)
	{
		float w = XMax - XMin;
		float h = YMax - YMin;
		float sh = h/w*s;
	
		XMin -= s;
		XMax += s;
		YMin -= sh;
		YMax += sh;

		if (XMin > XMax || YMin > YMax)
		{
			XMin += s;
			XMax -= s;
			YMin += sh;
			YMax -= sh;
		}

	}
	
	float XMin, YMin;
	float XMax, YMax;
};
//------------------------------------------------------------------------------
class Renderer
{
public:
    Renderer (const ParticleData& particles, const DisplayRectangle& rect,
		float pointSize, bool setBlack = false);
    ~Renderer ();
    void Render () const;
	void SetDisplayRectangle (const DisplayRectangle& rect);

private:
    const ParticleData& mParticles;
	GLuint* mActiveIDs;
    GLuint mProgram;
    GLuint mPositionsVBO;
	GLuint mColorsVBO;
	GLuint mActiveIDsVBO;
    GLuint mPositionsVAO;
};
//------------------------------------------------------------------------------
#endif /* end of include guard: RENDERER_H */
//------------------------------------------------------------------------------
