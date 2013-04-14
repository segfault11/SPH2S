#pragma once
#include "Solver.h"
#include "OpenGL.h"

class Renderer
{
public:
    Renderer (const ParticleData& particles, float xs, float ys, float xe, 
        float ye, float r, float g, float b, float a);
    ~Renderer ();
    void Render () const;

private:
    const ParticleData& mParticles;
	GLuint* mActiveIDs;
    GLuint mProgram;
    GLuint mPositionsVBO;
	GLuint mActiveIDsVBO;
    GLuint mPositionsVAO;
};
