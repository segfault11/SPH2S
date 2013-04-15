#include "Renderer.h"
#include <iostream>

//------------------------------------------------------------------------------
Renderer::Renderer (const ParticleData& particles, const DisplayRectangle& rect)
: mParticles(particles)
{
	//
	mActiveIDs = new GLuint[particles.NumParticles];

    // create program
    mProgram = glCreateProgram();
    GL::AttachShader(mProgram, "RendererVertexShader.glsl", 
        GL_VERTEX_SHADER);
    GL::AttachShader(mProgram, "RendererFragmentShader.glsl", 
        GL_FRAGMENT_SHADER);
    GL::BindAttribLocation(mProgram, "position", 0);
    GL::BindAttribLocation(mProgram, "colors", 1);
    GL::BindFragDataLocation(mProgram, "fragOutput", 0);
    GL::LinkProgram(mProgram);
    GL::DumpLog(mProgram);

    // init program
	SetDisplayRectangle(rect);

    // create vbo & vao
    GL::CreateBufferObject
    (
        mPositionsVBO,
        GL_ARRAY_BUFFER,
        mParticles.NumParticles*2*sizeof(float),
        NULL,
        GL_DYNAMIC_DRAW
    );

    GL::CreateBufferObject
    (
        mColorsVBO,
        GL_ARRAY_BUFFER,
        mParticles.NumParticles*sizeof(float),
        NULL,
        GL_DYNAMIC_DRAW
    );

	GL::CreateBufferObject
	(
		mActiveIDsVBO,
		GL_ELEMENT_ARRAY_BUFFER,
		mParticles.NumParticles*sizeof(unsigned int),
		NULL,
		GL_DYNAMIC_DRAW
	);	

    glGenVertexArrays(1, &mPositionsVAO);
    glBindVertexArray(mPositionsVAO);
    glBindBuffer(GL_ARRAY_BUFFER, mPositionsVBO);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, 0);
	glBindBuffer(GL_ARRAY_BUFFER, mColorsVBO);
    glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 0, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mActiveIDsVBO);
    glEnableVertexAttribArray(0);
	glEnableVertexAttribArray(1);
}
//------------------------------------------------------------------------------
Renderer::~Renderer ()
{
	delete[] mActiveIDs;
    glDeleteProgram(mProgram);
    glDeleteBuffers(1, &mPositionsVBO);
	glDeleteBuffers(1, &mActiveIDsVBO);
	glDeleteBuffers(1, &mColorsVBO);
	glDeleteVertexArrays(1, &mPositionsVAO);
}
//------------------------------------------------------------------------------
void Renderer::Render () const
{
    glUseProgram(mProgram);
    glBindVertexArray(mPositionsVAO);
    glBindBuffer(GL_ARRAY_BUFFER, mPositionsVBO);
    glBufferData
    (
        GL_ARRAY_BUFFER, 
        mParticles.NumParticles*2*sizeof(float), 
        mParticles.Positions, 
		GL_STATIC_DRAW
    );
    glBindBuffer(GL_ARRAY_BUFFER, mColorsVBO);
    glBufferData
    (
        GL_ARRAY_BUFFER, 
        mParticles.NumParticles*sizeof(float), 
        mParticles.Colors, 
		GL_STATIC_DRAW
    );

	//
	IDList::const_iterator i = mParticles.ActiveIDs.begin();
	IDList::const_iterator e = mParticles.ActiveIDs.end();
	unsigned int c = 0;

	for (; i != e; i++)
	{
		mActiveIDs[c] = *i;
		c++;
	}
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mActiveIDsVBO);
	glBufferData
	(
		GL_ELEMENT_ARRAY_BUFFER,
		mParticles.ActiveIDs.size()*sizeof(unsigned int),
		mActiveIDs,
		GL_STATIC_DRAW
	);

   	glDrawElements(GL_POINTS, mParticles.ActiveIDs.size(),
		GL_UNSIGNED_INT, 0);
	
}
//------------------------------------------------------------------------------
void Renderer::SetDisplayRectangle (const DisplayRectangle& rect)
{
    float width = rect.XMax - rect.XMin;
    float height = rect.YMax - rect.YMin;

    glUseProgram(mProgram);
    GLint loc;
    loc = glGetUniformLocation(mProgram, "xs");
    glUniform1f(loc, rect.XMin);
    loc = glGetUniformLocation(mProgram, "ys");
    glUniform1f(loc, rect.YMin);
    loc = glGetUniformLocation(mProgram, "width");
    glUniform1f(loc, width);
    loc = glGetUniformLocation(mProgram, "height");
    glUniform1f(loc, height);
}
//------------------------------------------------------------------------------
