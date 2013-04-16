#ifndef _OpenGL_h
#define _OpenGL_h

#include <GL/GLEW.h>


namespace GL 
{    
//=============================================================================
//  Transformations functions
//=============================================================================
void LookAt(float ex, float ey, float ez, float cx, float cy, float cz,
    float ux, float uy, float uz, GLfloat mat[16]);
void Frustum(float l, float r, float b, float t, float n, float f, 
    GLfloat mat[16]);
void Perspective(float fovy, float aspect, float n, float f, 
    GLfloat mat[16]);

//=============================================================================
//  Shader functions
//=============================================================================

void AttachShader(GLuint program, const char* fileName, GLenum type);
void BindAttribLocation(GLuint program, const char* attrName, 
                              GLuint index);
void BindFragDataLocation(GLuint program, const char* szColorBufName,
                                GLuint index);
void LinkProgram(GLuint program);
void DumpLog(GLuint program);


//=============================================================================
// 	Buffer Objects
//=============================================================================

void CreateBufferObject (GLuint& buffer, GLenum target, 
    GLsizeiptr size, const GLvoid* data, GLenum usage);
}
#endif
