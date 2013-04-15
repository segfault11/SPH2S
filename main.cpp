//------------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <SDL2/SDL.h>
#include "OpenGL.h"
#include <iostream>
#include "Solver.h"
#include "Renderer.h"
//------------------------------------------------------------------------------
#define WIDTH  800
#define HEIGHT 800
#define PI 3.14159265358979323846f
//------------------------------------------------------------------------------
static ParticleData* gsParticleData[2];
static ParticleData* gsBoundaryData;
static Renderer* gsRenderer[2];
static Renderer* gsBoundaryRenderer;
static Solver* gsSolver;
static DisplayRectangle gsDispRect(0.0f, 0.0f, 1.0f, 1.0f);
static bool gsPause = false;
//------------------------------------------------------------------------------
static void init ()
{
	// init particle data  
	gsParticleData[LOW] = ParticleData::CreateCube
	(
		Grid(Vector2f(0.21f, 0.21f), Vector2ui(26, 51), 0.005f)
	);
	
	gsParticleData[HIGH] = ParticleData::CreateCube
	(
		Grid(Vector2f(0.64f, 0.21f), Vector2ui(26, 51), 0.0025f)  
	);

	// init boundary data
	gsBoundaryData = ParticleData::CreateCanvas
	(
		Grid(Vector2f(0.2f, 0.2f), Vector2f(0.8f, 0.8f), 0.005f),
		3
	);	

    // init the renderer
    gsRenderer[LOW] = new Renderer
	(
		*(gsParticleData[LOW]),
		gsDispRect,
		0.005f
	);
    gsRenderer[HIGH] = new Renderer
	(
		*(gsParticleData[HIGH]), 
		gsDispRect,
		0.0025f
	);

	// init boundary renderer
	gsBoundaryRenderer = new Renderer
	(
		*gsBoundaryData, 
		gsDispRect, 
		0.010f
	);

    // init the solver configuration
    SolverConfiguration configuration;
    configuration.EffectiveRadius[LOW] = std::sqrt
    (
        0.125f*0.25f*30.0f/(PI*gsParticleData[LOW]->NumParticles)
    );
    configuration.EffectiveRadius[HIGH] = 
		configuration.EffectiveRadius[LOW]/2.0f;
    configuration.Domain[LOW] = Domain
    (
        Vector2f(0.0f, 0.0f), 
        Vector2f(1.0f, 1.0f),
        configuration.EffectiveRadius[LOW] 
    );
    configuration.Domain[HIGH] = Domain
    (
        Vector2f(0.0f, 0.0f), 
        Vector2f(1.0f, 1.0f),
        configuration.EffectiveRadius[HIGH] 
    );
    configuration.RestDensity = 1000.0f;
    configuration.FluidParticleMass[LOW] = configuration.RestDensity*0.25f*
		0.125f/gsParticleData[LOW]->NumParticles;
    configuration.FluidParticleMass[HIGH] = 
		configuration.FluidParticleMass[LOW]/4.0f;
    configuration.BoundaryParticleMass = configuration.FluidParticleMass[LOW];
    configuration.TensionCoefficient = 0.08f; 
    configuration.SpeedSound = 88.1472f;
    configuration.TaitCoefficient = 1119.0714f;
    configuration.Alpha = 0.04f;

    // init the solver
    gsSolver = new Solver(gsParticleData, gsBoundaryData, configuration);
	
	//
	// init gl
	//	
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
}
//------------------------------------------------------------------------------
static void release ()
{
	delete gsParticleData[LOW];
	delete gsParticleData[HIGH];
	delete gsRenderer[LOW];
	delete gsRenderer[HIGH];
	delete gsSolver;
}
//------------------------------------------------------------------------------
static void draw ()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	gsRenderer[LOW]->Render();
	gsRenderer[HIGH]->Render();
	gsBoundaryRenderer->Render();

	if (!gsPause)
	{
		gsSolver->Advance(0.0005f);
	}
}
//------------------------------------------------------------------------------
int main(int argc, const char *argv[])
{
    SDL_Window* mainWindow;
    SDL_GLContext mainContext;

    if (SDL_Init(SDL_INIT_VIDEO) < 0) 
    {
        
    }

    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 2);
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);

	mainWindow = SDL_CreateWindow
	(
		"SPH2SCALE", 
		SDL_WINDOWPOS_CENTERED, 
        SDL_WINDOWPOS_CENTERED, 
		WIDTH, 
		HEIGHT, 
		SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN
	);
   
	mainContext = SDL_GL_CreateContext(mainWindow);
	SDL_GL_SetSwapInterval(1);
	
    glewExperimental = GL_TRUE; 

    if (GLEW_OK != glewInit())
    {
        printf("glewInit failed\n");
    }
    else
    {
        printf("glewInit succeded\n");
    }
	
	init();

	bool running = true;
	SDL_Event event;

	while (running)
	{
		while (SDL_PollEvent(&event))
		{
			switch (event.type)
			{
				case SDL_KEYDOWN:
					switch (event.key.keysym.sym)
					{
						case SDLK_LEFT:
							gsDispRect.Translate(-0.01f, 0.0f);
							gsRenderer[LOW]->SetDisplayRectangle(gsDispRect);
							gsRenderer[HIGH]->SetDisplayRectangle(gsDispRect);
							gsBoundaryRenderer->SetDisplayRectangle(gsDispRect);
							break;
						case SDLK_RIGHT:
							gsDispRect.Translate(0.01f, 0.0f);
							gsRenderer[LOW]->SetDisplayRectangle(gsDispRect);
							gsRenderer[HIGH]->SetDisplayRectangle(gsDispRect);
							gsBoundaryRenderer->SetDisplayRectangle(gsDispRect);
							break;
						case SDLK_UP:
							gsDispRect.Translate(0.0f, 0.01f);
							gsRenderer[LOW]->SetDisplayRectangle(gsDispRect);
							gsRenderer[HIGH]->SetDisplayRectangle(gsDispRect);
							gsBoundaryRenderer->SetDisplayRectangle(gsDispRect);
							break;
						case SDLK_DOWN:
							gsDispRect.Translate(0.0f, -0.01f);
							gsRenderer[LOW]->SetDisplayRectangle(gsDispRect);
							gsRenderer[HIGH]->SetDisplayRectangle(gsDispRect);
							gsBoundaryRenderer->SetDisplayRectangle(gsDispRect);
							break;
						case SDLK_ESCAPE:
							running = 0;
							break;
						case 'q':
							gsDispRect.Scale(-0.01f);
							gsRenderer[LOW]->SetDisplayRectangle(gsDispRect);
							gsRenderer[HIGH]->SetDisplayRectangle(gsDispRect);
							gsBoundaryRenderer->SetDisplayRectangle(gsDispRect);
							break;	
						case 'a':
							gsDispRect.Scale(0.01f);
							gsRenderer[LOW]->SetDisplayRectangle(gsDispRect);
							gsRenderer[HIGH]->SetDisplayRectangle(gsDispRect);
							gsBoundaryRenderer->SetDisplayRectangle(gsDispRect);
							break;	
						case 'p':
							gsPause = gsPause == true ? false : true;
							break;
					}
					break;
		//		case SDL_KEYUP:
				case SDL_QUIT:
					running = 0;
			}
		}

		draw();
		SDL_GL_SwapWindow(mainWindow);
	}	

	release();
	SDL_GL_DeleteContext(mainContext);
    SDL_DestroyWindow(mainWindow);
    SDL_Quit();

    return 0;
}
//------------------------------------------------------------------------------
