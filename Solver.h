//------------------------------------------------------------------------------
#pragma once
//------------------------------------------------------------------------------
#include <list>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <stdio.h>
#include <memory.h>
//------------------------------------------------------------------------------
#define LOW 0
#define HIGH 1
//------------------------------------------------------------------------------
typedef std::list<unsigned int> IDList;
//------------------------------------------------------------------------------
template<typename T>
struct Tuple2
{
    T X;
    T Y;

    Tuple2 ()
    {
        X = T(0);
        Y = T(0);
    }

    Tuple2 (const T* dataPtr)
    {
        X = dataPtr[0];
        Y = dataPtr[1];
    }

    Tuple2 (T x, T y)
    : X(x), Y(y)
    {
    
    }
    
    Tuple2 (const Tuple2& original)
    {
        X = original.X;
        Y = original.Y;
    };
};
typedef Tuple2<float> Vector2f;
typedef Tuple2<unsigned int> Vector2ui;
//------------------------------------------------------------------------------
struct Domain
{
    Domain ()
    : 
        Origin(0.0f, 0.0f), 
        Dimensions(10, 10),
        Spacing(0.1f) 
    {
    
    }

    Domain (const Vector2f& origin, const Vector2f& end, float spacing)
    : Origin(origin), Spacing(spacing)
    {
        int i = static_cast<int>(ceil((end.X - origin.X)/spacing));
        int j = static_cast<int>(ceil((end.Y - origin.Y)/spacing));
        Dimensions.X = i;
        Dimensions.Y = j;
    }

    Domain (const Vector2f& origin, const Vector2ui& dimensions, float spacing)
    : Origin(origin), Dimensions(dimensions), Spacing(spacing)
    {

    }

    Vector2f Origin;
    Vector2ui Dimensions;
    float Spacing;
};
typedef Domain Grid;
//------------------------------------------------------------------------------
struct ParticleData
{
    ParticleData (unsigned int numParticles)
    : NumParticles(numParticles)
    {
        Positions = new float[2*numParticles];
        memset(Positions, 0, sizeof(float)*2*numParticles);

		for (unsigned int i = 0; i < numParticles; i++)
		{
			ActiveIDs.push_back(i);
		}
    }

    ~ParticleData ()
    {
        delete[] Positions;
    }

    static ParticleData* CreateCube (const Grid& grid)
    {
        ParticleData* data = new ParticleData
        (
            grid.Dimensions.X*grid.Dimensions.Y
        );

        for (unsigned int i = 0; i < grid.Dimensions.X; i++)
        {
            for (unsigned int j = 0; j < grid.Dimensions.Y; j++)
            {
                float x = grid.Origin.X + i*grid.Spacing;
                float y = grid.Origin.Y + j*grid.Spacing;

                data->Positions[2*(j*grid.Dimensions.X + i) + 0] = x;
                data->Positions[2*(j*grid.Dimensions.X + i) + 1] = y;
            }
        }

        data->NumParticles = grid.Dimensions.X*grid.Dimensions.Y;

        return data;
    }

	static ParticleData* CreateCanvas (const Grid& grid, int offset)
	{
		ParticleData* data = new ParticleData
        (
            (grid.Dimensions.X + 2*offset)*(grid.Dimensions.Y + 2*offset) -
				grid.Dimensions.X*grid.Dimensions.Y
        );

		unsigned int c = 0;


		for (int i = -offset; i < (int)grid.Dimensions.X + offset; i++)
        {
            for (int j = -offset; j < (int)grid.Dimensions.Y + offset; j++)
            {
				if 
				(
					(i < 0 || i >= (int)grid.Dimensions.X) ||
					(j < 0 || j >= (int)grid.Dimensions.Y)
				)
				{
                	float x = grid.Origin.X + i*grid.Spacing;
                	float y = grid.Origin.Y + j*grid.Spacing;

                	data->Positions[2*c + 0] = x;
                	data->Positions[2*c + 1] = y;
					c++;
				}
            }
        }

		return data;
	}

    float* Positions;
    unsigned int NumParticles;
	IDList ActiveIDs;

};
//------------------------------------------------------------------------------
struct SolverConfiguration
{
    SolverConfiguration ()
    {
    
    }

    float FluidParticleMass[2];
    float BoundaryParticleMass;
    float EffectiveRadius[2];
    float RestDensity;
    float TaitCoefficient;
    float SpeedSound;
    float Alpha;
    float TensionCoefficient;

    ::Domain Domain[2];
};
//------------------------------------------------------------------------------
class Solver
{
    class HashTable
    {
    public:
        HashTable (const float* positions, const Domain& domain);
        ~HashTable ();

        void Fill (const IDList& particleIDs);
        void Fill (unsigned int numParticles);
        void Query (const Vector2f& position);
        const IDList& GetResult () const;           
    
    private:
        const float* mPositions;
        std::list<unsigned int> mResult;
        std::list<unsigned int> *mBuckets;
        Domain mDomain;
    };

public:
    Solver (ParticleData* fluidParticles[2], ParticleData* boundaryParticles, 
        const SolverConfiguration configuration);
    ~Solver ();

    void Advance (float timeStep);

private:
    //==========================================================================
    // private methods used to advance the particle system
    //==========================================================================

    inline void computeDensity (unsigned char res);
    inline void computeAcceleration (unsigned char res);
    inline void integrate (unsigned char res, float timeStep);

    //==========================================================================
    // private member
    //==========================================================================
    
    // config of the solver
    SolverConfiguration mConfiguration;

    // fluid particle related data
    ParticleData* mFluidParticles[2];
    HashTable* mFluidHashTable[2];
    float* mVelocities[2];
    float* mAccelerations[2];
    float* mDensities[2];
    float* mPressures[2];
    
    ParticleData* mBoundaryParticles;
    HashTable* mBoundaryHashTable;
};
//------------------------------------------------------------------------------
