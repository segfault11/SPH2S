//------------------------------------------------------------------------------
#include "Solver.h"
#include "Timer.h"
//------------------------------------------------------------------------------

//==============================================================================
//  HashTable's Definition
//==============================================================================

//------------------------------------------------------------------------------
inline Vector2ui computeCoordinate
(
    const Vector2f& position, 
    const Domain& domain
)
{
    int i = static_cast<int>
    (
        std::ceil((position.X - domain.Origin.X)/domain.Spacing)
    );

    int j = static_cast<int>
    (
        std::ceil((position.Y - domain.Origin.Y)/domain.Spacing)
    );

    i = std::max<int>(0, std::min<int>(i, domain.Dimensions.X - 1));
    j = std::max<int>(0, std::min<int>(j, domain.Dimensions.Y - 1));

    return Vector2ui(i, j);
}
//------------------------------------------------------------------------------
inline unsigned int computeHash 
(
    const Vector2ui& coordinate,
    unsigned int domainWidth
)
{
    return coordinate.Y*domainWidth + coordinate.X; 
}
//------------------------------------------------------------------------------
inline unsigned int computeHash 
(
    unsigned int i,
    unsigned int j,
    unsigned int domainWidth
)
{
    return j*domainWidth + i; 
}
//------------------------------------------------------------------------------
inline unsigned int computeHash
(
    const Vector2f& position, 
    const Domain& domain
)
{
    return computeHash
	(
		computeCoordinate(position, domain), 
		domain.Dimensions.X
	);
}
//------------------------------------------------------------------------------
Solver::HashTable::HashTable 
(
    const float* positions, 
    const Domain& domain
)
: mPositions(positions), mDomain(domain)
{
    mBuckets = new std::list<unsigned int>[domain.Dimensions.X*
        domain.Dimensions.Y];
}
//------------------------------------------------------------------------------
Solver::HashTable::~HashTable ()
{
    delete[] mBuckets;
}
//------------------------------------------------------------------------------
void Solver::HashTable::Fill 
(
    const IDList& particleIDs
)
{
    // clear all buckets
    for (unsigned int i = 0; i < mDomain.Dimensions.X*mDomain.Dimensions.Y; i++)
    {
        mBuckets[i].clear();
    }

    // fill buckets
	IDList::const_iterator i = particleIDs.begin();
	IDList::const_iterator e = particleIDs.end();

    for (; i != e; i++)
    {
        unsigned int hash = computeHash
        (
            Vector2f(&mPositions[2*(*i)]), 
            mDomain
        );
        mBuckets[hash].push_back(*i);
    }

}
//------------------------------------------------------------------------------
void Solver::HashTable::Fill (unsigned int numParticles)
{
    // clear all buckets
    for (unsigned int i = 0; i < mDomain.Dimensions.X*mDomain.Dimensions.Y; i++)
    {
        mBuckets[i].clear();
    }

    // fill buckets
    for (unsigned int i = 0; i < numParticles; i++)
    {
        unsigned int hash = computeHash
        (
            Vector2f(&mPositions[2*i]), 
            mDomain
        );
        mBuckets[hash].push_back(i);
    }

}
//------------------------------------------------------------------------------
void Solver::HashTable::Query (const Vector2f& position)
{
    // reset previous result
    mResult.clear();
    
    Vector2f begf(position);
    begf.X -= mDomain.Spacing;
    begf.Y -= mDomain.Spacing;
    Vector2f endf(position);
    endf.X += mDomain.Spacing;
    endf.Y += mDomain.Spacing;

    Vector2ui beg = computeCoordinate(begf, mDomain);
    Vector2ui end = computeCoordinate(endf, mDomain);

    for (unsigned int i = beg.X; i <= end.X; i++)
    {
        for (unsigned int j = beg.Y; j <= end.Y; j++)
        {
            unsigned int hash = computeHash(i, j, mDomain.Dimensions.X);
            mResult.insert
            (
                mResult.end(), 
                mBuckets[hash].begin(), 
                mBuckets[hash].end()
            );
        }
    }   
}
//------------------------------------------------------------------------------
void Solver::HashTable::Query (const Vector2f& position, float range)
{
	// reset previous result
    mResult.clear();
    
    Vector2f begf(position);
    begf.X -= range;
    begf.Y -= range;
    Vector2f endf(position);
    endf.X += range;
    endf.Y += range;

    Vector2ui beg = computeCoordinate(begf, mDomain);
    Vector2ui end = computeCoordinate(endf, mDomain);

    for (unsigned int i = beg.X; i <= end.X; i++)
    {
        for (unsigned int j = beg.Y; j <= end.Y; j++)
        {
            unsigned int hash = computeHash(i, j, mDomain.Dimensions.X);
            mResult.insert
            (
                mResult.end(), 
                mBuckets[hash].begin(), 
                mBuckets[hash].end()
            );
        }
    }   
}
//------------------------------------------------------------------------------
const IDList& Solver::HashTable::GetResult () const
{
    return mResult;
}
//------------------------------------------------------------------------------

//==============================================================================
//  Solvers's Definition
//==============================================================================

//------------------------------------------------------------------------------
#define PI 3.14159265358979323846f
//------------------------------------------------------------------------------
inline float computeDotProduct (const float *x, const float *y)
{
    return x[0]*y[0] + x[1]*y[1];
}
//------------------------------------------------------------------------------
inline float computeDistance (const float *x, const float *y)
{
    return std::sqrt((x[0] - y[0])*(x[0] - y[0]) + (x[1] - y[1])*(x[1] - y[1]));
}
//------------------------------------------------------------------------------
inline float computeNorm (const float *x)
{
	return std::sqrt(x[0]*x[0] + x[1]*x[1]);
}
//------------------------------------------------------------------------------
inline float evaluateKernel (float dist, float h)
{
    float c = 20.0f/(14.0f*PI*h*h);
    float q = 2.0f*dist/h;

    if (q < 1.0f)
    {
        float a = 2.0f - q;
        float b = 1.0f - q;

        return c*(a*a*a - 4.0f*b*b*b);
    }
    if (q < 2.0f)
    {
        float a = 2.0f - q;

        return c*a*a*a;
    }

    return 0.0f;
}
//------------------------------------------------------------------------------
inline Vector2f evaluateKernelGradient 
(
	const Vector2f& xij, 
	float dist, 
	float h
)
{
    Vector2f grad(0.0f, 0.0f);
    float c = -120.0f/(14.0f*PI*h*h*h*dist);
    float q = 2.0f*dist/h;

    if (dist == 0.0f)
    {
        return grad;
    }

    if (q < 1.0f)
    {
        float a = 2.0f - q;
        float b = 1.0f - q;
        float d = c*(a*a - 4.0f*b*b);
        grad.X = d*xij.X;
        grad.Y = d*xij.Y;

        return grad;
    }

    if (q < 2.0f)
    {
        float a = 2.0f - q;
        float d = c*a*a;
        grad.X = d*xij.X;
        grad.Y = d*xij.Y;

        return grad;
    }
    return grad;
}
//------------------------------------------------------------------------------
inline float evaluateBoundaryWeight (float dist, float h, float speedSound)
{
	float q = 2.0f*dist/h;
	float c = 0.02f*speedSound*speedSound/(dist*dist);

	if (q < 2.0f/3.0f)
	{
		return c*2.0f/3.0f;
	}
	else if (q < 1.0f)
	{
		float a = 2.0f*q - 3.0f/2.0f*q*q;
		return c*a;
	}
	else if (q < 2.0f)
	{
		float a = (2.0f - q);
		return c*0.5f*a*a;
	}
	else
	{
		return 0.0f;
	}
}
//------------------------------------------------------------------------------
Solver::Solver 
(
    ParticleData* fluidParticles[2], 
    ParticleData* boundaryParticles, 
    const SolverConfiguration configuration
)
:
    mConfiguration(configuration)
{
	// assign particle data to member
   	mFluidParticles[LOW] = fluidParticles[LOW];
   	mFluidParticles[HIGH] = fluidParticles[HIGH];
	mBoundaryParticles = boundaryParticles;

	// allocate Hash tables
	mFluidHashTable[LOW] = new HashTable
	(
		fluidParticles[LOW]->Positions, 
		configuration.Domain[LOW]
	);
	mFluidHashTable[HIGH] = new HashTable
	(
		fluidParticles[HIGH]->Positions, 
		configuration.Domain[HIGH]
	);
	mBoundaryHashTable = new HashTable
	(
		boundaryParticles->Positions, 
		configuration.Domain[LOW]
	);
	
	// allocate and set additional particle data for the low and high res
	// of the fluid
	for (unsigned int i = 0; i < 2; i++)
	{ 
		// allocate additional particle data
    	mVelocities[i] = new float[2*fluidParticles[i]->NumParticles];
    	mAccelerations[i] = new float[2*fluidParticles[i]->NumParticles];
    	mDensities[i] = new float[fluidParticles[i]->NumParticles];
    	mPressures[i] = new float[fluidParticles[i]->NumParticles];

		// memset everthing to 0
	    memset
		(
			mVelocities[i], 
			0, 
			sizeof(float)*2*fluidParticles[i]->NumParticles
		);
    	memset
		(
			mAccelerations[i] , 
			0, 
			sizeof(float)*2*fluidParticles[i]->NumParticles
		);
    	memset
		(
			mDensities[i], 
			0, 
			sizeof(float)*fluidParticles[i]->NumParticles
		);
    	memset
		(
			mPressures[i], 
			0, 
			sizeof(float)*fluidParticles[i]->NumParticles
		);
	}

	// insert boundary particles into the hash map
	mBoundaryHashTable->Fill(mBoundaryParticles->NumParticles);

}
//------------------------------------------------------------------------------
Solver::~Solver ()
{
	// release everything
	for (unsigned int i = 0; i < 2; i++)
	{
    	delete[] mVelocities[i];
    	delete[] mAccelerations[i];
    	delete[] mDensities[i];
    	delete[] mPressures[i];
		delete mFluidHashTable[i];
	}

}
//------------------------------------------------------------------------------
void Solver::Advance (float timeStep) 										
{
	mFluidHashTable[LOW]->Fill(mFluidParticles[LOW]->ActiveIDs);                
	mFluidHashTable[HIGH]->Fill(mFluidParticles[HIGH]->ActiveIDs);
	computeDensity(LOW); 
	computeDensity(HIGH); 
	computeAcceleration(LOW);
	computeAcceleration(HIGH);
	integrate(LOW, timeStep);
    integrate(HIGH, timeStep);
	inject();
}
//------------------------------------------------------------------------------
void Solver::computeDensity (unsigned char res)
{
	IDList::const_iterator i = mFluidParticles[res]->ActiveIDs.begin();
	IDList::const_iterator e = mFluidParticles[res]->ActiveIDs.end();

    for (; i != e; i++)
    {
		//======================================================================
		// 	compute densities in the same domain
		//======================================================================

        float *pos = &mFluidParticles[res]->Positions[2*(*i)];
        
        // get neighbor ids
        mFluidHashTable[res]->Query(Vector2f(pos));
        const IDList& neighbors = mFluidHashTable[res]->GetResult();
        
        //
        float density = 0.0f;

		IDList::const_iterator j = neighbors.begin();
		IDList::const_iterator je = neighbors.end();
		
		// iterate through neighbors and add up their contribution to 
		// the density of particle i
        for (; j != je; j++)
        {
            float *posj = &mFluidParticles[res]->Positions[2*(*j)];
            float posij[2];
            posij[0] = pos[0] - posj[0];
            posij[1] = pos[1] - posj[1];
            float dist = std::sqrt(posij[0]*posij[0] + posij[1]*posij[1]);

            if (dist < mConfiguration.EffectiveRadius[res])
            {
                density += evaluateKernel
				(
					dist, 
					mConfiguration.EffectiveRadius[res]
				);
            }

        }		

        density *= mConfiguration.FluidParticleMass[res];

		//======================================================================
		// compute densities in the complementary domain
		//======================================================================

		unsigned char resc = (res + 1) % 2;
		float densityc = 0.0f;


		// get neighbor ids
		// (in the contemp. domain, we always have a search radius of the 
		// low effective radius)
        mFluidHashTable[resc]->Query
		( 
			Vector2f(pos),
			mConfiguration.EffectiveRadius[LOW] 
		);
        const IDList& neighborsc = mFluidHashTable[resc]->GetResult();

       	IDList::const_iterator jc = neighborsc.begin();
		IDList::const_iterator jce = neighborsc.end();

 		for (; jc != jce; jc++)
        {
            float *posjc = &mFluidParticles[resc]->Positions[2*(*jc)];
            float posijc[2];
            posijc[0] = pos[0] - posjc[0];
            posijc[1] = pos[1] - posjc[1];
            float distc = std::sqrt(posijc[0]*posijc[0] + posijc[1]*posijc[1]);

            if (distc < mConfiguration.EffectiveRadius[LOW])
            {
					
                float w = evaluateKernel
				(
					distc, 
					mConfiguration.EffectiveRadius[res]
				);

				float wc = evaluateKernel
				(
					distc,
					mConfiguration.EffectiveRadius[resc]
				);

				densityc += 0.5f*(w + wc);

            }

        }		
		
		densityc *= mConfiguration.FluidParticleMass[resc];

		//======================================================================
		// 	Update density and pressure of particle i
		//======================================================================

        mDensities[res][*i] = density + densityc;
        
        float a = density/mConfiguration.RestDensity;
        float a3 = a*a*a;
        mPressures[res][*i] = mConfiguration.TaitCoefficient*(a3*a3*a - 1.0f);
    } 

}
//------------------------------------------------------------------------------
void Solver::computeAcceleration (unsigned char res)
{
	IDList::const_iterator i = mFluidParticles[res]->ActiveIDs.begin();
	IDList::const_iterator e = mFluidParticles[res]->ActiveIDs.end();

    for (; i != e; i++)
    {
        float *pos = &mFluidParticles[res]->Positions[2*(*i)];
        float *vel = &mVelocities[res][2*(*i)]; 
        float *acc = &mAccelerations[res][2*(*i)];
        float den = mDensities[res][*i];
        float pre = mPressures[res][*i];
	

	//	std::cout << mDensities[LOW][*i] << std::endl;
    

	    // reset acc
        acc[0] = 0.0f;
        acc[1] = 0.0f;

        // reserve mem for tension force
        float accT[2];
        accT[0] = 0.0f;
        accT[1] = 0.0f;

		// reserve mem for boundary force
		float accB[2];
		accB[0] = 0.0f;
		accB[1] = 0.0f;
	
		//======================================================================
		// 	Compute force contribution from the same domain
		//======================================================================

        // get neighbor ids
        mFluidHashTable[res]->Query(Vector2f(pos));
        const IDList& neighbors = mFluidHashTable[res]->GetResult();

		IDList::const_iterator j = neighbors.begin();
		IDList::const_iterator je = neighbors.end();
        
        for (; j != je; j++)
        {
            float *posj = &mFluidParticles[res]->Positions[2*(*j)];
            float *velj = &mVelocities[res][2*(*j)];
            float denj = mDensities[res][*j];
            float prej = mPressures[res][*j];
            float dist = computeDistance(pos, posj);
            float xij[2];
            xij[0] = pos[0] - posj[0];
            xij[1] = pos[1] - posj[1];
            float vij[2];
            vij[0] = vel[0] - velj[0];
            vij[1] = vel[1] - velj[1];
        
            if (dist < mConfiguration.EffectiveRadius[res])
            {
                // compute SPH pressure coefficient
                float coeff = (pre/(den*den) + prej/(denj*denj));
                float vx = computeDotProduct(xij, vij);
                float h = mConfiguration.EffectiveRadius[res];
                
                // check if artificial viscosity shoud be applied
                if (vx < 0.0f)
                {
                    // add artificial velocity to the SPH Force coefficient
                    float x2 = dist*dist;
                    float nu = -2.0f*mConfiguration.Alpha*h*
                        mConfiguration.SpeedSound/(den + denj);
                    coeff += nu*vx/(x2 + 0.01f*h*h);
                }

                //evaluate kernel gradient and add contribution of particle j
                //to acc
				Vector2f grad = evaluateKernelGradient(Vector2f(xij), dist, h);
                acc[0] += coeff*grad.X;
                acc[1] += coeff*grad.Y;

                // compute tension force and add to the acc
                float kernelT = evaluateKernel(dist, h);
                accT[0] -= mConfiguration.TensionCoefficient*xij[0]*kernelT;
                accT[1] -= mConfiguration.TensionCoefficient*xij[1]*kernelT;
            }

        }

		//======================================================================
		// 	Compute force contribution from the contemp. domain
		//======================================================================

		float accC[2];
		accC[0] = 0.0f;
		accC[1] = 0.0f;

		float accCT[2];
		accCT[0] = 0.0f;
		accCT[1] = 0.0f;

		unsigned char resc = (res + 1) % 2;
        mFluidHashTable[resc]->Query
		(
			Vector2f(pos),
    		mConfiguration.EffectiveRadius[LOW] 
		);
       	const IDList& neighborsc = mFluidHashTable[resc]->GetResult();

		IDList::const_iterator jc = neighborsc.begin();
		IDList::const_iterator jce = neighborsc.end();

		for (; jc != jce; jc++)
       	{
          	float *posj = &mFluidParticles[resc]->Positions[2*(*jc)];
           	float *velj = &mVelocities[resc][2*(*jc)];
           	float denj = mDensities[resc][*jc];
           	float prej = mPressures[resc][*jc];
           	float dist = computeDistance(pos, posj);
           	float xij[2];
       	    xij[0] = pos[0] - posj[0];
   	        xij[1] = pos[1] - posj[1];
            float vij[2];
           	vij[0] = vel[0] - velj[0];
           	vij[1] = vel[1] - velj[1];
    
        	if (dist < mConfiguration.EffectiveRadius[LOW])
           	{
			    // compute SPH pressure coefficient
               	float coeff = (pre/(den*den) + prej/(denj*denj));
               	float vx = computeDotProduct(xij, vij);
               	float h = mConfiguration.EffectiveRadius[res];
				float hc = mConfiguration.EffectiveRadius[resc];	
				
		        if (vx < 0.0f)
               	{
                   	// add artificial velocity to the SPH Force coefficient
                   	float x2 = dist*dist;
                   	float nu = -2.0f*mConfiguration.Alpha*h*
                       	mConfiguration.SpeedSound/(den + denj);
                   	coeff += nu*vx/(x2 + 0.01f*h*h);
              	}
				
				Vector2f grad = evaluateKernelGradient
				(
					Vector2f(xij), 
					dist,
					h
				);
				Vector2f gradc = evaluateKernelGradient
				(
					Vector2f(xij), 
					dist, 
					hc
				);

                accC[0] += coeff*(grad.X + gradc.X)*0.5f;
                accC[1] += coeff*(grad.Y + gradc.Y)*0.5f;

                float wt = evaluateKernel(dist, h);
                float wct = evaluateKernel(dist, hc);

                accCT[0] -= mConfiguration.TensionCoefficient*xij[0]*
					(wt + wct)/2.0f;
                accCT[1] -= mConfiguration.TensionCoefficient*xij[1]*
					(wt + wct)/2.0f;
       		}

		}

		//======================================================================
		// 	compute penalty force of the boundary
		//======================================================================

		mBoundaryHashTable->Query(Vector2f(pos));
		const IDList& bneighbors = mBoundaryHashTable->GetResult();

		IDList::const_iterator k = bneighbors.begin();
		IDList::const_iterator ke = bneighbors.end();

   	    for (; k != ke; k++)
		{
            float *posk = &mBoundaryParticles->Positions[2*(*k)];
       	    float xik[2];
   	        xik[0] = pos[0] - posk[0];
            xik[1] = pos[1] - posk[1];
			float dist = computeNorm(xik);

			if (dist < mConfiguration.EffectiveRadius[LOW]) 
			{
			    float w = evaluateBoundaryWeight
				(
					dist, 
					mConfiguration.EffectiveRadius[LOW],
					mConfiguration.SpeedSound
				);
				float c = mConfiguration.BoundaryParticleMass/
					(mConfiguration.FluidParticleMass[LOW] + 
					mConfiguration.BoundaryParticleMass);

				float d = c*w;
				accB[0] += d*xik[0];		
				accB[1] += d*xik[1];		
			}
				
		}
		
		//======================================================================
		// 	Update acceleration of particle i
		//======================================================================

        acc[0] *= -mConfiguration.FluidParticleMass[res];
        acc[1] *= -mConfiguration.FluidParticleMass[res];
		accC[0] *= -mConfiguration.FluidParticleMass[resc];
		accC[1] *= -mConfiguration.FluidParticleMass[resc];
		accCT[0] *= (mConfiguration.FluidParticleMass[resc]/
			mConfiguration.FluidParticleMass[res]);
		accCT[1] *= (mConfiguration.FluidParticleMass[resc]/
			mConfiguration.FluidParticleMass[res]);
        acc[0] += accT[0];
        acc[1] += accT[1];
        acc[0] += accC[0];
        acc[1] += accC[1];
        acc[0] += accCT[0];
        acc[1] += accCT[1];
        acc[0] += accB[0];
        acc[1] += accB[1];
        acc[1] -= 9.81f;
    }
}
//------------------------------------------------------------------------------
void Solver::integrate (unsigned char res, float timeStep)
{
	IDList::iterator i = mFluidParticles[res]->ActiveIDs.begin();
	IDList::iterator e = mFluidParticles[res]->ActiveIDs.end();

    for (;i != e; i++)
    {
        float* pos = &mFluidParticles[res]->Positions[2*(*i)];
        float* vel = &mVelocities[res][2*(*i)];
        float* acc = &mAccelerations[res][2*(*i)];

        vel[0] += timeStep*acc[0];
        vel[1] += timeStep*acc[1];
        pos[0] += timeStep*vel[0];
        pos[1] += timeStep*vel[1];


		float aNorm = sqrt(acc[0]*acc[0] + acc[1]*acc[1]);

		mFluidParticles[res]->Colors[*i] =
			mConfiguration.FluidParticleMass[res]*timeStep*aNorm/
			(0.25f*mConfiguration.EffectiveRadius[res]);
	}
}
//------------------------------------------------------------------------------
void Solver::inject () 
{
	static float dir[] = 
	{
		 1.0f, 1.0f,
		-1.0f, 1.0f,
		 1.0f, -1.0f,
		-1.0f, -1.0f
	};
	float sq = sqrt(2.0f);

	//==========================================================================
	// 	check for each particle i, if it should be split
	//==========================================================================

	IDList::iterator i = mFluidParticles[LOW]->ActiveIDs.begin();
	IDList::iterator e = mFluidParticles[LOW]->ActiveIDs.end();

    for (;i != e; i++)
    {
        float *pos = &mFluidParticles[LOW]->Positions[2*(*i)];
        float *vel = &mVelocities[LOW][2*(*i)]; 

		if (pos[0] > 0.5f)
		{
			//==================================================================
			// insert high res particles 
			//==================================================================
			
			// compute distance from parent particle
			float r = sqrt(mConfiguration.FluidParticleMass[LOW]/
				(M_PI*mDensities[LOW][*i]));

			// compute id of first child particle
			unsigned int k = 4*(*i);

			// for all child particles
			for (unsigned int j = 0; j < 4; j++)
			{
				// init position and velocity
				mFluidParticles[HIGH]->Positions[2*(k + j) + 0] = 
					pos[0] + dir[2*j + 0]*r/sq;
				mFluidParticles[HIGH]->Positions[2*(k + j) + 1] = 
					pos[1] + dir[2*j + 1]*r/sq;

				mVelocities[HIGH][2*(k + j) + 0] = vel[0];
				mVelocities[HIGH][2*(k + j) + 1] = vel[1];
				
				// set active
				mFluidParticles[HIGH]->ActiveIDs.push_back(k + j);
			}
			
			//==================================================================
			// 	set particle i inactive / remove it from the active list
			//==================================================================
			mFluidParticles[LOW]->ActiveIDs.erase(i);	
		}

	}	

}
//------------------------------------------------------------------------------
