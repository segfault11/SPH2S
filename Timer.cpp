//-----------------------------------------------------------------------------
//  Timer.cpp
//  C(++) Programming Practice
//
//  Created by Arno in Wolde Lübke on 06.01.13.
//-----------------------------------------------------------------------------

#include <iostream>
#include "Timer.h"

using namespace std;

//-----------------------------------------------------------------------------
Timer::Timer()
    : mState(kTimerStopped)
{
}
//-----------------------------------------------------------------------------
Timer::~Timer()
{
}
//-----------------------------------------------------------------------------
void Timer::Start()
{
    gettimeofday(&mStart, NULL);
    
    if (mState == kTimerStopped)
    {
        mElapsed = 0.0;
    }
    
    mState = kTimerStarted;
}
//-----------------------------------------------------------------------------
void Timer::Stop()
{
    if (mState != kTimerStarted)
    {
        return;
    }
    
    gettimeofday(&mStop, NULL);
    mElapsed += (mStop.tv_sec - mStart.tv_sec)*1000.0;
    mElapsed += (mStop.tv_usec - mStart.tv_usec)/1000.0; 
    mState = kTimerStopped;
}
//-----------------------------------------------------------------------------
void Timer::Pause()
{
    this->Stop();
    mState = kTimerPaused;
}
//-----------------------------------------------------------------------------
double Timer::GetElapsed() const
{
    return mElapsed;
}
//-----------------------------------------------------------------------------
void Timer::DumpElapsed() const
{
    cout << mElapsed << endl;
}







