//-----------------------------------------------------------------------------
//  Timer.h
//  C(++) Programming Practice
//
//  Created by Arno in Wolde Lübke on 06.01.13.
//-----------------------------------------------------------------------------

#ifndef _TIMER_H
#define _TIMER_H

#include <sys/time.h>

class Timer
{
    enum
    {
        kTimerStarted,
        kTimerStopped,
        kTimerPaused
    };
    
public:
    Timer();
    ~Timer();
    
    void Start();
    void Pause();
    void Stop();
    double GetElapsed() const;
    void DumpElapsed() const;
        
private:
    int mState;
    timeval mStart;
    timeval mStop;
    double mElapsed;
};

#endif
