#include <string>
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include "OpenGL.h"

class VideoWriter
{
public:
    VideoWriter (const std::string& filename, unsigned int width, 
        unsigned int height);
    ~VideoWriter ();

    void SaveScreenshot (const std::string& filename) const;    
    void CaptureFrame () const;
    void Save (const std::string& filename) const;
    
private:
    std::string mFilename;
    unsigned int mWidth;
    unsigned int mHeight;

    CvVideoWriter* mWriter;
    IplImage* mCurrentFrame;

};
