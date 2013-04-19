#include "VideoWriter.h"

//-----------------------------------------------------------------------------
VideoWriter::VideoWriter (const std::string& filename, unsigned int width, 
    unsigned int height)
: mFilename(filename), mWidth(width), mHeight(height)
{
    mWriter = cvCreateVideoWriter(filename.c_str(), CV_FOURCC('D', 'I', 'V', 'X'),
        25, cvSize(mWidth, mHeight), 1);
    mCurrentFrame = cvCreateImage(cvSize(mWidth, mHeight), IPL_DEPTH_8U, 3);
}
//-----------------------------------------------------------------------------
VideoWriter::~VideoWriter ()
{
    cvReleaseVideoWriter(&mWriter);
    cvReleaseImage(&mCurrentFrame);
}
//-----------------------------------------------------------------------------
void VideoWriter::SaveScreenshot (const std::string& filename) const
{
    IplImage* img = cvCreateImage(cvSize(mWidth, mHeight), IPL_DEPTH_8U, 4);
    
    glReadPixels(0, 0, mWidth, mHeight, GL_BGRA, GL_UNSIGNED_BYTE, 
        reinterpret_cast<void*> (img->imageData));
    
    cvFlip(img, img);
    cvSaveImage(filename.c_str(), img);

    cvReleaseImage(&img);       
}
//-----------------------------------------------------------------------------
void VideoWriter::CaptureFrame () const
{
    // read framebuffer in current frame
    glReadPixels(0, 0, mWidth, mHeight, GL_BGR, GL_UNSIGNED_BYTE, 
        reinterpret_cast<void*> (mCurrentFrame->imageData));
    
    // flip current frame
    cvFlip(mCurrentFrame, mCurrentFrame);
    cvWriteFrame(mWriter, mCurrentFrame);       
}
//-----------------------------------------------------------------------------
