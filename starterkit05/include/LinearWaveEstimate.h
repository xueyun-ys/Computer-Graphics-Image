

#ifndef LWE_H
#define LWE_H

#include <string>
#include "fftimgproc.h"
#include "imgproc.h"

namespace img
{
    class LinearWaveEstimate
    {
      public:

        LinearWaveEstimate(const ImgProc& init, const double dispersion_factor);
        ~LinearWaveEstimate(){}

        void ingest(const ImgProc& I);//update A, B etc; take an image and transfor into F space

        //helper functions
        const FFTImgProc& getA() const {return A;}
        const FFTImgProc& getB() const {return B;}

        //use A,B to create a new F space amplitude for image;i,j coordinates in FS, at frame n
        //for multiple channels
        void value(int i, int j, int n, std::vector<std::complex<double> >& amplitude) const;

      protected:

        FFTImgProc A;
        FFTImgProc B;

      private:

        double alpha;
        int frame_count;//N, how far you are

        double dispersion(double kx, double ky) const;
    };

    //img: output
    void extract_image(const LinearWaveEstimate& l, int frame, ImgProc& img);
}

#endif
