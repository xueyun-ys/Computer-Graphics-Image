#include "LinearWaveEstimate.h"
#include "fftimgproc.h"
#include <cmath>

using namespace img;

LinearWaveEstimate::LinearWaveEstimate(const ImgProc& init, const double dispersion_factor) :
  alpha(dispersion_factor),
  frame_count(0)
  {
    A.clear(init.nx(), init.ny(), init.depth());
    B.clear(init.nx(), init.ny(), init.depth());
  }

double LinearWaveEstimate::dispersion(double kx, double ky) const
{
  double kmag = std::sqrt(kx*kx + ky*ky);
  double freq = alpha * std::sqrt(kmag);
  return freq;
}

void LinearWaveEstimate::ingest( const ImgProc& I)
{
  FFTImgProc Itilde;
  img::load_fft_new(I, Itilde);
  Itilde.fft_forward();

  for(int j=0;j<Itilde.ny();j++)
  {
    #pragma omp parallel for
    for(int i=0;i<Itilde.nx();i++)
    {
      std::vector<std::complex<double> > itilde;
      std::vector<std::complex<double> > a;
      std::vector<std::complex<double> > b;
      Itilde.value(i, j, itilde);
      A.value(i, j, a);
      B.value(i, j, b);
      std::vector<std::complex<double> > aupdate = a;
      std::vector<std::complex<double> > bupdate = b;
      std::complex<double> phase(0.0, frame_count * dispersion(Itilde.kx(i),Itilde.ky(j)));
      phase = std::exp(phase);
      double one_over_N = 1.0/(frame_count +1);
      for(size_t c=0;c<itilde.size();c++)
      {
        aupdate[c] += (itilde[c]/phase - b[c]/(phase*phase) -a[c] )*one_over_N;
        bupdate[c] += (itilde[c]*phase - a[c]*(phase*phase) -b[c] )*one_over_N;
      }
      A.set_value(i, j, aupdate);
      B.set_value(i, j, bupdate);
    }
  }
  frame_count++;
}

void LinearWaveEstimate::value(int i, int j, int n, std::vector< std::complex<double> >&amplitude) const
{
  std::complex<double> phase(0.0, n * dispersion(A.kx(i), A.ky(j)));
  phase = std::exp(phase);
  std::vector<std::complex<double> > a;
  std::vector<std::complex<double> > b;
  A.value(i, j, a);
  B.value(i, j, b);
  amplitude.resize(a.size());
  for(size_t c=0;c<a.size();c++)
  {
    amplitude[c] = a[c]*phase + b[c]/phase;
  }
}

void img::extract_image(const LinearWaveEstimate& l, int frame, ImgProc& img)
{
  std::cout<<"successfully get here0"<<std::endl;
  img.clear(l.getA().nx(), l.getA().ny(), l.getA().depth() );
  std::cout<<"successfully get here1"<<std::endl;
  FFTImgProc fftimg;
  fftimg.clear(img.nx(), img.ny(), img.depth() );
  std::cout<<"successfully get here2"<<std::endl;
  for(int j=0;j<img.ny();j++)
  {
    for(int i=0;i<img.nx();i++)
    {
      std::vector<std::complex<double> >v;
      l.value(i, j, frame, v);
      fftimg.set_value(i, j, v);
      //std::cout<<"nx:    "<<img.nx()<<std::endl;
      //std::cout<<"nx:    "<<i<<std::endl;
    }
    std::cout<<"successfully get here33"<<std::endl;
  }
  std::cout<<"successfully get here3"<<std::endl;
  fftimg.fft_backward();
  for(int j=0;j<img.ny();j++)
  {
    for(int i=0;i<img.nx();i++)
    {
      std::vector<std::complex<double> > v;
      fftimg.value(i, j, v);
      std::vector<float> iv(v.size());
      for(size_t c=0;c<v.size();c++)
      {
        iv[c] = v[c].real();
      }
      img.set_value(i, j, iv);
    }
  }
  std::cout<<"successfully get here4"<<std::endl;
}
