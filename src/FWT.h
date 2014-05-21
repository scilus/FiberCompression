
#ifndef FWT_H
#define FWT_H

#include <deque>
#include <vector>

#include <QString>

#include "Utils.h"


class FWT
{
public:
    FWT(){};
    ~FWT(){};
    
    // General Wavelet Method
    bool getDWT(const QVector<float>& array,
                const TransformationType name,
                const int dir,
                QVector<float>& result,
                bool* stop,
                const int origLength = 0);
    
private:
    void compute_wavelet_filters(TransformationType filter,
                             std::deque<float>& f,
                             std::deque<float>& g);

    void makeDyadic(std::vector<float>& x);

    std::vector<float> subsampling(const std::vector<float>& x,
                                   const int start = 0,
                                   const int p = 2);

    std::vector<float> upsampling(const std::vector<float>& x,
                                  const int p = 2);

    std::vector<float> circshift(const std::vector<float>& array,
                                 const int p);


    std::vector<float> cconv(const std::vector<float>& x,
                             const std::deque<float>& h);

    void removePoints(std::vector<float>& array,
                      const int N);


    // Forward Wavelet Transform
    bool fastWaveletTransf(const std::vector<float>& array,
                           const TransformationType filter,
                           std::vector<float>& result,
                           bool* stop);
    bool iFastWaveletTransf(const std::vector<float>& array,
                            const TransformationType filter,
                            const int origLength,
                            std::vector<float>& result,
                            bool* stop);


    bool forward_lifting_step(const std::vector<float>& x,
                              const std::deque<float>& h,
                              const int m,
                              std::vector<float>& result,
                              bool* stop);

    bool backward_lifting_step(const std::vector<float>& x,
                               const std::deque<float>& h,
                               const int m,
                               std::vector<float>& result,
                               bool* stop);

    // Wavelet transform (Biorthogonal)
    bool perform_wavelet_transform(const std::vector<float>& array,
                                   const TransformationType filter,
                                   const int dir,
                                   const int origLength,
                                   std::vector<float>& result,
                                   bool* stop,
                                   const int Jmin = 1);
};






                                        
#endif