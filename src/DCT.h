
#ifndef DCT_H
#define DCT_H

#include <QVector>

#include "../libs/FFTW/fftw3.h"

class DCT
{
public :
    DCT(){};
    ~DCT(){};
    
    bool getDCT(const QVector<float>& array,
                const int dir,
                QVector<float>& transformedArray,
                bool* stop);
    
private :
    bool discreteCosineTransf(const QVector<float>& input,
                              fftw_r2r_kind kind,
                              QVector<float>& transformedArray,
                              bool* stop);

};

#endif