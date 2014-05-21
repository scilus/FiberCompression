
#include "DCT.h"


/********************************************//**
\brief Compute DCT on a 1D array with ffwt methods
\param input : input array to transform with the DCT
\param kind : direction (FFTW_REDFT10 : forward, FFTW_REDFT01 : backward)
\param transformedArray : resulting transformed array
\param stop : flag indicating if we need to stop the process or not
\return boolean true if the execution is finished,
        false if the process has been stopped during the execution
***********************************************/
bool DCT::discreteCosineTransf(const QVector<float>& input,
                                         fftw_r2r_kind kind,
                                         QVector<float>& transformedArray,
                                         bool* stop)
{
    int n = input.size();
    
    // Allocate resources
    double* in = (double*)fftw_malloc(sizeof(double) * n);
    double* out = (double*)fftw_malloc(sizeof(double) * n);
    
    // Matrix data collection
    for(int i = 0; i < n; i++)
    {
        if(*stop)
        {
            return false;
        }
        in[i] = input[i];
    }
    
    // Create plan
    fftw_plan p;
    p = fftw_plan_r2r_1d(n, in, out, kind, FFTW_ESTIMATE);
    
    if(*stop)
    {
        return false;
    }
    
    // Execute plan
    fftw_execute(p);
    
    if(*stop)
    {
        return false;
    }
    
    // Set output
    QVector<float> output(n, 0.0f);
    for(int i = 0; i < n; i++)
    {
        output[i] = out[i];
        if(kind == FFTW_REDFT01)
        {
            output[i] /= 2.0f * (float)n;
        }
    }
    
    // Free resources
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
    
    transformedArray = output;
    return true;
}


/********************************************//**
\brief Compute DCT on a 1D array
\param array : input array to transform with the DCT
\param dir : direction (1 : forward, -1 : backward)
\param transformedArray : resulting transformed array
\param stop : flag indicating if we need to stop the process or not
\return boolean true if the execution is finished,
        false if the process has been stopped during the execution
***********************************************/
bool DCT::getDCT(const QVector<float>& array,
                 const int dir,
                 QVector<float>& transformedArray,
                 bool* stop)
{
    transformedArray.clear();
    
    if( dir == 1) // Forward
    {
        return discreteCosineTransf(array, FFTW_REDFT10, transformedArray, stop);
    }
    else if( dir == -1) // Backward
    {
        return discreteCosineTransf(array, FFTW_REDFT01, transformedArray, stop);
    }
    
    return true;
};
