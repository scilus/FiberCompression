
#include "FWT.h"


/********************************************//**
\brief Compute wavelet filters f and g according to name
\param filter : name of the filter to use (from TransformationType enum)
\param f : f filter
\param g : g filter
***********************************************/
void FWT::compute_wavelet_filters(const TransformationType filter,
                                  std::deque<float>& f,
                                  std::deque<float>& g)
{
    f.clear();
    g.clear();
    
    switch(filter)
    {
        case DB4_TRANSFORMATION :
        {
            f.push_back(.482962913145f);
            f.push_back(.836516303738f);
            f.push_back(.224143868042f);
            f.push_back(-.129409522551f);
            break;
        }
        case DB6_TRANSFORMATION :
        {
            f.push_back(.332670552950f);
            f.push_back(.806891509311f);
            f.push_back(.459877502118f);
            f.push_back(-.135011020010f);
            f.push_back(-.085441273882f);
            f.push_back(.035226291882f);
            break;
        }
        case DB8_TRANSFORMATION :
        {
            f.push_back(.230377813309f);
            f.push_back(.714846570553f);
            f.push_back(.630880767930f);
            f.push_back(-.027983769417f);
            f.push_back(-.187034811719f);
            f.push_back(.030841381836f);
            f.push_back(.032883011667f);
            f.push_back(-.01059740178f);
            break;
        }
        case BIOR_5_3_TRANSFORMATION :
        {
            f.push_back(0.5f);
            f.push_back(0.25f);
            f.push_back(sqrtf(2.0f));
            break;
        }
        case BIOR_9_7_TRANSFORMATION :
        {
            f.push_back(1.586134342f);
            f.push_back(-.05298011854f);
            f.push_back(-.8829110762f);
            f.push_back(.4435068522f);
            f.push_back(1.149604398f);
            break;
        }
        default :
        {
            break;
        }
    }
    
    if(DB4_TRANSFORMATION || DB6_TRANSFORMATION || DB8_TRANSFORMATION)
    {
        // Compute f's norm
        float f_norm = 0.0f;
        for(unsigned int i = 0; i < f.size(); i++)
        {
            f_norm += powf(f[i], 2);
        }
        f_norm = sqrtf(f_norm);
        
        // Normalize f
        for(unsigned int i = 0; i < f.size(); i++)
        {
            f[i] /= f_norm;
        }
        
        if(f.size() % 2 == 0)
        {
            // Add zero at the beginning of f
            f.push_front(0.0f);
                
            // Add zero at the end of g
            g.push_back(0.0f);
            int k = 2;
            for(unsigned int i = f.size() - 1; i >= 1; i--)
            {
                g.push_back(f[i] * powf(-1.0f, k));
                k++;
            }
        }
        else
        {
            int k = 1;
            for(int i = (int)f.size() - 1; i >= 0; i--)
            {
                g.push_back(f[i] * powf(-1.0f, k));
                k++;
            }
        }
    }
}


/********************************************//**
\brief Take an input array and modify it so it became dyadic
\param x : input array
***********************************************/
void FWT::makeDyadic(std::vector<float>& x)
{
    int H = x.size();
    int toAdd = pow(2, ceil(log2((float)H))) - H;
    x.insert(x.end(), toAdd, 0.0f);
}


/********************************************//**
\brief Subsample an input array
\param x : input array
\param start : starting index
\param p : keeps 1 over p samples
\return subsampled array
***********************************************/
std::vector<float> FWT::subsampling(const std::vector<float>& x,
                                    const int start,
                                    const int p)
{
    std::vector<float> y;
    for(unsigned int i = start; i < x.size(); i += p)
    {
        y.push_back(x[i]);
    }
    return y;
}


/********************************************//**
\brief Upsample an input array
\param x : input array
\param p : pad 1 over p samples
\return upsampled array
***********************************************/
std::vector<float> FWT::upsampling(const std::vector<float>& x,
                                   const int p)
{
    std::vector<float> y(p * x.size());
    for(unsigned int i = 0; i < x.size(); i++)
    {
        y[p * i] = x[i];
    }
    return y;
}


/********************************************//**
\brief Circular shift of a 1D array
\param x : input 1D array
\param p : shift by p
\return shifted array
***********************************************/
std::vector<float> FWT::circshift(const std::vector<float>& array,
                                  const int p)
{
    std::vector<float> result = array;
    unsigned int m = array.size();
    for(unsigned int i = 0; i < m; i++)
    {
        int idx = ((i - (p + 1)) % m);
        if(idx < 0) idx += m;
        result[i] = array[idx];
    }
    return result;
}


/********************************************//**
\brief Circular conolution of a 1D array with a 1D filter
\param x : input array (1D)
\param h : filter array (1D)
\return convolved array
***********************************************/
std::vector<float> FWT::cconv(const std::vector<float>& x,
                              const std::deque<float>& h)
{
    int p = h.size();
    int pc = (p+1) / 2;
    if(p % 2 == 0) pc = p / 2;
    
    std::vector<float> y(x.size(), 0.0f);
    for(unsigned int i = 0; i < h.size(); ++i)
    {
        // Circular shift of x
        std::vector<float> circX = circshift(x, i-pc);
        // Multiply circX by h[i]
        std::transform(circX.begin(), circX.end(), circX.begin(), std::bind2nd(std::multiplies<float>(), h[i]));
        // Add result to y
        std::transform(circX.begin(), circX.end(), y.begin(), y.begin(), std::plus<float>());
    }
    return y;
}


/********************************************//**
\brief Remove N points at the end of the input array
\param array : input array (1D)
\param N : Number of points to remove at the end
\return modified array
***********************************************/
void FWT::removePoints(std::vector<float>& array,
                       const int N)
{
    array.resize(array.size() - N);
}


/********************************************//**
\brief Compute Fast Wavelet Transform (FWT) on 1D array
\param array : input array (1D)
\param filter : name of the filter to use (from TransformationType enum)
\param result : resulting transformed array
\param stop : flag indicating if we need to stop the process or not
\return boolean true if the execution is finished,
        false if the process has been stopped during the execution
***********************************************/
bool FWT::fastWaveletTransf(const std::vector<float>& array,
                            const TransformationType filter,
                            std::vector<float>& result,
                            bool* stop)
{
    std::vector<float> fw = array;//(array.begin(), array.end());
    
    // Get filters
    std::deque<float> h, g;
    compute_wavelet_filters(filter, h, g);
    
    // Make sure input is dyadic
    makeDyadic(fw);
    
    // Compute FWT
    int Jmax = log2((float)fw.size()) - 1;
    int Jmin = 0;
    for(int j = Jmax; j >= Jmin; --j)
    {
        if(*stop)
        {
            return false;
        }
        std::vector<float> A(&fw[0], &fw[pow(2.,j+1)]);
        std::vector<float> Coarse = subsampling(cconv(A, h));
        std::vector<float> Detail = subsampling(cconv(A, g));
        std::vector<float> coeff = Coarse;
        coeff.insert(coeff.end(), Detail.begin(), Detail.end());
        
        for(int k = 0; k < pow(2., j+1); ++k)
        {
            if(*stop)
            {
                return false;
            }
            fw[k] = coeff[k];
        }
    }
    result = fw;
    return true;
}


/********************************************//**
\brief Compute Inverse Fast Wavelet Transform (iFWT) on 1D array
\param array : input array (1D)
\param filter : name of the filter to use (from TransformationType enum)
\param origLength : original length of the array so the result is cut at this length
\param result : resulting transformed array
\param stop : flag indicating if we need to stop the process or not
\return boolean true if the execution is finished,
        false if the process has been stopped during the execution
***********************************************/
bool FWT::iFastWaveletTransf(const std::vector<float>& array,
                             const TransformationType filter,
                             const int origLength,
                             std::vector<float>& result,
                             bool* stop)
{
    // Get filters
    std::deque<float> h, g;
    compute_wavelet_filters(filter, h, g);
    
    // Reverse filters
    std::deque<float> hRev = h;
    std::reverse(hRev.begin(), hRev.end());
    std::deque<float> gRev = g;
    std::reverse(gRev.begin(), gRev.end());
    
    std::vector<float> inv(array.begin(), array.end());
    
    // iFWT
    int Jmax = log2((float)inv.size()) - 1;
    int Jmin = 0;
    for(int j = Jmin; j <= Jmax; ++j)
    {
        if(*stop)
        {
            return false;
        }
        std::vector<float> Coarse(&inv[0], &inv[pow(2., j)]);
        std::vector<float> Detail(&inv[pow(2., j)], &inv[pow(2.,j+1)]);
        
        // Get Coarse and Detail
        Coarse = cconv(upsampling(Coarse), hRev);
        Detail = cconv(upsampling(Detail), gRev);
        
        for(int k = 0; k < pow(2., j+1); ++k)
        {
            if(*stop)
            {
                return false;
            }
            inv[k] = Coarse[k] + Detail[k];
        }
        
    }
    
    // Remove points so result is origLength length
    int N = pow(2, ceil(log2((float)origLength))) - origLength;
    removePoints(inv, N);
    result = inv;
    return true;
}


/********************************************//**
\brief Forward Lifting Step for biorthogonal wavelets
\param x : input array (1D)
\param h : h filter
\param m :
\param result : result of forward lifting step
\param stop : flag indicating if we need to stop the process or not
\return boolean true if the execution is finished,
        false if the process has been stopped during the execution
***********************************************/
bool FWT::forward_lifting_step(const std::vector<float>& x,
                               const std::deque<float>& h,
                               const int m,
                               std::vector<float>& result,
                               bool* stop)
{
    // Split x
    std::vector<float> d = subsampling(x, 1);
    std::vector<float> y = subsampling(x, 0);
    
    for(int i = 0; i < m; i++)
    {
        if(*stop)
        {
            return false;
        }
        // Add y[1:end-1 end-1] to y in yTemp
        std::vector<float> yTemp = y;
        for(unsigned int k = 0; k < y.size() - 1; k++)
        {
            if(*stop)
            {
                return false;
            }
            yTemp[k] = y[k] + y[k + 1];
        }
        yTemp.back() = 2 * y.back();
        
        // Multiply yTemp by h[2*i]
        std::transform(yTemp.begin(), yTemp.end(), yTemp.begin(), std::bind2nd(std::multiplies<float>(), h[2*i]));
        
        // d = d - yTemp
        std::transform(d.begin(), d.end(), yTemp.begin(), d.begin(), std::minus<float>());
        
        // Add d[0 0:end-2] to d in dTemp
        std::vector<float> dTemp = d;
        for(unsigned int k = 0; k < d.size(); k++)
        {
            if(*stop)
            {
                return false;
            }
            if(k == 0)
            {
                dTemp[k] = 2 * d[k];
            }
            else
            {
                dTemp[k] = d[k] + d[k - 1];
            }
        }
        
        // Multiply dTemp by h[2*i+1]
        std::transform(dTemp.begin(), dTemp.end(), dTemp.begin(), std::bind2nd(std::multiplies<float>(), h[2*i+1]));
        
        // y = y + dTemp
        std::transform(dTemp.begin(), dTemp.end(), y.begin(), y.begin(), std::plus<float>());
    }
    
    // y = [y*h[end-1]  d/h[end-1]]
    std::vector<float> yTemp, dTemp;
    for(unsigned int k = 0; k < y.size(); k++)
    {
        if(*stop)
        {
            return false;
        }
        yTemp.push_back(y[k] * h.back());
        dTemp.push_back(d[k] / h.back());
    }
    
    y = yTemp;
    y.insert(y.end(), dTemp.begin(), dTemp.end());
    
    result = y;
    return true;
}


/********************************************//**
\brief Backward Lifting Step for biorthogonal wavelets
\param x : input array (1D)
\param h : h filter
\param m :
\param result : result of backward lifting step
\param stop : flag indicating if we need to stop the process or not
\return boolean true if the execution is finished,
        false if the process has been stopped during the execution
***********************************************/
bool FWT::backward_lifting_step(const std::vector<float>& x,
                                const std::deque<float>& h,
                                const int m,
                                std::vector<float>& result,
                                bool* stop)
{
    std::vector<float> y;
    std::vector<float> d;
    
    // Retrieve detail coeffs
    for(unsigned int k = (x.size() -1) / 2 + 1; k < x.size(); k++)
    {
        if(*stop)
        {
            return false;
        }
        d.push_back(x[k] * h.back());
    }
    
    for(unsigned int k = 0; k < x.size() / 2; k++)
    {
        if(*stop)
        {
            return false;
        }
        y.push_back(x[k] / h.back());
    }
    
    for(int i = m - 1; i >= 0; i--)
    {
        if(*stop)
        {
            return false;
        }
        std::vector<float> dTemp = d;
        for(unsigned int k = 0; k < d.size(); k++)
        {
            if(k == 0)
            {
                dTemp[k] = 2 * d[k];
            }
            else
            {
                dTemp[k] = d[k] + d[k - 1];
            }
        }
        
        // Multiply dTemp by h[2*i+1]
        std::transform(dTemp.begin(), dTemp.end(), dTemp.begin(), std::bind2nd(std::multiplies<float>(), h[2*i+1]));
        
        // y = y - dTemp
        std::transform(y.begin(), y.end(), dTemp.begin(), y.begin(), std::minus<float>());
        
        // Add y[1:end-1 end-1] to y in yTemp
        std::vector<float> yTemp = y;
        for(unsigned int k = 0; k < y.size() - 1; k++)
        {
            if(*stop)
            {
                return false;
            }
            yTemp[k] = y[k] + y[k + 1];
        }
        yTemp[yTemp.size() - 1] = 2 * y[y.size() -1];
        
        // Multiply yTemp by h[2*i]
        std::transform(yTemp.begin(), yTemp.end(), yTemp.begin(), std::bind2nd(std::multiplies<float>(), h[2*i]));
        
        // d = d + xTemp
        std::transform(yTemp.begin(), yTemp.end(), d.begin(), d.begin(), std::plus<float>());
    }
    // Merge
    std::vector<float> y1(y.size() * 2, 0.0f);
    int k = 0;
    for(unsigned int i = 0; i < y1.size(); i+=2)
    {
        if(*stop)
        {
            return false;
        }
        y1[i] = y[k];
        y1[i+1] = d[k];
        k++;
    }
    result = y1;
    
    return true;
}


/********************************************//**
\brief Biorthogonal wavelet transformation
\param array : input array (1D)
\param filter : name of the filter to use (from TransformationType enum)
\param dir : direction (1:forward, -1:backward)
\param origLength : original Length useful when dir is -1 to cut result
\return result of lifting step
\param stop : flag indicating if we need to stop the process or not
\param Jmin : minimum step
\return boolean true if the execution is finished,
        false if the process has been stopped during the execution
***********************************************/
bool FWT::perform_wavelet_transform(const std::vector<float>& array,
                                    const TransformationType filter,
                                    const int dir,
                                    const int origLength,
                                    std::vector<float>& result,
                                    bool* stop,
                                    const int Jmin)
{
    std::vector<float> data = array;//(array.begin(), array.end());
    std::vector<float> res;
    
    // Get filters
    std::deque<float> h, g;
    compute_wavelet_filters(filter, h, g);
    
    // Make sure input is dyadic
    makeDyadic(data);
    
    // Number of lifting steps
    int n = data.size();
    int m = (h.size() - 1) / 2;
    int Jmax = log2((float)n) - 1;
    
    switch(dir)
    {
        case 1 :    // Forward
        {
            for(int j = Jmax; j >= Jmin; j--)
            {
                res.clear();
                std::vector<float> x(&data[0], &data[pow(2., j+1)]);
                ;
                if(!forward_lifting_step(x, h, m, res, stop))
                {
                    return false;
                }
                for(int k = 0; k < pow(2., j+1); k++)
                {
                    if(*stop)
                    {
                        return false;
                    }
                    data[k] = res[k];
                }
            }
            break;
        }
        case -1 :   // Backward
        {
            for(int j = Jmin; j <= Jmax; j++)
            {
                std::vector<float> x(&data[0], &data[pow(2., j+1)]);
                if(!backward_lifting_step(x, h, m, res, stop))
                {
                    return false;
                }
                for(int k = 0; k < pow(2., j+1); k++)
                {
                    if(*stop)
                    {
                        return false;
                    }
                    data[k] = res[k];
                }
            }
            
            // Remove N points at the end so it becames the same length as original
            int N = pow(2, ceil(log2((float)origLength))) - origLength;
            removePoints(data, N);
            break;
        }
    }
    result = data;
    return true;
}


/********************************************//**
\brief Discrete Wavelet Transform including Daubechies and Biorthogonal
\param array : input array (1D)
\param filter : name of the filter to use (from TransformationType enum)
\param dir : direction (1:forward, -1:backward)
\param transformedArray : resulting transformed array
\param stop : flag indicating if we need to stop the process or not
\param origLength : original Length useful when dir is -1 to cut result
\return boolean true if the execution is finished,
        false if the process has been stopped during the execution
***********************************************/
bool FWT::getDWT(const QVector<float>& array,
                 const TransformationType filter,
                 const int dir,
                 QVector<float>& transformedArray,
                 bool* stop,
                 const int origLength)
{
    std::vector<float> result;
    if(filter == DB4_TRANSFORMATION || filter == DB6_TRANSFORMATION || filter == DB8_TRANSFORMATION)
    {
        if(dir == 1) return fastWaveletTransf(array.toStdVector(), filter, result, stop);
        else return iFastWaveletTransf(array.toStdVector(), filter, origLength, result, stop);
    }
    else if(filter == BIOR_5_3_TRANSFORMATION || filter == BIOR_9_7_TRANSFORMATION)
    {
        return perform_wavelet_transform(array.toStdVector(), filter, dir, origLength, result, stop);
    }

    transformedArray.fromStdVector(result);
    return true;
}
