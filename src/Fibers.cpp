#include <assert.h>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdio.h>

#include <QFile>
#include <QPair>
#include <QStringList>

#include "DCT.h"
#include "Encoding.h"
#include "Fibers.h"
#include "FWT.h"

#include "../libs/nifti/nifti1_io.h"

bool sortDescending(float a, float b)
{
    return a > b;
}


/********************************************//**
* \brief Default Constructor
***********************************************/
Fibers::Fibers() :
m_colorArrayLoadedFromFile(false),
m_countLines(0),
m_fiberType(DEFAULT_FIBER_TYPE),
m_globalMin(std::numeric_limits<float>::max()),
m_globalMax(0.0f)
{
    m_pStop = new bool(false);
}


/********************************************//**
* \brief Default Destructor
***********************************************/
Fibers::~Fibers()
{
    if(m_pStop != NULL)
    {
        delete m_pStop;
        m_pStop = NULL;
    }
}


/********************************************//**
\brief Initialization of all parameters
\param compress : boolean that indicates if we compress or not
\param inputPath : input file to compress
\param outputPath : compressed output file
\param outputStatsFilename : text file where to output statistics during process
\param errorMax : maximum error threshold for linearization and compression (in mm)
\param transformationType : transformation type
\param quantizationType : quantization type
\param precisionQuantization : precision of the rounding (10^n position for uniform and n most
significant digits for non-uniform)
\param encodingType : encoding type
***********************************************/
void Fibers::init(const bool compress,
                  const QString& inputPath,
                  const QString& outputPath,
                  const QString& outputStatsFilename,
                  const float errorMax,
                  const TransformationType transformationType,
                  const QuantizationType quantizationType,
                  const int precisionQuantization,
                  const EncodingType encodingType)
{
    m_compress = compress;
    m_inputPath = inputPath;
    m_outputPath = outputPath;
    m_outputStatsFilename = outputStatsFilename;
    m_errorMax = errorMax;
    
    bool outputExtSameAsInput = (getExtension(m_inputPath) == getExtension(m_outputPath));
    
    m_transformationType = transformationType;
    m_quantizationType = quantizationType;
    m_precisionQuantization = precisionQuantization;
    m_encodingType = encodingType;
    
    if(outputExtSameAsInput)
    {
        m_transformationType = NO_TRANSFORMATION;
        m_quantizationType = NO_QUANTIZATION;
        m_encodingType = NO_ENCODING;
    }
    
    if(m_errorMax < 0.2)
    {
        m_precisionQuantization = precisionQuantization;
    }
    else
    {
        m_precisionQuantization = -1;
    }
}


/********************************************//**
\brief Clean up pointers
***********************************************/
void Fibers::cleanUp()
{
    if(m_pStop != NULL)
    {
        delete m_pStop;
        m_pStop = NULL;
    }
}


/********************************************//**
\brief Reset member variables to default values
***********************************************/
void Fibers::reset()
{
    m_colorArrayLoadedFromFile = false;
    m_countLines = 0;

    m_pOriginal = NULL;
    
    m_pStop = new bool(false);
    
    m_outputStatsFilename = "";
    m_inputPath = "";
    m_outputPath = "";
    m_errorMax = DEFAULT_ERROR_MAX;
    m_transformationType = DEFAULT_TRANSF_TYPE;
    m_quantizationType = DEFAULT_QUANT_TYPE;
    m_precisionQuantization = DEFAULT_PRECISION_QUANT;
    m_encodingType = DEFAULT_ENCODING_TYPE;
    m_compress = false;
    
    m_fiberType = DEFAULT_FIBER_TYPE;
    m_huffmanDict.cleanUp();
    m_colorHuffmanDict.cleanUp();
    m_arithmeticDict.cleanUp();
    m_colorArithmeticDict.cleanUp();
    
    //m_pointsArray.clear();
    m_rColorArray.clear();
    m_gColorArray.clear();
    m_bColorArray.clear();
    m_xPointArray.clear();
    m_yPointArray.clear();
    m_zPointArray.clear();
    m_origLineLength.clear();
    m_xEncoded = "";
    m_yEncoded = "";
    m_zEncoded = "";
    m_rEncodedColor = "";
    m_gEncodedColor = "";
    m_bEncodedColor = "";
}


/********************************************//**
\brief Process compression or decompression pipeline
***********************************************/
void Fibers::process()
{
    if(m_compress)
    {
        if(!compress())
        {
            emit isStopped();
        }
    }
    else
    {
        if(!decompress())
        {
            emit isStopped();
        }
    }
}


/********************************************//**
\brief Stop current process
***********************************************/
void Fibers::stop()
{
    (*m_pStop) = true;
}


/********************************************//**
\brief Compression pipeline
\return boolean true if everything went well, false otherwise
***********************************************/
bool Fibers::compress()
{
    // 1. Load fibers file
    if(!load(m_inputPath))
    {
        return false;
    }
    
    // 2. Keep a copy of the original fiber (suitable only to compute errors in transformation step)
    if(m_transformationType != NO_TRANSFORMATION)
    {
        m_pOriginal = this;
    }
    
    // 3. Get Quantization maximum error
    displayMessage(QString("Linearization..."));
    float linError = m_errorMax;
    if(m_quantizationType != NO_QUANTIZATION) // If we do not quantize, use full error max
    {
        float qMax = getQuantizationMaxError(m_quantizationType, m_precisionQuantization);
        linError = m_errorMax - qMax;
    }
    
    // 4. Linearization
    saveResults("Linearization maximum error : " + QString::number(linError) + "\n");
    if(linError > 0.0f)
    {
        if(!linearize(linError))
        {
            return false;
        }
    }
    
    // 5. Approximation
    if(m_transformationType != NO_TRANSFORMATION)
    {
        displayMessage(QString("Transformation..."));
        if(!approximate(m_transformationType, m_errorMax))
        {
            return false;
        }
        
    }
    
    // 6. Quantization
    if(m_quantizationType != NO_QUANTIZATION)
    {
        displayMessage(QString("Quantization..."));
        if(!quantization(m_xPointArray, m_quantizationType, m_precisionQuantization))
        {
            return false;
        }
        if(!quantization(m_yPointArray, m_quantizationType, m_precisionQuantization))
        {
            return false;
        }
        if(!quantization(m_zPointArray, m_quantizationType, m_precisionQuantization))
        {
            return false;
        }
    }
    
    // Encode or save directly in original format
    if(m_encodingType != NO_ENCODING)
    {
        // 7. Encode
        displayMessage(QString("Encoding..."));
        if(!encode(m_xPointArray, m_yPointArray, m_zPointArray, m_encodingType))
        {
            return false;
        }
        
        // 8. Save compressed file
        displayMessage(QString("Saving..."));
        if(!saveCompressed(m_encodingType, m_transformationType, m_outputPath))
        {
            return false;
        }
    }
    else
    {
        // 7. Save with original format
        displayMessage(QString("Saving..."));
        if(!save(m_outputPath))
        {
            return false;
        }
    }
    
    #if _SAVE_STATS
        float cr = getCompressionRatio();
        displayMessage("Compression Ratio : " + QString::number(cr));
    #endif
    
    emit finished();
    
    return true;
}


/********************************************//**
\brief Decompression pipeline
\return boolean true if everything went well, false otherwise
***********************************************/
bool Fibers::decompress()
{
    m_encodingType = NO_ENCODING;
    m_transformationType = NO_TRANSFORMATION;
    
    displayMessage(QString("Loading..."));
    
    // 1. Load fibers
    if(!loadCompressed(m_inputPath, m_encodingType, m_transformationType))
    {
        return false;
    }
    
    // 2. Validate fibertype and output extension
    bool match = false;
    QString fiberTypeExt;
    QString outputExt = getExtension(m_outputPath);
    switch(m_fiberType)
    {
        case FIBERTYPE_MRTRIX :
            fiberTypeExt = "tck";
            match = (outputExt == fiberTypeExt);
            break;
        case FIBERTYPE_TRK :
            fiberTypeExt = "trk";
            match = (outputExt == fiberTypeExt);
            break;
        default :
            break;

    }
    if(!match)
    {
        showMessage(QString("Output file extension provided (."
                            + outputExt + ") doesn't match the input file extension used to compress this file (."
                            + fiberTypeExt + ")."));
        return false;
    }
    
    // 3. Decode
    displayMessage(QString("Decoding..."));
    if(!decode(m_encodingType))
    {
        return false;
    }
    
    // 4. Inverse transformation on each dimension
    if(m_transformationType != NO_TRANSFORMATION)
    {
        displayMessage(QString("Inverse Transformation..."));
        QList< QVector<float> > transformed;
        if(!getTransformation(m_transformationType, m_xPointArray, -1, transformed))
        {
            return false;
        }
        setPointsArrayDim(transformed, 0);
        if(!getTransformation(m_transformationType, m_yPointArray, -1, transformed))
        {
            return false;
        }
        setPointsArrayDim(transformed, 1);
        if(!getTransformation(m_transformationType, m_zPointArray, -1, transformed))
        {
            return false;
        }
        setPointsArrayDim(transformed, 2);
    }
    
    // 5. Save uncompressed file
    displayMessage(QString("Saving..."));
    if(!save(m_outputPath))
    {
        return false;
    }
    
    emit finished();
    return true;
}


/********************************************//**
\brief display a message
\param message : message to display
***********************************************/
void Fibers::displayMessage(const QString& message)
{
    emit showMessage(message);
}


/********************************************//**
\brief Linearization of the entire dataset
\param thresholdValue : maximum error threshold in mm
\return boolean true if the execution is finished, 
        false if the process has been stopped during the execution
***********************************************/
bool Fibers::linearize(const double thresholdValue)
{
    int countPoints = 0;
    assert(thresholdValue >= 0.0);
    
    QVector<float> linearizedFiberX, linearizedFiberY, linearizedFiberZ;
    QVector<unsigned char> colorFiberX, colorFiberY, colorFiberZ;
    
    // Linearize each fiber separately
    for (int f = 0; f < getCountLines(); f++)
    {
        if(!linearizeFiber(f, thresholdValue, colorFiberX, colorFiberY, colorFiberZ, linearizedFiberX,
                           linearizedFiberY, linearizedFiberZ))
        {
            return false;
        }
        
        setFiberInPointArray(linearizedFiberX, f, 0);
        setFiberInPointArray(linearizedFiberY, f, 1);
        setFiberInPointArray(linearizedFiberZ, f, 2);
        countPoints += linearizedFiberX.size();
        
        if(m_colorArrayLoadedFromFile)
        {
            setColorArray(colorFiberX, f, 0);
            setColorArray(colorFiberY, f, 1);
            setColorArray(colorFiberZ, f, 2);
        }
    }
    saveResults("Total number of points after linearization : " + QString::number(countPoints) + "\n");
    
    return true;
}


/********************************************//**
\brief Linearization of a single fiber
\param fibidx : index of the current fiber to linearize
\param thresholdValue : maximum error threshold in mm
\param colorArray : color array that corresponds to the fiber to update at the same time if necessary
\param fiberResult : the linearized fiber
\return boolean true if the execution is finished,
        false if the process has been stopped during the execution
***********************************************/
bool Fibers::linearizeFiber(const int fibIdx,
                            const double thresholdValue,
                            QVector<unsigned char>& rColorArray,
                            QVector<unsigned char>& gColorArray,
                            QVector<unsigned char>& bColorArray,
                            QVector<float>& fiberResultX,
                            QVector<float>& fiberResultY,
                            QVector<float>& fiberResultZ)
{
    Vector Q1, Q2, P;
    fiberResultX.resize(0);
    fiberResultY.resize(0);
    fiberResultZ.resize(0);
    rColorArray.resize(0);
    gColorArray.resize(0);
    bColorArray.resize(0);
    
    int pointIdx = 0;
    int segmentSecondPtIdx = 2;
    
    Q1 = Vector(m_xPointArray[fibIdx][0], m_xPointArray[fibIdx][0], m_zPointArray[fibIdx][0]);
    fiberResultX.push_back(m_xPointArray[fibIdx][0]);
    fiberResultY.push_back(m_yPointArray[fibIdx][0]);
    fiberResultZ.push_back(m_zPointArray[fibIdx][0]);
    
    while(segmentSecondPtIdx < getLineSize(fibIdx))
    {
        if(*m_pStop)
        {
            return false;
        }
        float distance = 0.0f;
        
        // Ignore as many points as possible according to error threshold
        while((distance < thresholdValue) && segmentSecondPtIdx < getLineSize(fibIdx))
        {
            if( pointIdx == 0)
            {
                pointIdx = 1;
            }
            Q2 = Vector(m_xPointArray[fibIdx][segmentSecondPtIdx],
                        m_xPointArray[fibIdx][segmentSecondPtIdx],
                        m_zPointArray[fibIdx][segmentSecondPtIdx]);
            
            float d = 0.0f;
            for (int k = pointIdx; k < segmentSecondPtIdx; k++)
            {
                P = Vector(m_xPointArray[fibIdx][k], m_xPointArray[fibIdx][k], m_zPointArray[fibIdx][k]);
                d = getPointToSegmentDistance(P, Q1, Q2);
                distance = std::max(distance, d);
            }
            if(distance < thresholdValue)
            {
                segmentSecondPtIdx++;
            }
        }
        // Add the last point that can't be removed
        pointIdx = segmentSecondPtIdx;
        Q1 = Vector(m_xPointArray[fibIdx][segmentSecondPtIdx-1],
                    m_xPointArray[fibIdx][segmentSecondPtIdx-1],
                    m_zPointArray[fibIdx][segmentSecondPtIdx-1]);
        fiberResultX.push_back(m_xPointArray[fibIdx][segmentSecondPtIdx-1]);
        fiberResultY.push_back(m_yPointArray[fibIdx][segmentSecondPtIdx-1]);
        fiberResultZ.push_back(m_zPointArray[fibIdx][segmentSecondPtIdx-1]);
        
        // Make sure that last point will always be added if necessary
        if(segmentSecondPtIdx == getLineSize(fibIdx) - 1)
        {
            Q1 = Vector(m_xPointArray[fibIdx][segmentSecondPtIdx],
                        m_xPointArray[fibIdx][segmentSecondPtIdx],
                        m_zPointArray[fibIdx][segmentSecondPtIdx]);
            fiberResultX.push_back(m_xPointArray[fibIdx][segmentSecondPtIdx]);
            fiberResultY.push_back(m_yPointArray[fibIdx][segmentSecondPtIdx]);
            fiberResultZ.push_back(m_zPointArray[fibIdx][segmentSecondPtIdx]);
        }
        
        // Sync color array if necessary
        if(m_colorArrayLoadedFromFile)
        {
            rColorArray.push_back(m_rColorArray[fibIdx][segmentSecondPtIdx-1] * 255);
            gColorArray.push_back(m_gColorArray[fibIdx][segmentSecondPtIdx-1] * 255);
            bColorArray.push_back(m_bColorArray[fibIdx][segmentSecondPtIdx-1] * 255);
        }
        distance = 0.0f;
        segmentSecondPtIdx++;
    }
    return true;
}


/********************************************//**
\brief Transformation and non-linear approximation of an entire dataset
\param ttype : transformation type to use (0:DCT, 1:db4, 2:db6, 3:5-3, 4:9-7)
\param errorMaxT : maximum error threshold in mm
\param xPoints : list of fibers in x dimension approximated with the transformation
\param yPoints : list of fibers in y dimension approximated with the transformation
\param zPoints : list of fibers in z dimension approximated with the transformation
\param originalArray original dataset to compare with when findind the best threshold in non-linear
approximation
\return boolean true if the execution is finished,
        false if the process has been stopped during the execution
***********************************************/
bool Fibers::approximate(TransformationType ttype,
                         const float errorMaxT)
{
    // Make sure dimension fits
    assert(m_xPointArray.size() > 0);
    assert(m_xPointArray.size() == m_yPointArray.size() &&
           m_xPointArray.size() == m_zPointArray.size());
    
    std::vector<int> cpt;
    std::vector<float> histo;
    std::vector<float> std;
    
    // Treat each fiber separately
    for(int k = 0; k < m_xPointArray.size(); k++)
    {
        if(*m_pStop)
        {
            return false;
        }
        int size = m_xPointArray[k].size();
        assert(size > 0);
        
        // Compute forward transformations on each dimension of the current fiber
        QVector<float> transf_x;
        if(!getArrayTransformation(ttype, m_xPointArray[k], 1, size, transf_x))
        {
            return false;
        }
        QVector<float> transf_y;
        if(!getArrayTransformation(ttype, m_yPointArray[k], 1, size, transf_y))
        {
            return false;
        }
        QVector<float> transf_z;
        if(!getArrayTransformation(ttype, m_zPointArray[k], 1, size, transf_z))
        {
            return false;
        }
        assert(transf_x.size() == transf_y.size() && transf_x.size() == transf_z.size());
        
        
        // Find number of coefficients for non-linear approximation according to errorMaxT
        int lower = 1;
        int upper = transf_x.size();
        int m = lower + ((upper - lower) / 2);
        int res_m = m;
        
        bool notFound = true;
        // Find best m
        while(notFound)
        {
            // Check flag to know if we need to stop the process or not
            if(*m_pStop)
            {
                return false;
            }
            m = res_m;
            
            // Approximate with m coefficients
            QVector<float> aTransf_x = nonLinearApproximation(transf_x, m);
            QVector<float> aTransf_y = nonLinearApproximation(transf_y, m);
            QVector<float> aTransf_z = nonLinearApproximation(transf_z, m);
            
            // Get backward transformations
            if(!getArrayTransformation(ttype, aTransf_x , -1, size, aTransf_x))
            {
                return false;
            }
            if(!getArrayTransformation(ttype, aTransf_y , -1, size, aTransf_y))
            {
                return false;
            }
            if(!getArrayTransformation(ttype, aTransf_z , -1, size, aTransf_z))
            {
                return false;
            }
            
            // Compute error
            float errorMax, errorMean, errorStd;
            QVector<Vector> mergeVec = mergeDimensions(aTransf_x,
                                                       aTransf_y,
                                                       aTransf_z);
            QVector<Vector> origMergeVec = mergeDimensions(m_xPointArray[k],
                                                           m_yPointArray[k],
                                                           m_zPointArray[k]);
            if(!errorsBetweenFibers(origMergeVec, mergeVec, errorMax, errorMean, errorStd))
            {
                return false;
            }
            
            // Update bounds
            if(errorMax > errorMaxT)
            {
                lower = m;
            }
            else
            {
                upper = m;
            }
            res_m = lower + ((upper - lower) / 2);
            if(upper - lower == 1)
            {
                notFound = false;
                res_m = upper;
            }
        }
    
        // Approximate with best m
        if(res_m < transf_x.size())
        {
            transf_x = nonLinearApproximation(transf_x, res_m);
            transf_y = nonLinearApproximation(transf_y, res_m);
            transf_z = nonLinearApproximation(transf_z, res_m);
        }
        
        // Save the result
        setFiberInPointArray(transf_x, k, 0);
        setFiberInPointArray(transf_y, k, 1);
        setFiberInPointArray(transf_z, k, 2);
    }
    
    for(int k = 0; k < m_xPointArray.size(); k++)
    {
        int size = m_xPointArray[k].size();
        assert(size > 0);
        
        // Compute forward transformations on each dimension of the current fiber
        QVector<float> transf_x;
        if(!getArrayTransformation(ttype, m_xPointArray[k], 1, size, transf_x))
        {
            return false;
        }
        QVector<float> transf_y;
        if(!getArrayTransformation(ttype, m_yPointArray[k], 1, size, transf_y))
        {
            return false;
        }
        QVector<float> transf_z;
        if(!getArrayTransformation(ttype, m_zPointArray[k], 1, size, transf_z))
        {
            return false;
        }
        
    }
    
    return true;
}


/********************************************//**
\brief Transformation of several 1D arrays (list of fibers in x, y or z dimension)
\param type : transformation type to use (0:DCT, 1:db4, 2:db6, 3:5-3, 4:9-7)
\param array : dataset to transform
\param dir : direction of the transformation (1:forward, -1:backward)
\param transformed : input array transformed according to transformtion type
\return boolean true if the execution is finished,
        false if the process has been stopped during the execution
***********************************************/
bool Fibers::getTransformation(const TransformationType type,
                               const QList< QVector<float> >& array,
                               const int dir,
                               QList< QVector<float> >& transformed)
{
    transformed.clear();
    
    // Get transformation on each fiber separately
    for(int f = 0; f < array.size(); f++)
    {
        if(*m_pStop)
        {
            return false;
        }
        QVector<float> transfresult;
        if(!getArrayTransformation(type, array[f], dir, m_xPointArray[f].size(), transfresult))
        {
            return false;
        }
        transformed.push_back(transfresult);
    }
    return true;
}


/********************************************//**
\brief Apply a transformation on a 1D array
\param type : transformation type to use (0:DCT, 1:db4, 2:db6, 3:5-3, 4:9-7)
\param array : 1D array to transform
\param dir : direction of the transformation (1:forward, -1:backward)
\param origLength : originalLength of the array. This is useful for inverse transformation of wavelets
(type 1,2,3 or 4 and dir -1) because of the dyadic padding to make sure that the
inverse transformation have the same length than the original array.
\param transformedArray : array transformed with transformation type
\return boolean true if the execution is finished,
        false if the process has been stopped during the execution
***********************************************/
bool Fibers::getArrayTransformation(const TransformationType type,
                                    const QVector<float>& array,
                                    const int dir,
                                    const int origLength,
                                    QVector<float>& transformedArray)
{
    transformedArray.clear();
    
    // Transform the 1D array according to the transformation type and the direction
    switch(type)
    {
        case DCT_TRANSFORMATION :
        {
            DCT dct;
            return dct.getDCT(array, dir, transformedArray, m_pStop);
            break;
        }
        case DB4_TRANSFORMATION || DB6_TRANSFORMATION || DB8_TRANSFORMATION ||
            BIOR_5_3_TRANSFORMATION || BIOR_9_7_TRANSFORMATION:
        {
            FWT wavelets;
            return wavelets.getDWT(array, type, dir, transformedArray, m_pStop, origLength);
            break;
        }
        default:
            return false;
    }
    return true;
}


/********************************************//**
\brief Non-linear approximation of a 1D array
\param array : 1D array to approximate
\param nbKeptCoeffs : number of coefficients to keep in the approximation
\return input array approximated with nbKeptCoeffs coefficients
***********************************************/
QVector<float> Fibers::nonLinearApproximation(const QVector<float>& array,
                                              const int nbKeptCoeffs)
{
    assert(nbKeptCoeffs > 0);
    assert(array.size() > 0);
    
    QVector<float> result = array;
    QVector<float> toSort;
    int size = (int)result.size();
    
    // Find threshold corresponding to nbKeptCoeffs
    for(int i = 0; i < size; i++)
    {
        toSort.push_back(fabs(result[i]));
    }
    std::sort(toSort.begin(), toSort.end(), sortDescending);
    float threshold = toSort[std::min(nbKeptCoeffs, size) - 1];
    
    // Apply threshold to array
    for(int i = 0; i < size; i++)
    {
        if(fabs(result[i]) < threshold)
        {
            result[i] = 0.0f;
        }
    }
    return result;
}


/********************************************//**
\brief compute quantization maximum error
\param quantizationType : rounding type (0:uniform, 1:non-uniform)
\param precision : precision of the rounding (10^n position for uniform
and n most significant digits for non-uniform)
\return maximum quantization error
***********************************************/
float Fibers::getQuantizationMaxError(QuantizationType type, int precision)
{
    float maxError = 0.0f;
    
    if(type == UNIFORM_QUANTIZATION)
    {
        maxError = sqrt(3 * powf(powf(10, precision), 2));
    }
    else
    {
        float max_num = std::max(fabs(m_globalMin), fabs(m_globalMax));
        
        // Find number of digits at the left of .
        QString num = QString::number(fabs(max_num));
        QStringList list = num.split(QChar('.'));
        
        int i = 0;
        int p = -9999999;
        while(i < list[0].size())
        {
            if(list[0][i] != '0')
            {
                p = list[0].size() - i - 1;
                break;
            }
            i++;
        }
        if(p == -9999999)
        {
            i = 0;
            while(i < list[1].size())
            {
                if(list[1][i] != '0')
                {
                    p = i - 1;
                    break;
                }
                i++;
            }
        }
        p = p - precision + 1;
        maxError = sqrt(3 * powf(powf(10, p), 2));
    }
    saveResults("Quantization maximum error : " + QString::number(maxError) + "\n");
    return maxError;
}


/********************************************//**
\brief Rounding of several 1D arrays (list of fibers in x, y or z dimension)
\param array : list of 1D arrays to round
\param type : rounding type (0:uniform, 1:non-uniform)
\param p : precision of the rounding (10^n position for uniform and n most significant digits for non-
uniform)
\return boolean true if the execution is finished,
        false if the process has been stopped during the execution
***********************************************/
bool Fibers::quantization(QList< QVector<float> >& array,
                          const QuantizationType type,
                          const int p)
{
    // Apply rounding according to quantization type and precision p
    switch(type)
    {
        case UNIFORM_QUANTIZATION :
            return uniformQuantization(array, p);
        case NON_UNIFORM_QUANTIZATION :
            return nonUniformQuantization(array, p);
        default:
            std::cout << "Invalid quantization type. Must be 0 for uniform quantization or 1 for non-uniform quantization" << std::endl;
            return false;
    }
}


/********************************************//**
\brief Uniform rounding of several 1D arrays (list of fibers in x, y or z dimension)
\param array : list of 1D arrays to round
\param p : precision of the rounding (at 10^n position)
\return boolean true if the execution is finished,
        false if the process has been stopped during the execution
***********************************************/
bool Fibers::uniformQuantization(QList< QVector<float> >& array,
                                 const int p)
{
    int num = 0;
    for (int i = 0; i < array.size(); i++)
    {
        // Check flag to know if we need to stop the process or not
        if(*m_pStop)
        {
            return false;
        }
        
        for (int j = 0; j < array[i].size(); j++)
        {
            float value = array[i][j];
            
            // Find number of digits at the left of .
            int length = QString::number(int(fabs(value))).size();
            
            // Find corresponding number of significant digits to use non-uniform quantization
            if(p >= 0)
            {
                int m = QString::number(int(pow(10.0, p))).size();
                num = length - m + 1;
            }
            else
            {
                num = length + abs(p);
            }
            
            // Use non-uniform quantization
            float result = nonUniformRounding(value, num);
            array[i][j] = result;
        }
    }
    return true;
}


/********************************************//**
\brief Non-uniform rounding of several 1D arrays (list of fibers in x, y or z dimension)
\param array : list of 1D arrays to round
\param p : precision of the rounding (with p most significant digits)
\return boolean true if the execution is finished,
        false if the process has been stopped during the execution
***********************************************/
bool Fibers::nonUniformQuantization(QList< QVector<float> >& array,
                                    const int p)
{
    assert(p > 0);
    
    for (int i = 0; i < array.size(); i++)
    {
        // Check flag to know if we need to stop the process or not
        if(*m_pStop)
        {
            return false;
        }
        
        for (int j = 0; j < array[i].size(); j++)
        {
            float value = nonUniformRounding(array[i][j], p);
            array[i][j] = value;
        }
    }
    return true;
}



/********************************************//**
\brief Non-uniform rounding of a floating point value
\param value : value to round
\param p : precision of the rounding (with p most significant digits)
\return rounded value
***********************************************/
float Fibers::nonUniformRounding(const float value,
                                 const int p)
{
    // Find the power
    if(value == 0.0f) return value;
    float d = ceil(log10(value < 0 ? -value : value));
    int power = p - (int) d;
    
    // Find the magnitude
    float magnitude = powf(10, power);
    assert(magnitude != std::numeric_limits<float>::infinity() && "Precision is too high !");
    long shifted = floor((value * magnitude) + 0.5);
    
    // Return the rounded value
    return shifted/magnitude;
}


/********************************************//**
\brief compute max, mean, mean-max, mean-std and max-std errors between an entire dataset and an
original dataset
\param original : original instance useful to compare with the original dataset
\param xPointsArray : list of fibers in x dimension to compute errors with
\param yPointsArray : list of fibers in y dimension to compute errors with
\param zPointsArray : list of fibers in z dimension to compute errors with
\param maxError : maximum error in mm
\param meanError : mean error in mm
\param meanMaxError : mean of maximum errors in mm
\param meanStdError : mean standard deviation error in mm
\param maxStdError : maximum standard deviation error in mm
\return boolean true if the execution is finished,
        false if the process has been stopped during the execution
***********************************************/
bool Fibers::getErrors(Fibers* original,
                       const QList< QVector<float> >& xPointsArray,
                       const QList< QVector<float> >& yPointsArray,
                       const QList< QVector<float> >& zPointsArray,
                       float& maxError,
                       float& meanError,
                       float& meanMaxError,
                       float& meanStdError,
                       float& maxStdError)
{
    maxError = meanError = meanMaxError = meanStdError = maxStdError = 0.0f;
    
    // Make sure original and final fibers have the same number of fibers
    assert(original->m_xPointArray.size() == xPointsArray.size());
    assert(original->m_yPointArray.size() == yPointsArray.size());
    assert(original->m_zPointArray.size() == zPointsArray.size());
    
    for( int f = 0; f < xPointsArray.size(); f++)
    {
        if(*m_pStop)
        {
            return false;
        }
        // Merge dimensions to get final fiber
        QVector<Vector> finalFiber = mergeDimensions(xPointsArray[f],
                                                     yPointsArray[f],
                                                     zPointsArray[f]);
        QVector<Vector> origFiber = mergeDimensions(original->m_xPointArray[f],
                                                    original->m_yPointArray[f],
                                                    original->m_zPointArray[f]);
        
        // Get errors between original and final fibers
        float currentMaxError, currentMeanError, currentStdError;
        if(!errorsBetweenFibers(origFiber,
                                finalFiber,
                                currentMaxError,
                                currentMeanError,
                                currentStdError))
        {
            return false;
        }
        
        // Update global errors from current error
        maxError = std::max(maxError, currentMaxError);
        meanMaxError += maxError;
        meanError += currentMeanError;
        meanStdError += currentStdError;
        maxStdError = std::max(maxStdError, currentStdError);
    }
    meanMaxError /= (float)xPointsArray.size();
    meanError /= (float)xPointsArray.size();
    meanStdError /= (float)xPointsArray.size();
    
    return true;
}


/********************************************//**
\brief compute max, mean, mean-max, mean-std and max-std errors between a modified and an original fiber
\param original : original fiber (1D array of 3D points)
\param final : modified fiber (1D array of 3D points)
\param maxError : maximum error in mm
\param meanError : mean error in mm
\param stdError : standard deviation error in mm
\return boolean true if the execution is finished,
        false if the process has been stopped during the execution
***********************************************/
bool Fibers::errorsBetweenFibers(const QVector<Vector>& original,
                                 const QVector<Vector>& final,
                                 float& maxError,
                                 float& meanError,
                                 float& stdError)
{
    QVector<float> dist;
    maxError = meanError = stdError = 0.0f;
    
    // Find the distance between each original points with each segments of the final fiber
    for(int j = 0; j < original.size(); j++)
    {
        float dist_min_pt = std::numeric_limits<float>::max();
        
        // Find the closest final segment of original point
        for( int i = 0; i < final.size() - 1; i++)
        {
            if(*m_pStop)
            {
                return false;
            }
            Vector P = original[j];     // original point
            Vector P0 = final[i];       // first point of final segment
            Vector P1 = final[i+1];     // second point of final segment
            
            // Get the distance between original point (P) and final segment (P0P1)
            float d = getPointToSegmentDistance(P, P0, P1);
            dist_min_pt = std::min(d, dist_min_pt);
        }
        // The maximum error for the original point will be the max distance with all final segments
        dist.push_back(dist_min_pt);
        maxError = std::max(dist_min_pt, maxError);
        meanError += dist_min_pt;
    }
    
    if(*m_pStop)
    {
        return false;
    }
    
    // Update mean and std errors
    meanError /= (float)original.size();
    
    for(int j = 0; j < dist.size(); j++)
    {
        stdError += powf(dist[j] - meanError, 2);
    }
    stdError /= (float)dist.size();

    return true;
}



/********************************************//**
\brief encode arrays in 0,1 bytes
\param xArray : list of fibers in x dimension to encode
\param yArray : list of fibers in y dimension to encode
\param zArray : list of fibers in z dimension to encode
\encodingType encoding type (0:huffman, 1:arithmetic)
\return boolean true if the execution is finished,
        false if the process has been stopped during the execution
***********************************************/
bool Fibers::encode(const QList< QVector<float> >& xArray,
                    const QList< QVector<float> >& yArray,
                    const QList< QVector<float> >& zArray,
                    const EncodingType encodingType)
{
    assert((xArray.size() == yArray.size()) && (xArray.size() == zArray.size()));
    
    // Encode according to encodingType
    switch(encodingType)
    {
        case HUFFMAN_ENCODING :
        {
            return huffmanEncoding(xArray,
                                   yArray,
                                   zArray,
                                   m_huffmanDict,
                                   m_xEncoded,
                                   m_yEncoded,
                                   m_zEncoded,
                                   m_pStop);
        }
        case ARITHMETIC_ENCODING :
        {
            return arithmeticEncoding(xArray,
                                      yArray,
                                      zArray,
                                      m_arithmeticDict,
                                      m_xEncoded,
                                      m_pStop);
        }
        case NO_ENCODING :
        {
            return false;
        }
    }
    return true;
}


/********************************************//**
\brief decode arrays of 0 and 1 bytes to original symbols
\encodingType : decoding type (0:huffman, 1:arithmetic)
\return boolean true if the execution is finished,
        false if the process has been stopped during the execution
***********************************************/
bool Fibers::decode(const EncodingType decodingType)
{
    // Check flag to know if we need to stop the process or not
    if(*m_pStop)
    {
        return false;
    }
    
    // Decode according to decodingType
    switch(decodingType)
    {
        case HUFFMAN_ENCODING :
        {
            if(!huffmanDecoding<float>(m_huffmanDict, m_xEncoded, m_yEncoded, m_zEncoded, m_xPointArray, m_yPointArray, m_zPointArray, m_countLines, m_origLineLength, m_pStop))
            {
                return false;
            }
            break;
        }
        case ARITHMETIC_ENCODING :
        {
            if(!arithmeticDecoding<float>(m_arithmeticDict, m_xEncoded, m_xPointArray, m_yPointArray, m_zPointArray, m_countLines, m_origLineLength, m_pStop))
            {
                return false;
            }
            break;
        }
        default :
        {
            return false;
        }
    }

    return true;
}


/********************************************//**
\brief compute compression ratio between input and output files
\return compression ratio
***********************************************/
float Fibers::getCompressionRatio()
{
    int origSizeBytes = getFileSizeBytes(m_inputPath);
    int compSizeBytes = getFileSizeBytes(m_outputPath);
    float cr = 100.0f * ((float)(origSizeBytes - compSizeBytes) / (float)origSizeBytes);
    saveResults(QString("Compression Ratio : ") + QString::number(cr) + QString("\n"));
    return cr;
}


void Fibers::addFiberToPointArray(const QVector<Vector> array)
{
    assert(array.size() > 0);
    
    // Update m_xPointArray, m_yPointArray and m_zPointArray
    int lineSize = array.size();
    QVector<float> xLine(lineSize), yLine(lineSize), zLine(lineSize);
    for(int i = 0; i < lineSize; i++)
    {
        xLine[i] = array[i].x;
        yLine[i] = array[i].y;
        zLine[i] = array[i].z;
        updateGlobalMinMax(xLine[i]);
        updateGlobalMinMax(yLine[i]);
        updateGlobalMinMax(zLine[i]);
    }
    m_xPointArray.push_back(xLine);
    m_yPointArray.push_back(yLine);
    m_zPointArray.push_back(zLine);
    
    assert(m_xPointArray.size() == m_yPointArray.size() && m_xPointArray.size() == m_zPointArray.size());
    
    m_countLines = m_xPointArray.size();
}


void Fibers::setFiberInPointArray(const QVector<float> array, const int fibIdx, const int dimIdx)
{
    assert(array.size() > 0);
    
    // Update m_pointsArray and synchronize m_xPointArray, m_yPointArray and m_zPointArray with it
    int lineSize = array.size();
    for(int i = 0; i < lineSize; i++)
    {
        updateGlobalMinMax(array[i]);
    }
    
    if(dimIdx == 0)
    {
        m_xPointArray[fibIdx] = array;
        m_countLines = m_xPointArray.size();
    }
    if(dimIdx == 1)
    {
        m_yPointArray[fibIdx] = array;
        m_countLines = m_yPointArray.size();
    }
    if(dimIdx == 2)
    {
        m_zPointArray[fibIdx] = array;
        m_countLines = m_zPointArray.size();
    }
}


void Fibers::setColorArray(const QVector<unsigned char> array, const int fibIdx, const int dimIdx)
{
    assert(array.size() > 0);
    
    // Update m_pointsArray and synchronize m_xPointArray, m_yPointArray and m_zPointArray with it
    int lineSize = array.size();
    for(int i = 0; i < lineSize; i++)
    {
        if(dimIdx == 0)
        {
            m_rColorArray[fibIdx] = array;
        }
        else if(dimIdx == 1)
        {
            m_gColorArray[fibIdx] = array;
        }
        else if(dimIdx == 2)
        {
            m_bColorArray[fibIdx] = array;
        }
    }
}




/********************************************//**
\brief Manipulation method to assign a new value in m_xPointArray, m_yPointArray or m_zPointArray. 
       This method must be used if you need to assign a new value.
\param value : new value to assign
\param fibIdx : index of the fiber concerned
\param ptIdx : index of the point in the concerned fiber
\param ptIdx : index of the dimension concerned (0:x, 1:y, 2:z)
***********************************************/
void Fibers::setPointArray(const float value,
                            const int fibIdx,
                            const int ptIdx,
                            const int dimIdx)
{
    assert(fibIdx >= 0 && fibIdx < getCountLines());
    assert(ptIdx >= 0 && ptIdx < getLineSize(fibIdx));
    assert(dimIdx >= 0 && dimIdx < 3);
    
    // Update m_pointsArray and m_xPointArray, m_yPointArray or m_zPointArray at the right position
    if(dimIdx == 0)
    {
        m_xPointArray[fibIdx][ptIdx] = value;
    }
    if(dimIdx == 1)
    {
        m_yPointArray[fibIdx][ptIdx] = value;
    }
    if(dimIdx == 2)
    {
        m_zPointArray[fibIdx][ptIdx] = value;
    }
    updateGlobalMinMax(value);
}


/********************************************//**
\brief Manipulation method to assign new values for a specific dimension in m_xPointArray, m_yPointArray or m_zPointArray. This method must ALWAYS be used if you need to assign new values for a specific dimension.
\param array : list of 1D arrays that contains new values to assign and that corresponds to a position
\param i : index of the dimension where to assign the new values (0:x, 1:y, 2:z)
***********************************************/
void Fibers::setPointsArrayDim(const QList< QVector<float> > array,
                               const int i)
{
    assert(array.size() > 0);
    assert(i >= 0 && i < 3);
    
    // Update m_xPointArray, m_yPointArray or m_zPointArray according to i
    if(i == 0)
    {
        assert(array.size() == m_yPointArray.size() && array.size() == m_zPointArray.size());
        m_xPointArray = array;
    }
    else if(i == 1)
    {
        assert(array.size() == m_xPointArray.size() && array.size() == m_zPointArray.size());
        m_yPointArray = array;
    }
    else if(i == 2)
    {
        assert(array.size() == m_xPointArray.size() && array.size() == m_yPointArray.size());
        m_zPointArray = array;
    }
    assert(m_xPointArray.size() == m_yPointArray.size() && m_xPointArray.size() == m_zPointArray.size());

    
    // Update min / max values
    for(int f = 0; f < array.size(); f++)
    {
        for(int j = 0; j < array[f].size(); j++)
        {
            updateGlobalMinMax(array[f][j]);
        }
    }
    m_countLines = m_xPointArray.size();
}


/********************************************//**
\brief To Get the number of points of a specific fiber
\param l : index of the fiber concerned
\return number of points of the fiber l
***********************************************/
int Fibers::getLineSize(const int l)
{
    assert(m_xPointArray.size() >= l);
    assert(m_xPointArray[l].size() == m_yPointArray[l].size() &&
           m_xPointArray[l].size() == m_zPointArray[l].size());
    return m_xPointArray[l].size();
}


/********************************************//**
\brief To Get the number of fibers
\return total number of fibers
***********************************************/
int Fibers::getCountLines()
{
    assert(m_xPointArray.size() == m_yPointArray.size() && m_xPointArray.size() == m_zPointArray.size());
    return m_xPointArray.size();
}


/********************************************//**
\brief Clear all arrays before doing anything
***********************************************/
void Fibers::clearArrays()
{
    m_rColorArray.clear();
    m_gColorArray.clear();
    m_bColorArray.clear();
    m_xEncoded.clear();
    m_yEncoded.clear();
    m_zEncoded.clear();
    m_xPointArray.clear();
    m_yPointArray.clear();
    m_zPointArray.clear();
    m_origLineLength.clear();
}


void Fibers::updateGlobalMinMax(float f)
{
    m_globalMin = std::min(m_globalMin, f);
    m_globalMax = std::max(m_globalMax, f);
}

/********************************************//**
\brief General method to load a set of fibers
\param filename : path of the input file to load
\return boolean true if everything went well, else otherwise
***********************************************/
bool Fibers::load(const QString& filename)
{
    bool res = false;
    QString ext = getExtension(filename);
    
    // Load fibers
    if(ext == QString("tck"))
    {
        displayMessage(QString("Loading..."));
        m_fiberType = FIBERTYPE_MRTRIX;
        res = loadTCK( filename );
    }
    if(ext == QString("trk"))
    {
        displayMessage(QString("Loading..."));
        m_fiberType = FIBERTYPE_TRK;
        res = loadTRK( filename );
    }
    if(m_fiberType == FIBERTYPE_NOT_SET)
    {
        displayMessage(QString("Error : unknown input file type !"));
    }
    return res;
}


/********************************************//**
\brief Method to load a DMRI file
\param filename : path of the input file to load
\return boolean true if everything went well, else otherwise
***********************************************/
bool Fibers::loadDmri(const QString& filename)
{
    m_origLineLength.clear();
    
    FILE* pFile;
    
    // Open input file
    pFile = fopen(filename.toUtf8().constData(), "r");
    
    if(pFile == NULL)
    {
        return false;
    }
    
    char *pS1 = new char[10];
    char *pS2 = new char[10];
    char *pS3 = new char[10];
    char *pS4 = new char[10];
    float f1, f2, f3, f4, f5;
    int res, countLines;
    
    // Read header
    res = fscanf(pFile, "%f %s", &f1, pS1);
    res = fscanf(pFile, "%f %s %s %s %s", &f1, pS1, pS2, pS3, pS4);
    res = fscanf(pFile, "%f", &f1);
    res = fscanf(pFile, "%f %f %f %f %f", &f1, &f2, &f3, &f4, &f5);
    res = fscanf(pFile, "%f %f %f %f %f", &f1, &f2, &f3, &f4, &f5);
    res = fscanf(pFile, "%f %f %f %f %f", &f1, &f2, &f3, &f4, &f5);
    res = fscanf(pFile, "%d %f", &countLines, &f2);
    
    delete[] pS1;
    delete[] pS2;
    delete[] pS3;
    delete[] pS4;
    
    pS1 = NULL;
    pS2 = NULL;
    pS3 = NULL;
    pS4 = NULL;
    
    // Read points
    int countPoints = 0;
    float back, front;
    
    for( int i = 0; i < countLines; i++ )
    {
        QVector<Vector> fiberArray;
        res = fscanf( pFile, "%f %f %f", &back, &front, &f1 );
        int nbpoints = back + front;
        
        if( back != 0 && front != 0 )
        {
            nbpoints--;
        }
        if( nbpoints > 0 )
        {
            // back points
            for( int j = back - 1; j >= 0; j-- )
            {
                res = fscanf( pFile, "%f %f %f %f", &f1, &f2, &f3, &f4 );
                Vector point(f1, f2, f3);
                fiberArray.push_back(point);
            }
            if( back != 0 && front != 0 )
            {
                //repeated pts
                res = fscanf( pFile, "%f %f %f %f", &f1, &f2, &f3, &f4 );
            }
            // front points
            for( int j = back; j < nbpoints; j++ )
            {
                res = fscanf( pFile, "%f %f %f %f", &f1, &f2, &f3, &f4 );
                Vector point(f1, f2, f3);
                fiberArray.push_back(point);
            }
            countPoints += fiberArray.size();
        }
        m_origLineLength.push_back(fiberArray.size());
        addFiberToPointArray(fiberArray);
    }
    saveResults("Original total number of points : " + QString::number(countPoints) + "\n");
    
    fclose( pFile );
    return true;
}

#include <arpa/inet.h>

/********************************************//**
\brief Method to load a TCK file
\param filename : path of the input file to load
\return boolean true if everything went well, else otherwise
***********************************************/
bool Fibers::loadTCK(const QString& filename)
{
    long int pc = 0;
    size_t data_offset(0);
    
    converterByteFloat cbf;
    m_origLineLength.clear();
    
    //Open file
    FILE *pFs = fopen(filename.toStdString().c_str(), "r");
    
    // Read header
    char lineBuffer[200];
    std::string readLine("");

    // TODO check datatype
    while(readLine.find( "END" ) == std::string::npos)
    {
        if(*m_pStop)
        {
            return false;
        }
        
        fgets(lineBuffer, 200, pFs);
        readLine = std::string(lineBuffer);
        
        // Get data offset.
        if(readLine.find("file") != std::string::npos)
        {
            sscanf(lineBuffer, "file: . %ld", &pc);
            data_offset = pc;
        }
    }
    fclose(pFs);

    
    // Get size of the file.
    FILE *pDataFile = fopen(filename.toStdString().c_str(), "rb");
    fseek(pDataFile, 0, SEEK_END);
    size_t nSize = ftell(pDataFile);
    
    fseek(pDataFile, data_offset, SEEK_SET);
    
    int countPoints = 0; // number of points
    size_t remainingSize(nSize - data_offset);
    
    
    float x, y, z;
    unsigned char *pBuffer = new unsigned char[12];
    //QList< QVector<Vector> > pointsArray;
    QVector<Vector> fiberArray;
    
    while(remainingSize > 0)
    {
        fiberArray.resize(0);
        
        fread((void*) pBuffer, 1, 12, pDataFile);
        
        // read first point of the track
        pc = 0;
        cbf.b[0] = pBuffer[pc++];
        cbf.b[1] = pBuffer[pc++];
        cbf.b[2] = pBuffer[pc++];
        cbf.b[3] = pBuffer[pc++];
        x = cbf.f;
        cbf.b[0] = pBuffer[pc++];
        cbf.b[1] = pBuffer[pc++];
        cbf.b[2] = pBuffer[pc++];
        cbf.b[3] = pBuffer[pc++];
        y = cbf.f;
        cbf.b[0] = pBuffer[pc++];
        cbf.b[1] = pBuffer[pc++];
        cbf.b[2] = pBuffer[pc++];
        cbf.b[3] = pBuffer[pc++];
        z = cbf.f;
        
        remainingSize -= 12;
        
        // If points is Inf, it is the end of the file
        if(x == std::numeric_limits<float>::infinity() &&
           y == std::numeric_limits<float>::infinity() &&
           z == std::numeric_limits<float>::infinity())
        {
            break;
        }
        
        Vector point(x, y, z);
        fiberArray.push_back(point);
        
        //Read points (x,y,z) until x2 equals NaN (0x0000C07F), meaning end of the tract.
        while(remainingSize > 0 && !(cbf.b[0] == 0x00 && cbf.b[1] == 0x00 && cbf.b[2] == 0xC0 && cbf.b[3] == 0x7F))
        {
            if(*m_pStop)
            {
                return false;
            }
            
            fread((void*) pBuffer, 1, 12, pDataFile);
            
            if(*m_pStop)
            {
                return false;
            }
            
            // read next point
            pc = 0;
            cbf.b[0] = pBuffer[pc++];
            cbf.b[1] = pBuffer[pc++];
            cbf.b[2] = pBuffer[pc++];
            cbf.b[3] = pBuffer[pc++];
            x = cbf.f;
            cbf.b[0] = pBuffer[pc++];
            cbf.b[1] = pBuffer[pc++];
            cbf.b[2] = pBuffer[pc++];
            cbf.b[3] = pBuffer[pc++];
            y = cbf.f;
            cbf.b[0] = pBuffer[pc++];
            cbf.b[1] = pBuffer[pc++];
            cbf.b[2] = pBuffer[pc++];
            cbf.b[3] = pBuffer[pc++];
            z = cbf.f;
            
            remainingSize -= 12;
            
            // Add point if not NaN
            if(!(cbf.b[0] == 0x00 && cbf.b[1] == 0x00 && cbf.b[2] == 0xC0 && cbf.b[3] == 0x7F))
            {
                Vector point(x, y, z);
                fiberArray.push_back(point);
            }
            cbf.f = x;
        }
        addFiberToPointArray(fiberArray);
        m_origLineLength.push_back(fiberArray.size());
        countPoints += fiberArray.size();
    }
    fclose(pDataFile);
    delete[] pBuffer;
    pBuffer = NULL;
    
    saveResults("Original total number of points : " + QString::number(countPoints) + "\n");
    
    return true;
}


/********************************************//**
\brief Method to load a TRK file
\param filename : path of the input file to load
\return boolean true if everything went well, else otherwise
***********************************************/
bool Fibers::loadTRK(const QString& filename)
{
    converterByteFloat cbf;
    
    // READ HEADER
    m_pTrkHdr = new TRK_hdr();
    FILE* fp = fopen(filename.toStdString().c_str(), "rb");
    if(fp != NULL)
    {
        int n = fread(m_pTrkHdr, sizeof(*m_pTrkHdr), 1, fp);
        if(n != 1)
        {
            fclose(fp);
            return false;
        }
    }
    
    int countPoints = 0;
    QVector<Vector> fiberArray;
    m_rColorArray.clear();
    m_gColorArray.clear();
    m_bColorArray.clear();
    
    QVector<unsigned char> rArray, gArray, bArray;
    QVector< QVector<float> > pArray;
    
    for( int i = 0; i < m_pTrkHdr->nbCount; ++i )
    {
        fiberArray.clear();
        rArray.clear();
        gArray.clear();
        bArray.clear();
        pArray.clear();
        
        //Number of points in this track. [4 bytes]
        int nbPoints;
        int n = fread(&nbPoints, sizeof(int), 1, fp);
        if(n != 1)
        {
            fclose(fp);
            return false;
        }
        
        //Read data of one track.
        size_t ptsSize = 3 + m_pTrkHdr->nbScalars;
        size_t tractSize = 4 * (nbPoints * ptsSize + m_pTrkHdr->nbProperties);
        char* pBuffer = new char[tractSize];
        
        // Read tract
        n = fread(pBuffer, 1, tractSize, fp);
        if(n != tractSize)
        {
            fclose(fp);
            return false;
        }
        
        for( unsigned int j = 0; j != nbPoints; ++j )
        {
            QVector<float> color;
            QVector<float> point;
            QVector<float> prop;
            
            //Read coordinates (x,y,z) and scalars associated to each point.
            for( unsigned int k = 0; k < ptsSize; ++k )
            {
                memcpy(cbf.b, &pBuffer[4 * ( j * ptsSize + k )], 4);
                
                if( k > 6) // Properties at each point
                {
                    prop.push_back(cbf.f);
                }
                if( k >= 3 ) //RGB color of each point.
                {
                    m_colorArrayLoadedFromFile = true;
                    color.push_back(cbf.f);
                }
                else // Point
                {
                    point.push_back(cbf.f);
                }
            }
            if(color.size() == 3)
            {
                rArray.push_back(color[0]);
                gArray.push_back(color[1]);
                bArray.push_back(color[2]);
            }
            if(point.size() == 3)
            {
                fiberArray.push_back(Vector(point[0], point[1], point[2]));
            }
            if(prop.size() > 0)
            {
                pArray.push_back(prop);
            }
        }
        countPoints += fiberArray.size();
        addFiberToPointArray(fiberArray);
        m_rColorArray.push_back(rArray);
        m_gColorArray.push_back(gArray);
        m_gColorArray.push_back(bArray);
        m_propertiesArray.push_back(pArray);
        
        delete[] pBuffer;
        pBuffer = NULL;
    }
    
    fclose(fp);
    
    saveResults("Original total number of points : " + QString::number(countPoints) + "\n");

    return true;
}


/********************************************//**
\brief Method to load a compressed file
\param filename : path of the input file to load
\param encodingType : decoding type that will be read during the loading process
\param ttype : transformation type that will be read during the loading process
\return boolean true if everything went well, else otherwise
***********************************************/
bool Fibers::loadCompressed(const QString& filename,
                            const EncodingType& decodingType,
                            const TransformationType& ttype)
{
    clearArrays();
    // Open file
    std::ifstream myfile(filename.toStdString().c_str(), std::ios::in | std::ios::binary);
    
    if(!myfile.is_open())
    {
        std::cout << "Can't open input file !" << std::endl;
        return false;
    }
    
    // Read fibertype, number of points and number of lines
    int fiberType, totalLines;
    myfile.read((char*)&(fiberType), sizeof(fiberType));
    myfile.read((char*)&totalLines, sizeof(totalLines));
    m_fiberType = (FiberType)fiberType;
    m_countLines = totalLines;
    
    myfile.read((char*)&(decodingType), sizeof(decodingType));
    myfile.read((char*)&(ttype), sizeof(ttype));
    
    // Read encoded signal or points
    switch(decodingType)
    {
        case HUFFMAN_ENCODING :
        {
            // Read size of lines
            for( int i = 0; i < m_countLines; i++ )
            {
                int line;
                myfile.read((char*)&line, sizeof(line));
                m_origLineLength.push_back(line);
            }
            
            if(!readHuffmanEncodedSignal<float>(myfile, m_xEncoded, m_yEncoded, m_zEncoded, m_huffmanDict, m_pStop))
            {
                return false;
            }
            break;
        }
        case ARITHMETIC_ENCODING :
        {
            // Read size of lines
            for( int i = 0; i < m_countLines; i++ )
            {
                int line;
                myfile.read((char*)&line, sizeof(line));
                m_origLineLength.push_back(line);
            }
            
            if(!readArithmeticEncodedSignal<float>(myfile, m_xEncoded, m_arithmeticDict, m_pStop))
            {
                return false;
            }
            break;
        }
        case NO_ENCODING :
        {
            // Read points
            myfile.close();
            converterByteFloat cbf;
            size_t data_offset(0);
            
            FILE *pDataFile = fopen(filename.toStdString().c_str(), "rb");
            fseek(pDataFile, 0, SEEK_END);
            size_t nSize = ftell(pDataFile);
            fseek(pDataFile, data_offset, SEEK_SET);
            
            unsigned char *pBuffer = new unsigned char[nSize - data_offset];
            fread((void*) pBuffer, sizeof(unsigned char), nSize - data_offset, pDataFile);
            fclose(pDataFile);
            pDataFile = NULL;
            
            long int pc = 0;
            int countPoints = 0; // number of points
            size_t remainingSize(nSize - data_offset);
            
            while(remainingSize > 0)
            {
                QVector<Vector> fiberArray;
                
                // read one tract
                cbf.b[0] = pBuffer[pc++];
                cbf.b[1] = pBuffer[pc++];
                cbf.b[2] = pBuffer[pc++];
                cbf.b[3] = pBuffer[pc++];
                float x = cbf.f;
                cbf.b[0] = pBuffer[pc++];
                cbf.b[1] = pBuffer[pc++];
                cbf.b[2] = pBuffer[pc++];
                cbf.b[3] = pBuffer[pc++];
                float y = cbf.f;
                cbf.b[0] = pBuffer[pc++];
                cbf.b[1] = pBuffer[pc++];
                cbf.b[2] = pBuffer[pc++];
                cbf.b[3] = pBuffer[pc++];
                float z = cbf.f;
                
                Vector point(x, y, z);
                fiberArray.push_back(point);
                
                float x2 = x;
                cbf.f = x2;
                remainingSize -= 12;
                
                //Read points (x,y,z) until x2 equals NaN (0x0000C07F), meaning end of the tract.
                while(remainingSize > 0 && !(cbf.b[0] == 0x00 && cbf.b[1] == 0x00 && cbf.b[2] == 0xC0 && cbf.b[3] == 0x7F))
                {
                    cbf.b[0] = pBuffer[pc++];   // get next float
                    cbf.b[1] = pBuffer[pc++];
                    cbf.b[2] = pBuffer[pc++];
                    cbf.b[3] = pBuffer[pc++];
                    float x2 = cbf.f;
                    remainingSize -= 4;
                    
                    cbf.b[0] = pBuffer[pc++];
                    cbf.b[1] = pBuffer[pc++];
                    cbf.b[2] = pBuffer[pc++];
                    cbf.b[3] = pBuffer[pc++];
                    float y2 = cbf.f;
                    remainingSize -= 4;
                    
                    cbf.b[0] = pBuffer[pc++];
                    cbf.b[1] = pBuffer[pc++];
                    cbf.b[2] = pBuffer[pc++];
                    cbf.b[3] = pBuffer[pc++];
                    float z2 = cbf.f;
                    remainingSize -= 4;
                    
                    x = x2;
                    y = y2;
                    z = z2;
                    
                    Vector point(x, y, z);
                    fiberArray.push_back(point);
                    cbf.f = x2;
                }
                
                // Remove last point, which is the NaN
                fiberArray.pop_back();
                addFiberToPointArray(fiberArray);
                m_origLineLength.push_back(fiberArray.size());
                countPoints += fiberArray.size();
            }
            delete[] pBuffer;
            pBuffer = NULL;
            break;
        }
    }
    
    // Read colors
    bool readColors = false;
    myfile.read((char*)&readColors, sizeof(readColors));
    if(readColors)
    {
        // Decode according to decodingType
        switch(decodingType)
        {
            case HUFFMAN_ENCODING :
            {
                // Read encoded signal
                if( !readHuffmanEncodedSignal<unsigned char>(myfile, m_rEncodedColor, m_gEncodedColor,
                                                           m_bEncodedColor, m_colorHuffmanDict, m_pStop))
                {
                    return false;
                }
                // Huffman decoding in m_colorArray
                if(!huffmanDecoding<unsigned char>(m_colorHuffmanDict, m_rEncodedColor, m_gEncodedColor, m_bEncodedColor, m_rColorArray, m_gColorArray, m_bColorArray, m_countLines, m_origLineLength, m_pStop))
                {
                    return false;
                }
                
                break;
            }
            case ARITHMETIC_ENCODING :
            {
                // Read encoded signal
                if( !readArithmeticEncodedSignal<unsigned char>(myfile, m_rEncodedColor, m_colorArithmeticDict, m_pStop) )
                {
                    return false;
                }
                // Arithmetic decoding in m_colorArray
                if(!arithmeticDecoding<unsigned char>(m_colorArithmeticDict, m_rEncodedColor, m_rColorArray, m_gColorArray, m_bColorArray, m_countLines, m_origLineLength, m_pStop))
                {
                    return false;
                }
                
                break;
            }
            case NO_ENCODING :
            {
                // Read points
                myfile.close();
                converterByteFloat cbf;
                size_t data_offset(0);
                
                FILE *pDataFile = fopen(filename.toStdString().c_str(), "rb");
                fseek(pDataFile, 0, SEEK_END);
                size_t nSize = ftell(pDataFile);
                fseek(pDataFile, data_offset, SEEK_SET);
                
                unsigned char *pBuffer = new unsigned char[nSize - data_offset];
                fread((void*) pBuffer, sizeof(unsigned char), nSize - data_offset, pDataFile);
                fclose(pDataFile);
                pDataFile = NULL;
                
                long int pc = 0;
                int countPoints = 0; // number of points
                size_t remainingSize(nSize - data_offset);
                
                while(remainingSize > 0)
                {
                    QVector<Vector> colorArray;
                    
                    // read one tract
                    cbf.b[0] = pBuffer[pc++];
                    cbf.b[1] = pBuffer[pc++];
                    cbf.b[2] = pBuffer[pc++];
                    cbf.b[3] = pBuffer[pc++];
                    float x = cbf.f;
                    cbf.b[0] = pBuffer[pc++];
                    cbf.b[1] = pBuffer[pc++];
                    cbf.b[2] = pBuffer[pc++];
                    cbf.b[3] = pBuffer[pc++];
                    float y = cbf.f;
                    cbf.b[0] = pBuffer[pc++];
                    cbf.b[1] = pBuffer[pc++];
                    cbf.b[2] = pBuffer[pc++];
                    cbf.b[3] = pBuffer[pc++];
                    float z = cbf.f;
                    
                    Vector point(x, y, z);
                    colorArray.push_back(point);
                    
                    float x2 = x;
                    cbf.f = x2;
                    remainingSize -= 12;
                    
                    //Read points (x,y,z) until x2 equals NaN (0x0000C07F), meaning end of the tract.
                    while(remainingSize > 0 && !(cbf.b[0] == 0x00 && cbf.b[1] == 0x00 && cbf.b[2] == 0xC0 && cbf.b[3] == 0x7F))
                    {
                        cbf.b[0] = pBuffer[pc++];   // get next float
                        cbf.b[1] = pBuffer[pc++];
                        cbf.b[2] = pBuffer[pc++];
                        cbf.b[3] = pBuffer[pc++];
                        float x2 = cbf.f;
                        remainingSize -= 4;
                        
                        cbf.b[0] = pBuffer[pc++];
                        cbf.b[1] = pBuffer[pc++];
                        cbf.b[2] = pBuffer[pc++];
                        cbf.b[3] = pBuffer[pc++];
                        float y2 = cbf.f;
                        remainingSize -= 4;
                        
                        cbf.b[0] = pBuffer[pc++];
                        cbf.b[1] = pBuffer[pc++];
                        cbf.b[2] = pBuffer[pc++];
                        cbf.b[3] = pBuffer[pc++];
                        float z2 = cbf.f;
                        remainingSize -= 4;
                        
                        x = x2;
                        y = y2;
                        z = z2;
                        
                        Vector point(x, y, z);
                        colorArray.push_back(point);
                        cbf.f = x2;
                    }
                    
                    // Remove last point, which is the NaN
                    colorArray.pop_back();
                    addFiberToPointArray(colorArray);
                    m_origLineLength.push_back(colorArray.size());
                    countPoints += colorArray.size();
                }
                delete[] pBuffer;
                pBuffer = NULL;
                break;
            }
        }
        
    }
    m_colorArrayLoadedFromFile = readColors;
    
    myfile.close();
    return true;
}


/********************************************//**
\brief Utility method to save output statistics in text file
\param str : string to save in the text file
\return boolean true if everything went well, else otherwise
***********************************************/
bool Fibers::saveResults(const QString& str)
{
    #if _SAVE_STATS
        if(!m_outputStatsFilename.isEmpty())
        {
            // Open file
            m_resultFile.setFileName(m_outputStatsFilename);
            if(!m_resultFile.open(QIODevice::WriteOnly | QIODevice::Append))
            {
                return false;
            }
            // Add text
            QTextStream out(&m_resultFile);
            out << str;
            
            // Close file
            m_resultFile.close();
        }
    #endif
    return true;
}


/********************************************//**
\brief Method to save a compressed file
\param encodingType : encoding type that will be saved in the compressed file
\param ttype : transformation type that will be saved in the compressed file
\param filename : path of the output file
\return boolean true if everything went well, else otherwise
***********************************************/
bool Fibers::saveCompressed(const EncodingType encodingType,
                            const TransformationType ttype,
                            const QString& filename)
{
    QVector<int> linesToSave;
    QString name = filename;
    QString ext = getExtension(name);
	if(ext != QString("zfib"))
    {
        name += QString(".zfib");
    }
    
    // Open file
    std::ofstream myfile(name.toStdString().c_str(), std::ios::out | std::ios::binary);
    if(!myfile.is_open())
    {
        std::cout << "Can't open input file !" << std::endl;
        return false;
    }
    
    // Write fibertype, total number of points and total number of lines
    int countLines = getCountLines();
    myfile.write((const char*) &(m_fiberType), sizeof(m_fiberType));
    myfile.write((const char*) &(countLines), sizeof(countLines));
    
    // Write encoded signal or points
    myfile.write((const char*) &(encodingType), sizeof(encodingType));
    myfile.write((const char*) &(ttype), sizeof(ttype));
    
    // Check flag to know if we need to stop the process or not
    if(*m_pStop)
    {
        return false;
    }
    
    switch(encodingType)
    {
            // Write encoded points
        case HUFFMAN_ENCODING :
        {
            // Write size of each line
            for(int i = 0; i < countLines; ++i)
            {
                int lineSize = getLineSize(i);
                myfile.write((const char*) &(lineSize), sizeof(lineSize));
            }
            writeHuffmanEncodedSignal<float>(myfile, m_xEncoded, m_yEncoded,
                                             m_zEncoded, m_huffmanDict, m_pStop);
            break;
        }
        case ARITHMETIC_ENCODING :
        {
            // Write size of each line
            for(int i = 0; i < countLines; ++i)
            {
                int lineSize = getLineSize(i);
                myfile.write((const char*) &(lineSize), sizeof(lineSize));
            }
            writeArithmeticEncodedSignal<float>(myfile, m_xEncoded, m_arithmeticDict, m_pStop);
            break;
        }
        case NO_ENCODING :
        {
            // Write points directly
            myfile.close();
            FILE *pFs = fopen(name.toStdString().c_str(), "ab") ;
            
            converterByteFloat cbf;
            for(int f = 0; f < getCountLines(); ++f)
            {
                // Check flag to know if we need to stop the process or not
                if(*m_pStop)
                {
                    return false;
                }
                for( int j = 0; j < getLineSize(f); ++j )
                {
                    cbf.f = m_xPointArray[f][j];
                    fputc(cbf.b[0], pFs);
                    fputc(cbf.b[1], pFs);
                    fputc(cbf.b[2], pFs);
                    fputc(cbf.b[3], pFs);
                    
                    cbf.f = m_xPointArray[f][j];
                    fputc(cbf.b[0], pFs);
                    fputc(cbf.b[1], pFs);
                    fputc(cbf.b[2], pFs);
                    fputc(cbf.b[3], pFs);
                    
                    cbf.f = m_xPointArray[f][j];
                    fputc(cbf.b[0], pFs);
                    fputc(cbf.b[1], pFs);
                    fputc(cbf.b[2], pFs);
                    fputc(cbf.b[3], pFs);
                }
                // End a tract with 3 NaNs.
                cbf.b[0] = 0x00;
                cbf.b[1] = 0x00;
                cbf.b[2] = 0xC0;
                cbf.b[3] = 0x7F;
                fputc(cbf.b[0], pFs);
                fputc(cbf.b[1], pFs);
                fputc(cbf.b[2], pFs);
                fputc(cbf.b[3], pFs);
                fputc(cbf.b[0], pFs);
                fputc(cbf.b[1], pFs);
                fputc(cbf.b[2], pFs);
                fputc(cbf.b[3], pFs);
                fputc(cbf.b[0], pFs);
                fputc(cbf.b[1], pFs);
                fputc(cbf.b[2], pFs);
                fputc(cbf.b[3], pFs);
            }
            break;
        }
    }
    
    // Write encoded colors if needed
    myfile.write((char*)&m_colorArrayLoadedFromFile, sizeof(m_colorArrayLoadedFromFile));
    
    if(m_colorArrayLoadedFromFile)
    {
        // Check flag to know if we need to stop the process or not
        if(*m_pStop)
        {
            return false;
        }
        
        std::string rEncodedColor, gEncodedColor, bEncodedColor;
        switch(encodingType)
        {
                // Encode colors
            case HUFFMAN_ENCODING :
            {
                if(!huffmanEncoding(m_rColorArray, m_gColorArray,m_bColorArray,
                                    m_colorHuffmanDict, rEncodedColor, gEncodedColor,
                                    bEncodedColor, m_pStop))
                {
                    return false;
                }
                writeHuffmanEncodedSignal<signed char>(myfile,m_rEncodedColor, m_gEncodedColor,m_bEncodedColor, m_colorHuffmanDict, m_pStop);
                break;
            }
            case ARITHMETIC_ENCODING :
            {
                if(!arithmeticEncoding(m_rColorArray, m_gColorArray, m_bColorArray, m_colorArithmeticDict, m_rEncodedColor, m_pStop))
                {
                    return false;
                }
                writeArithmeticEncodedSignal<unsigned char>(myfile, m_rEncodedColor, m_colorArithmeticDict, m_pStop);
                break;
            }
            case NO_ENCODING :
            {
                // Write colors directly
                myfile.close();
                FILE *pFs = fopen(name.toStdString().c_str(), "ab") ;
                
                converterByteFloat cbf;
                for(int f = 0; f < m_rColorArray.size(); ++f)
                {
                    // Check flag to know if we need to stop the process or not
                    if(*m_pStop)
                    {
                        return false;
                    }
                    for( int j = 0; j < m_rColorArray[f].size(); ++j )
                    {
                        cbf.f = m_rColorArray[f][j];
                        fputc(cbf.b[0], pFs);
                        fputc(cbf.b[1], pFs);
                        fputc(cbf.b[2], pFs);
                        fputc(cbf.b[3], pFs);
                        
                        cbf.f = m_gColorArray[f][j];
                        fputc(cbf.b[0], pFs);
                        fputc(cbf.b[1], pFs);
                        fputc(cbf.b[2], pFs);
                        fputc(cbf.b[3], pFs);
                        
                        cbf.f = m_bColorArray[f][j];
                        fputc(cbf.b[0], pFs);
                        fputc(cbf.b[1], pFs);
                        fputc(cbf.b[2], pFs);
                        fputc(cbf.b[3], pFs);
                    }
                    // End a tract with 3 NaNs.
                    cbf.b[0] = 0x00;
                    cbf.b[1] = 0x00;
                    cbf.b[2] = 0xC0;
                    cbf.b[3] = 0x7F;
                    fputc(cbf.b[0], pFs);
                    fputc(cbf.b[1], pFs);
                    fputc(cbf.b[2], pFs);
                    fputc(cbf.b[3], pFs);
                    fputc(cbf.b[0], pFs);
                    fputc(cbf.b[1], pFs);
                    fputc(cbf.b[2], pFs);
                    fputc(cbf.b[3], pFs);
                    fputc(cbf.b[0], pFs);
                    fputc(cbf.b[1], pFs);
                    fputc(cbf.b[2], pFs);
                    fputc(cbf.b[3], pFs);
                }
                break;
            }
        }
    }
    
    // Close file
	myfile.close();
    
    return true;
}


/********************************************//**
\brief General method to save a file
\param filename : path of the output file
***********************************************/
bool Fibers::save(const QString& filename)
{
    // Save fibers file
    switch(m_fiberType)
    {
        case FIBERTYPE_MRTRIX :
            return saveTCK(filename);
        case FIBERTYPE_TRK :
            return saveTRK(filename);
        default :
            showMessage(QString("Error : unknown output file type !"));
    }
    return true;
}


/********************************************//**
\brief Method to save a DMRI file
\param filename : path of the output file
***********************************************/
bool Fibers::saveDMRI(const QString& filename)
{
    std::ofstream myfile;
    QString ext = getExtension(filename);
    
	if(ext != QString("fib"))
    {
        return false;
    }
    
    // Open file
    myfile.open(filename.toStdString().c_str(), std::ios::out);
    
    // Write header
	myfile << "1 FA\n4 min max mean var\n1\n4 0 0 0 0\n4 0 0 0 0\n4 0 0 0 0\n";
	myfile << getCountLines() << " " << 0.5 << "\n";
    
    // Write points
    for(int f = 0; f < getCountLines(); ++f)
    {
        // Check flag to know if we need to stop the process or not
        if(*m_pStop)
        {
            myfile.close();
            std::remove(filename.toStdString().c_str());
            return false;
        }
        myfile << getLineSize(f) << " 1\n1\n";
        
        for( int j = 0; j < getLineSize(f); ++j )
        {
            myfile << m_xPointArray[f][j] << " " << m_xPointArray[f][j] << " " << m_zPointArray[f][j] << " 0\n";
        }
        myfile << m_xPointArray[f][0] << " " << m_yPointArray[f][0] << " " << m_zPointArray[f][0] << " 0\n";
    }
    
    // Close file
    myfile.close();
    return true;
}


/********************************************//**
\brief Method to save a TCK file
\param filename : path of the output file
***********************************************/
bool Fibers::saveTCK(const QString& filename)
{
    std::ofstream myfile;
    
    QString ext = getExtension(filename);
    
	if(ext != QString("tck"))
    {
        return false;
    }
    
    // Open file
    myfile.open(filename.toStdString().c_str(), std::ios::out);
    
    // Create header in ostringstream to be able to get its length.
    std::ostringstream header;
    header << "mrtrix tracks\n";
    header << "count: " << m_countLines << "\n";
    header << "datatype: Float32LE\n";
    header << "file: . ";
    
    // To compute the offset of the data in the file, which is used in the
    // "file" field of the header, we need to get the length of the header.
    // +5 is for "\nEND\n"
    size_t header_length = header.tellp();
    header_length += 5;
    
    // We also need to get the length of the header length string,
    // to be able to add it to the offset.
    std::ostringstream temp;
    temp << header_length;
    size_t header_length_str_length = temp.tellp();
    size_t offset = header_length + header_length_str_length;
    
    // We also need to make sure that the final offset takes into account the
    // fact that the addition can raise the length of the offset by one
    // decimal position.
    std::ostringstream temp1;
    temp1 << offset;
    size_t offset_str_length = temp1.tellp();
    
    if(offset_str_length != header_length_str_length)
    {
        offset += 1;
    }
    
    header << offset << "\n";
    header << "END\n";
    
    myfile << header.str();
    myfile.close();
    
    
    // Reopen the file in binary, to write the bytes.
    FILE *pFs = fopen(filename.toStdString().c_str(), "ab") ;
    
    converterByteFloat cbf;
    for(int f = 0; f < getCountLines(); ++f)
    {
        // Check flag to know if we need to stop the process or not
        if(*m_pStop)
        {
            fclose(pFs);
            std::remove(filename.toStdString().c_str());
            return false;
        }
        for( int j = 0; j < getLineSize(f); ++j )
        {
            cbf.f = m_xPointArray[f][j];
            fputc(cbf.b[0], pFs);
            fputc(cbf.b[1], pFs);
            fputc(cbf.b[2], pFs);
            fputc(cbf.b[3], pFs);
            
            cbf.f = m_yPointArray[f][j];
            fputc(cbf.b[0], pFs);
            fputc(cbf.b[1], pFs);
            fputc(cbf.b[2], pFs);
            fputc(cbf.b[3], pFs);
            
            cbf.f = m_zPointArray[f][j];
            fputc(cbf.b[0], pFs);
            fputc(cbf.b[1], pFs);
            fputc(cbf.b[2], pFs);
            fputc(cbf.b[3], pFs);
        }
        // End a tract with 3 NaNs.
        cbf.b[0] = 0x00;
        cbf.b[1] = 0x00;
        cbf.b[2] = 0xC0;
        cbf.b[3] = 0x7F;
        fputc(cbf.b[0], pFs);
        fputc(cbf.b[1], pFs);
        fputc(cbf.b[2], pFs);
        fputc(cbf.b[3], pFs);
        fputc(cbf.b[0], pFs);
        fputc(cbf.b[1], pFs);
        fputc(cbf.b[2], pFs);
        fputc(cbf.b[3], pFs);
        fputc(cbf.b[0], pFs);
        fputc(cbf.b[1], pFs);
        fputc(cbf.b[2], pFs);
        fputc(cbf.b[3], pFs);
    }
    // End the file with 3 Infs.
    cbf.b[0] = 0x00;
    cbf.b[1] = 0x00;
    cbf.b[2] = 0x80;
    cbf.b[3] = 0x7F;
    fputc(cbf.b[0], pFs);
    fputc(cbf.b[1], pFs);
    fputc(cbf.b[2], pFs);
    fputc(cbf.b[3], pFs);
    fputc(cbf.b[0], pFs);
    fputc(cbf.b[1], pFs);
    fputc(cbf.b[2], pFs);
    fputc(cbf.b[3], pFs);
    fputc(cbf.b[0], pFs);
    fputc(cbf.b[1], pFs);
    fputc(cbf.b[2], pFs);
    fputc(cbf.b[3], pFs);
    
    fclose(pFs);
    return true;
}



/********************************************//**
\brief Method to save a TRK file
\param filename : path of the output file
***********************************************/
bool Fibers::saveTRK(const QString& filename)
{
    QString ext = getExtension(filename);
    
	if(ext != QString("trk"))
    {
        return false;
    }
    
    // Open file
    FILE *pFs = fopen(filename.toStdString().c_str(), "wb+") ;
    if(pFs == NULL)
    {
        return false;
    }
    
    // Write header
    fwrite((const char*)m_pTrkHdr->id_string, sizeof(m_pTrkHdr->id_string), 1, pFs);
    fwrite((const char*)m_pTrkHdr->dim, sizeof(m_pTrkHdr->dim), 1, pFs);
    fwrite((const char*)m_pTrkHdr->voxelSize, sizeof(m_pTrkHdr->voxelSize), 1, pFs);
    fwrite((const char*)m_pTrkHdr->origin, sizeof(m_pTrkHdr->origin), 1, pFs);
    fwrite((const char*)&m_pTrkHdr->nbScalars, sizeof(m_pTrkHdr->nbScalars), 1, pFs);
    fwrite((const char*)m_pTrkHdr->scalarNames, sizeof(m_pTrkHdr->scalarNames), 1, pFs);
    fwrite((const char*)&m_pTrkHdr->nbProperties, sizeof(m_pTrkHdr->nbProperties), 1, pFs);
    fwrite((const char*)m_pTrkHdr->propertyNames, sizeof(m_pTrkHdr->propertyNames), 1, pFs);
    fwrite((const char*)m_pTrkHdr->voxToRas, sizeof(m_pTrkHdr->voxToRas), 1, pFs);
    fwrite((const char*)m_pTrkHdr->reserved, sizeof(m_pTrkHdr->reserved), 1, pFs);
    fwrite((const char*)m_pTrkHdr->voxelOrder, sizeof(m_pTrkHdr->voxelOrder), 1, pFs);
    fwrite((const char*)m_pTrkHdr->pad2, sizeof(m_pTrkHdr->pad2), 1, pFs);
    fwrite((const char*)m_pTrkHdr->imageOrientationPatient, sizeof(m_pTrkHdr->imageOrientationPatient), 1, pFs);
    fwrite((const char*)m_pTrkHdr->pad1, sizeof(m_pTrkHdr->pad1), 1, pFs);
    fwrite((const char*)&m_pTrkHdr->invertX, sizeof(m_pTrkHdr->invertX), 1, pFs);
    fwrite((const char*)&m_pTrkHdr->invertY, sizeof(m_pTrkHdr->invertY), 1, pFs);
    fwrite((const char*)&m_pTrkHdr->invertZ, sizeof(m_pTrkHdr->invertZ), 1, pFs);
    fwrite((const char*)&m_pTrkHdr->swapXY, sizeof(m_pTrkHdr->swapXY), 1, pFs);
    fwrite((const char*)&m_pTrkHdr->swapYZ, sizeof(m_pTrkHdr->swapYZ), 1, pFs);
    fwrite((const char*)&m_pTrkHdr->swapZX, sizeof(m_pTrkHdr->swapZX), 1, pFs);
    fwrite((const char*)&m_pTrkHdr->nbCount, sizeof(m_pTrkHdr->nbCount), 1, pFs);
    fwrite((const char*)&m_pTrkHdr->version, sizeof(m_pTrkHdr->version), 1, pFs);
    fwrite((const char*)&m_pTrkHdr->hdrSize, sizeof(m_pTrkHdr->hdrSize), 1, pFs);
    
    
    //converterByteFloat cbf;
    for(int f = 0; f < getCountLines(); ++f)
    {
        // Check flag to know if we need to stop the process or not
        if(*m_pStop)
        {
            fclose(pFs);
            std::remove(filename.toStdString().c_str());
            return false;
        }
        
        // Write number of points of this track
        int nbPoints = getLineSize(f);
        fwrite(&nbPoints, sizeof(nbPoints), 1, pFs);
        
        for( int j = 0; j < nbPoints; ++j )
        {
            // Write point
            fwrite(&m_xPointArray[f][j], sizeof(m_xPointArray[f][j]), 1, pFs);
            fwrite(&m_yPointArray[f][j], sizeof(m_yPointArray[f][j]), 1, pFs);
            fwrite(&m_zPointArray[f][j], sizeof(m_zPointArray[f][j]), 1, pFs);
                    
            // Write RGB color
            if(m_pTrkHdr->nbScalars >= 3)
            {
                fwrite(&m_rColorArray[f][j], sizeof(m_rColorArray[f][j]), 1, pFs);
                fwrite(&m_gColorArray[f][j], sizeof(m_gColorArray[f][j]), 1, pFs);
                fwrite(&m_bColorArray[f][j], sizeof(m_bColorArray[f][j]), 1, pFs);
            }
            
            // Write properties
            if(m_pTrkHdr->nbScalars > 6)
            {
                for(int i = 0; i < m_propertiesArray[f][j].size(); i++)
                {
                    fwrite(&m_propertiesArray[f][j][i], sizeof(m_propertiesArray[f][j][i]), 1, pFs);
                }
            }
        }
    }
        
    fclose(pFs);
    return true;
}
