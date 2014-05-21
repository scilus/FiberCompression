#ifndef FIBERS_H
#define FIBERS_H

#include <QElapsedTimer>
#include <QFile>
#include <QList>
#include <QMatrix4x4>
#include <QObject>
#include <QString>
#include <QStatusBar>
#include <QVector>

#include "ArithmeticCoder.h"
#include "Huffman.h"
#include "Utils.h"
#include "Vector.h"


bool sortDescending(float a, float b);

class Fibers : public QObject
{
    Q_OBJECT
    
public:
    // Default constructor
    Fibers();
    
    // Default destructor
    ~Fibers();
    
    // Compress method
    void init(const bool compress,
              const QString& inputPath,
              const QString& outputPath,
              const QString& outputStatsFilename,
              const float errorMax = DEFAULT_ERROR_MAX,
              const TransformationType transformationType = DEFAULT_TRANSF_TYPE,
              const QuantizationType quantizationType = DEFAULT_QUANT_TYPE,
              const int precisionQuantization = DEFAULT_PRECISION_QUANT,
              const EncodingType encodingType = DEFAULT_ENCODING_TYPE);
    
public slots :
    void stop();
    void cleanUp();
    void reset();
    void process();
    
signals:
    void finished();
    void showMessage(const QString&);
    void isStopped();
    
private:
    
    // Compression methods
    bool compress();
    bool decompress();
    
    
    // Load methods
    bool load(const QString& filename);
    bool loadCompressed(const QString& filename,
                        const EncodingType& encodingType,
                        const TransformationType& ttype);
    bool loadDmri(const QString& filename);
    bool loadTRK(const QString& filename);
    bool loadTCK(const QString& filename);
    
    // Save methods
    bool save(const QString& filename);
    bool saveDMRI(const QString& filename);
    bool saveTCK(const QString& filename);
    bool saveTRK(const QString& filename);
    bool saveCompressed(const EncodingType etype,
                        const TransformationType ttype,
                        const QString& filename);
    bool saveResults(const QString& str);
    
    
    // Display methods
    void displayMessage(const QString& message);
    
    // Linearization methods
    bool linearize(const double thresholdValue);
    bool linearizeFiber(const int fibIdx,
                        const double thresholdValue,
                        QVector<unsigned char>& colorArrayX,
                        QVector<unsigned char>& colorArrayY,
                        QVector<unsigned char>& colorArrayZ,
                        QVector<float>& fiberResultX,
                        QVector<float>& fiberResultY,
                        QVector<float>& fiberResultZ);
    
    // Transformation methods
    bool getTransformation(const TransformationType type,
                           const QList< QVector<float> >& array,
                           const int dir,
                           QList< QVector<float> >& transformed);
    bool getArrayTransformation(const TransformationType type,
                                const QVector<float>& array,
                                const int dir,
                                const int origLength,
                                QVector<float>& result);
    
    // Approximation methods
    bool approximate(TransformationType ttype,
                     const float errorMaxT);
    QVector<float> nonLinearApproximation(const QVector<float>& array,
                                          const int nbKeptCoeffs);
    
    // Quantization methods
    float getQuantizationMaxError(QuantizationType type, int p);
    bool quantization(QList< QVector<float> >& array,
                      const QuantizationType type,
                      const int precision);
    bool uniformQuantization(QList< QVector<float> >& array,
                             const int precision);
    bool nonUniformQuantization(QList< QVector<float> >& array,
                                const int precision);
    float nonUniformRounding(const float value,
                             const int precision);
    
    // Max error methods
    bool getErrors(Fibers* original,
                   const QList< QVector<float> >& xPointsArray,
                   const QList< QVector<float> >& yPointsArray,
                   const QList< QVector<float> >& zPointsArray,
                   float& errorMax,
                   float& errorMean,
                   float& errorMeanMax,
                   float& meanStdError,
                   float& maxStdError);
    bool errorsBetweenFibers(const QVector<Vector>& original,
                             const QVector<Vector>& final,
                             float& errorMax,
                             float& errorMean,
                             float& errorStd);
    
    // Encoding methods
    bool encode(const QList< QVector<float> >& xArray,
                const QList< QVector<float> >& yArray,
                const QList< QVector<float> >& zArray,
                const EncodingType encodingType);
    
    
    // Decoding methods
    bool decode(const EncodingType encodingType);
    
    // Compression ratio methods
    float getCompressionRatio();
    
    // Manipulation methods
    void clearArrays();
    void setPointArray(const float value,
                       const int fibIdx,
                       const int ptIdx,
                       const int dimIdx);
    void setPointsArrayDim(const QList <QVector<float> > array,
                           const int i);
    void addFiberToPointArray(const QVector<Vector> array);
    void setFiberInPointArray(const QVector<float> array, const int fibIdx, const int dimIdx);
    void setColorArray(const QVector<unsigned char> array, const int fibIdx, const int dimIdx);
    
    // General methods
    void updateGlobalMinMax(float f);
    int getCountLines();
    int getLineSize(const int l);
    
    
    // Private Variables
    bool*           m_pStop;
    bool            m_colorArrayLoadedFromFile; // Flag that specifies if there was a color array loaded
    int             m_countLines;               // Number of fibers
    
    QString            m_outputStatsFilename;   // Text file to output stats during process
    QString            m_inputPath;             // Input file path
    QString            m_outputPath;            // Output file path
    float              m_errorMax;              // Maximum error
    TransformationType m_transformationType;    // Transformation type
    QuantizationType   m_quantizationType;      // Quantization type
    int                m_precisionQuantization; // Precision of quantization
    EncodingType       m_encodingType;          // Encoding type
    bool               m_compress;              // Flag indicating if we are compressing or decompressing
    
    FiberType                 m_fiberType;           // Type of fibers (DMRI, TCK, ...)
    Dictionnary               m_huffmanDict;         // Huffman dictionnary
    ColorDictionnary          m_colorHuffmanDict;    // Huffman color dictionnary
    ADictionnary<float>       m_arithmeticDict;      // Arithmetic dictionnary
    ADictionnary<unsigned char> m_colorArithmeticDict; // Arithmetic color dictionnary
    
    QList< QVector<unsigned char> > m_rColorArray;   // List of colors associated with fibers points
    QList< QVector<unsigned char> > m_gColorArray;   // List of colors associated with fibers points
    QList< QVector<unsigned char> > m_bColorArray;   // List of colors associated with fibers points
    QList< QVector< QVector<float> > > m_propertiesArray; // List of properties associated to fibers points (with trk)
    QList< QVector<float> >   m_xPointArray;         // List of fibers points in x dimension
    QList< QVector<float> >   m_yPointArray;         // List of fibers points in y dimension
    QList< QVector<float> >   m_zPointArray;         // List of fibers points in z dimension
    std::string               m_xEncoded;            // String containing all x coordinates encoded
    std::string               m_yEncoded;            // String containing all y coordinates encoded 
    std::string               m_zEncoded;            // String containing all z coordinates encoded 
    std::string               m_rEncodedColor;       // String containing all red colors encoded
    std::string               m_gEncodedColor;       // String containing all green colors encoded
    std::string               m_bEncodedColor;       // String containing all blue colors encoded
    QVector<int>              m_origLineLength;      // 1D vector containing the size of each fibers
    QFile                     m_resultFile;          // Statistics text file
    Fibers*                   m_pOriginal;
    
    TRK_hdr*                  m_pTrkHdr;
    
    float                     m_globalMin;
    float                     m_globalMax;
};


#endif