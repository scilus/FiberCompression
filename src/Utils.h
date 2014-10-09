#ifndef UTILS_H
#define UTILS_H

#include <sstream>

#include <QFile>
#include <QString>
#include <QVector>

#include "Vector.h"

// Enums
enum FiberType
{
    FIBERTYPE_VTK = 0,
    FIBERTYPE_DMRI,
    FIBERTYPE_PTK,
    FIBERTYPE_CAMINO,
    FIBERTYPE_TRK,
    FIBERTYPE_MRTRIX,
    FIBERTYPE_NOT_SET,
};

enum EncodingType
{
    HUFFMAN_ENCODING = 0,
    ARITHMETIC_ENCODING,
    NO_ENCODING
};


enum TransformationType
{
    DCT_TRANSFORMATION = 0,
    DB4_TRANSFORMATION,
    DB6_TRANSFORMATION,
    DB8_TRANSFORMATION,
    BIOR_5_3_TRANSFORMATION,
    BIOR_9_7_TRANSFORMATION,
    NO_TRANSFORMATION
};

enum QuantizationType
{
    UNIFORM_QUANTIZATION = 0,
    NON_UNIFORM_QUANTIZATION,
    NO_QUANTIZATION
};


union converterByteINT16
{
    unsigned char b[2];
    unsigned int i;
};

union converterByteINT32
{
    unsigned char b[4];
    int i;
};

union converterByteFloat
{
    unsigned char b[4];
    float f;
};

template<typename T>
std::string valToString(const T& t)
{
    std::ostringstream os;
    os << t;
    return os.str();
}

float log2( float n );

const float DEFAULT_ERROR_MAX = 0.5;
const TransformationType DEFAULT_TRANSF_TYPE = NO_TRANSFORMATION;
const QuantizationType DEFAULT_QUANT_TYPE = UNIFORM_QUANTIZATION;
const int DEFAULT_PRECISION_QUANT = 0;
const int DEFAULT_ENCODE = 1;
const EncodingType DEFAULT_ENCODING_TYPE = HUFFMAN_ENCODING;
const FiberType DEFAULT_FIBER_TYPE = FIBERTYPE_NOT_SET;

const int MIN_UNIFORM = -3;
const int MAX_UNIFORM = 4;
const int DEFAULT_UNIFORM = 0;
const int MIN_NON_UNIFORM = 1;
const int MAX_NON_UNIFORM = 4;
const int DEFAULT_NON_UNIFORM = 3;
const float DEFAULT_MAX_ERROR = 0.5;

struct Param
{
    QString filename;
    bool compress;
    QString inputPath;
    QString outputPath;
    double errorMax;
    TransformationType transfType;
    QuantizationType quantizType;
    int precisionQuantiz;
    bool encode;
    EncodingType encodingType;
    QString statsOutputPath;
};

struct TRK_hdr
{
    char id_string[6];
    short dim[3];
    float voxelSize[3];
    float origin[3];
    short nbScalars;
    char scalarNames[10][20];
    short nbProperties;
    char propertyNames[10][20];
    float voxToRas[4][4];
    char reserved[444];
    char voxelOrder[4];
    char pad2[4];
    float imageOrientationPatient[6];
    char pad1[2];
    unsigned char invertX;
    unsigned char invertY;
    unsigned char invertZ;
    unsigned char swapXY;
    unsigned char swapYZ;
    unsigned char swapZX;
    int nbCount;
    int version;
    int hdrSize;
};

QVector<QString> getValidExtensions();
QString getStringValidExtensions();


// String from enums
QString getEncodingStringFromEnum(EncodingType type);
QString getTransformationStringFromEnum(TransformationType type);
QString getQuantizationStringFromEnum(QuantizationType type);

// Enum from string
EncodingType getEncodingEnumFromString(QString name);
TransformationType getTransformationEnumFromString(QString name);
QuantizationType getQuantizationEnumFromString(QString name);

bool validateParameters(Param& param, QString& errorTitle, QString& errorMsg);
void saveParameters(QFile& file, const Param& param);

// Utilities methods
QString getExtension(const QString filename);
QString getName(const QString filename);
int getFileSizeBytes(const QString filename);
float getCompressionRatio(const QString inputPath,
                          const QString outputPath);
float getPointToSegmentDistance(const Vector c,
                                const Vector a,
                                const Vector b);
QVector<Vector> mergeDimensions(const QVector<float>& arrayX,
                                const QVector<float>& arrayY,
                                const QVector<float>& arrayZ);


// Binary methods
QString int2binarystring(unsigned int x, unsigned int m);
int binary2int(const QString& s);
char bool2char(bool b);


#endif //UTILS_H