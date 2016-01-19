//
//  Copyright(c) 2013 Caroline Presseau
//  Copyright(c) 2013 SCIL @t University of Sherbrooke.
//
#include <assert.h>

#include <bitset>
#include <set>

#include <QFile>
#include <QFileInfo>
#include <QMessageBox>
#include <QTextStream>
#include <QVectorIterator>

#include "Utils.h"

// Supported extensions
QString validExtensions[3] = 
{
    QString(".tck"),
    QString(".trk"),
    QString(".zfib")
};
// Size of validExtensions (useful on manipulation)
int validExtensionsSize = 3;


float log2( float n )
{
    return log(n) / log(2.0);
}


/********************************************//**
\brief get vector containing all string extensions
\return vector of strings
***********************************************/
QVector<QString> getValidExtensions()
{
    QVector<QString> validExt;
    for(int i = 0; i < validExtensionsSize; i++)
    {
        validExt.push_back(validExtensions[i]);
    }
    return validExt;
}


/********************************************//**
\brief get string containing all valid extensions
\return string
***********************************************/
QString getStringValidExtensions()
{
    QString validExt;
    for(int i = 0; i < validExtensionsSize - 1; i++)
    {
        validExt += validExtensions[i] + ", ";
    }
    validExt += validExtensions[validExtensionsSize - 1];
    return validExt;
}



/********************************************//**
\brief get string corresponding to a specific transformation enum
\param type : transformationType to convert
\return string
***********************************************/
QString getTransformationStringFromEnum(TransformationType type)
{
    switch(type)
    {
        case DCT_TRANSFORMATION: return QString("DCT");
        case DB4_TRANSFORMATION: return QString("Daubechies 4");
        case DB6_TRANSFORMATION: return QString("Daubechies 6");
        case DB8_TRANSFORMATION: return QString("Daubechies 8");
        case BIOR_5_3_TRANSFORMATION: return QString("Biorthogonal 5/3");
        case BIOR_9_7_TRANSFORMATION: return QString("Biorthogonal 9/7");
        case NO_TRANSFORMATION: return QString("No transformation");
        default: return QString("");
    }
}


/********************************************//**
\brief get transformation enum corresponding to a specific string
\param name : string to convert in enum
\return TransformationType enum
***********************************************/
TransformationType getTransformationEnumFromString(QString name)
{
    if(name ==  getTransformationStringFromEnum(DCT_TRANSFORMATION)) return DCT_TRANSFORMATION;
    else if(name ==  getTransformationStringFromEnum(DB4_TRANSFORMATION)) return DB4_TRANSFORMATION;
    else if(name ==  getTransformationStringFromEnum(DB6_TRANSFORMATION)) return DB6_TRANSFORMATION;
    else if(name ==  getTransformationStringFromEnum(DB8_TRANSFORMATION)) return DB8_TRANSFORMATION;
    else if(name ==  getTransformationStringFromEnum(BIOR_5_3_TRANSFORMATION)) return BIOR_5_3_TRANSFORMATION;
    else if(name ==  getTransformationStringFromEnum(BIOR_9_7_TRANSFORMATION)) return BIOR_9_7_TRANSFORMATION;
    else return NO_TRANSFORMATION;
}


/********************************************//**
\brief get string corresponding to a specific quantization enum
\param type : quantizationType to convert
\return string
***********************************************/
QString getQuantizationStringFromEnum(QuantizationType type)
{
    switch(type)
    {
        case UNIFORM_QUANTIZATION: return QString("Uniform");
        case NON_UNIFORM_QUANTIZATION: return QString("Non-uniform");
        case NO_QUANTIZATION: return QString("No quantization");
        default: return QString("");
    }
}


/********************************************//**
\brief get quantization enum corresponding to a specific string
\param name : string to convert in enum
\return QuantizationType enum
***********************************************/
QuantizationType getQuantizationEnumFromString(QString name)
{
    if(name ==  getQuantizationStringFromEnum(UNIFORM_QUANTIZATION)) return UNIFORM_QUANTIZATION;
    else if(name ==  getQuantizationStringFromEnum(NON_UNIFORM_QUANTIZATION)) return NON_UNIFORM_QUANTIZATION;
    else return NO_QUANTIZATION;
}


/********************************************//**
\brief get string corresponding to a specific encoding enum
\param type : encodingType to convert
\return string
***********************************************/
QString getEncodingStringFromEnum(EncodingType type)
{
    switch(type)
    {
        case HUFFMAN_ENCODING: return QString("Huffman");
        case ARITHMETIC_ENCODING: return QString("Arithmetic");
        case NO_ENCODING: return QString("No encoding");
        default: return QString("");
    }
}


/********************************************//**
\brief get encoding enum corresponding to a specific string
\param name : string to convert in enum
\return EncodingType enum
***********************************************/
EncodingType getEncodingEnumFromString(QString name)
{
    if(name ==  getEncodingStringFromEnum(HUFFMAN_ENCODING)) return HUFFMAN_ENCODING;
    else if(name ==  getEncodingStringFromEnum(ARITHMETIC_ENCODING)) return ARITHMETIC_ENCODING;
    else return NO_ENCODING;
}


/********************************************//**
\brief validate parameters used in the process before starting pipeline
\param param : structure containing parameters to save
\return boolean true if parameters are all valids
***********************************************/
bool validateParameters(Param& param, QString& errorTitle, QString& errorMsg)
{
    // Check if input path and all extensions are valid
    QString inputExt = getExtension(param.inputPath);
    QString outputExt = getExtension(param.outputPath);
    
    if(param.inputPath.isEmpty())
    {
        errorTitle = QString("Input file is empty");
        errorMsg = QString("Input file is empty. Please provide an input file.");
        return false;
    }
    if(param.outputPath.isEmpty())
    {
        errorTitle = QString("Output file is empty");
        errorMsg = QString("Output file is empty. Please provide an output file.");
        return false;
    }
    
    if(outputExt.isEmpty())
    {
        if(param.compress)
        {
            param.outputPath += QString(".zfib");
        }
    }
    
    // Validate if input file exists and if extensions are valid
    if(!QFile::exists(param.inputPath))
    {
        errorTitle = QString("File doesn't exist");
        errorMsg = QString("Input file doesn't exist. Please provide a valid file.");
        return false;
    }
    
    if(param.compress)
    {
        if(param.encodingType == NO_ENCODING && inputExt != outputExt)
        {
            errorTitle = QString("Invalid extension");
            errorMsg = QString("If you don't want to encode the file, output extension must be the same as input.");
            return false;
        }
        else if(param.transfType != NO_TRANSFORMATION && outputExt != QString("zfib"))
        {
            errorTitle = QString("Invalid extension");
            errorMsg = QString("If you want to use a transformation, output extension must be .zfib.");
            return false;
        }
        else
        {
            if(inputExt != QString("tck") && inputExt != QString("trk"))
            {
                errorTitle = QString("Invalid extension");
                errorMsg = QString("Invalid input file's extension.");
                return false;
            }
            if(inputExt != outputExt && inputExt != QString("zfib") && outputExt != QString("zfib"))
            {
                errorTitle = QString("Invalid extension");
                errorMsg = QString("Invalid output file's extension. If you don't want to save in .zfib, you must provide the same extension as input file.");
                return false;
            }
            if(param.encodingType != NO_ENCODING && outputExt != QString("zfib") && inputExt != outputExt)
            {
                errorTitle = QString("Invalid extension");
                errorMsg = QString("Invalid output file's extension. Output file's extension must be .zfib.");
                return false;
            }
        }
    }
    else
    {
        if(inputExt != QString("zfib"))
        {
            errorTitle = QString("Invalid extension");
            errorMsg = QString("Invalid input file's extension. Input file's extension must be .zfib.");
            return false;
        }
        if(outputExt != QString("tck") && outputExt != QString("trk"))
        {
            errorTitle = QString("Invalid extension");
            errorMsg = QString("Invalid output file's extension.");
            return false;
        }
    }
    return true;
}


/********************************************//**
\brief save parameters used in the process in a text file
\param myFile : instance of QFile to write to
\param param : structure containing parameters to save
***********************************************/
void saveParameters(QFile& myFile, const Param& param)
{
    QTextStream out(&myFile);
    
    if(myFile.open(QIODevice::WriteOnly | QIODevice::Text))
    {
        if(param.compress)
        {
            out << "Compression\n";
            out << "Input file : " + param.inputPath + "\n";
            out << "Output file : " + param.outputPath + "\n";
            out << "Error max : " + QString::number(param.errorMax) + " mm\n";
            out << "Transformation type : " + getTransformationStringFromEnum(param.transfType) + "\n";
            out << "Quantization type : " +
            getQuantizationStringFromEnum(param.quantizType) + "  with " +
            QString::number(param.precisionQuantiz) + " digits\n";
            out << "Encoding type : " + getEncodingStringFromEnum(param.encodingType) + "\n";
        }
        else
        {
            out << "Decompression\n";
            out << "Input file : " + param.inputPath + "\n";
            out << "Output file : " + param.outputPath + "\n";
        }
        out << "\n";
    }
    
    myFile.close();
}


/********************************************//**
\brief Get extension of a file path
\param filename : input file
\return extension without the .
***********************************************/
QString getExtension(const QString filename)
{
    QFileInfo file(filename);
    return file.suffix();
}


/********************************************//**
\brief Get name of a file path
\param filename : input file
\return name before .extension
***********************************************/
QString getName(const QString filename)
{
    QFileInfo file(filename);
	QString path = file.absoluteFilePath();
	path.remove("." + file.suffix());
    return path;
}

/********************************************//**
\brief Compute compression ratio between two files
\param inputPath : first file path
\param outputPath : second file path
\return compression ratio
***********************************************/
float getCompressionRatio(const QString inputPath,
                          const QString outputPath)
{
    int origSizeBytes = getFileSizeBytes(inputPath);
    int compSizeBytes = getFileSizeBytes(outputPath);
    return (float)(100.0f *((origSizeBytes - compSizeBytes) / (float)origSizeBytes));
}


/********************************************//**
\brief Get size of a file in bytes
\param filename : file path
\return size in bytes
***********************************************/
int getFileSizeBytes(const QString filename)
{
    int fileSizeBytes = 0;
    QFile myfile(filename);
    
    if(myfile.open(QIODevice::ReadOnly))
    {
        fileSizeBytes = (int)myfile.size();
    }
    myfile.close();
    
    return fileSizeBytes;
}


/********************************************//**
\brief Compute distance from a point (c) to a segment (ab)
\param c : point
\param a : first point of the segment
\param b : second point of the segment
\return distance
***********************************************/
float getPointToSegmentDistance(const Vector c,
                                const Vector a,
                                const Vector b)
{
    // Compute the distance from a point c to a line segment ab
    Vector ab = b - a;
    Vector ac = c - a;
    
    float ab_magn = ab.getMagnitude();
    if(ab_magn == 0.0f) // Distance point to point
    {
        // TODO HERE CHECK THIS AND THE CONDITION
        Vector d = ab - ac;
        return d.getMagnitude();
    }
    
    // TODO CHECK THIS FORMULA
    return (ab.cross(ac)).getMagnitude() / ab_magn;
}


/********************************************//**
\brief Merge 3 1D vectors in a 1D vector with 3 dimensions
\param arrayX : input 1D array corresponding to the first position
\param arrayY : input 1D array corresponding to the second position
\param arrayZ : input 1D array corresponding to the third position
\return merge vector
***********************************************/
QVector<Vector> mergeDimensions(const QVector<float>& arrayX,
                                const QVector<float>& arrayY,
                                const QVector<float>& arrayZ)
{
    assert(arrayX.size() == arrayY.size() && arrayX.size() == arrayZ.size());
    
    QVector<Vector> result(arrayX.size());
    for(int i = 0; i < arrayX.size(); i++)
    {
        result[i] = Vector(arrayX[i], arrayY[i], arrayZ[i]);
    }
    return result;
}


/********************************************//**
\brief Convert unsigned integer to binary string
\param x : unsigned integer to convert
\param m : unsigned integer which is the total length 
           wanted for the result
\return binary string
***********************************************/
QString int2binarystring(unsigned int x, unsigned int m)
{
    QString res;
    std::string s;
    // Convert x to binary string
    do
    {
        s.push_back('0' + (x & 1));
    } while (x >>= 1);
    std::reverse(s.begin(), s.end());
    
    while(s.size() < m)
    {
        s = '0' + s;
    }
    return res.fromStdString(s);
}


/********************************************//**
\brief Convert binary string to integer
\param s : binary string to convert
\return integer
***********************************************/
int binary2int(const QString& s)
{
    int result = 0;
    int k = 0;
    
    // Convert binary string to integer
    for(int i = s.length() - 1; i >= 0; i--)
    {
        QChar bit = s[i];
        if((bit != '0') && (bit != '1'))
        {
            return -1;
        }
        result += (bit == '1') * pow(2.0, k);
        k++;
    }
    return result;
}


/********************************************//**
\brief Convert boolean to char
\param b : boolean value to convert
\return char
***********************************************/
char bool2char(bool b)
{
    if(b)
    {
        return '1';
    }
    else
    {
        return '0';
    }
}