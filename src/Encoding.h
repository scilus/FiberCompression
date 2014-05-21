#ifndef ENCODING_H
#define ENCODING_H

#include <algorithm>
#include <iostream>

#include "ArithmeticCoder.h"
#include "Huffman.h"
#include "HuffmanUtilsIO.h"


/********************************************//**
\brief Encode data with Huffman algorithm
\param array : input array to encode
\param dict : dictionnary that will be filled
\return encoded array with Huffman
***********************************************/
template<typename Type, typename DictionnaryType>
bool huffmanEncoding(const QList< QVector<Type> >& xArray,
                     const QList< QVector<Type> >& yArray,
                     const QList< QVector<Type> >& zArray,
                     DictionnaryType& dict,
                     std::string& xEncoded,
                     std::string& yEncoded,
                     std::string& zEncoded,
                     bool* stop)
{
    xEncoded = yEncoded = zEncoded = "";
    assert(xArray.size() == yArray.size() && xArray.size() == zArray.size());

    // Fill Huffman Dictionnary
    if(!fillHuffmanSymbolsAndProbabilities(xArray, yArray, zArray, dict, stop))
    {
        return false;
    }

    // Create tree
    Hufftree<Type, float> hufftree(dict.frequencies.begin(), dict.frequencies.end(), stop);
    
    for(int i = 0; i < xArray.size(); i++)
    {
        // Get x indexes and convert them to bytes
        QVector<int> indexToBeEncoded;
        if(!symbolListToIndexList(dict, xArray[i], indexToBeEncoded, stop))
        {
            return false;
        }
        std::string encoded;
        if(!hufftree.encode(indexToBeEncoded.begin(), indexToBeEncoded.end(), encoded, stop))
        {
            return false;
        }
        xEncoded += encoded;
        
        // Get y indexes and convert them to bytes
        indexToBeEncoded.clear();
        if(!symbolListToIndexList(dict, yArray[i], indexToBeEncoded, stop))
        {
            return false;
        }
        encoded = "";
        if(!hufftree.encode(indexToBeEncoded.begin(), indexToBeEncoded.end(), encoded, stop))
        {
            return false;
        }
        yEncoded += encoded;
        
        // Get z indexes and convert them to bytes
        indexToBeEncoded.clear();
        if(!symbolListToIndexList(dict, zArray[i], indexToBeEncoded, stop))
        {
            return false;
        }
        encoded = "";
        if(!hufftree.encode(indexToBeEncoded.begin(), indexToBeEncoded.end(), encoded, stop))
        {
            return false;
        }
        zEncoded += encoded;
    }

    return true;
}


/********************************************//**
\brief Decode data encoded with Huffman algorithm
\param dict : dictionnary used to decode
\param encoded : input array to decode
\return decoded array with Huffman
***********************************************/
template<typename Type, typename DictionnaryType>
bool huffmanDecoding(DictionnaryType& dict,
                     const std::string& xEncoded,
                     const std::string& yEncoded,
                     const std::string& zEncoded,
                     QList< QVector<Type> >& xDecoded,
                     QList< QVector<Type> >& yDecoded,
                     QList< QVector<Type> >& zDecoded,
                     int countLines,
                     QVector<int> lineLength,
                     bool* stop)
{
    // Create tree
    Hufftree<Type, float> hufftree(dict.frequencies.begin(), dict.frequencies.end(), stop);
    if(hufftree.isStopped())
    {
        return false;
    }
    
    // Decode index symbols
    QVector<int> xIndexList, yIndexList, zIndexList;
    if(!hufftree.decode(xEncoded, std::back_inserter(xIndexList), stop))
    {
        return false;
    }
    if(!hufftree.decode(yEncoded, std::back_inserter(yIndexList), stop))
    {
        return false;
    }
    if(!hufftree.decode(zEncoded, std::back_inserter(zIndexList), stop))
    {
        return false;
    }
    
    // Convert index lists to symbols list
    QVector<Type> xDecodedList, yDecodedList, zDecodedList;
    if(!indexListToSymbolList<Type>(dict, xIndexList, xDecodedList, stop))
    {
        return false;
    }
    if(!indexListToSymbolList<Type>(dict, yIndexList, yDecodedList, stop))
    {
        return false;
    }
    if(!indexListToSymbolList<Type>(dict, zIndexList, zDecodedList, stop))
    {
        return false;
    }
    
    assert(xDecodedList.size() > 0 && yDecodedList.size() > 0 && zDecodedList.size() > 0);
    xDecoded.clear();
    yDecoded.clear();
    zDecoded.clear();
    
    // Fit 1D vector of each dimension in m_xPointArray, m_yPointArray and m_zPointArray
    int idx = 0;
    for(int l = 0; l < countLines; l++)
    {
        // Check flag to know if we need to stop the process or not
        if(*stop)
        {
            return false;
        }
        
        int length = lineLength[l];
        QVector<Type> xArray(length), yArray(length), zArray(length);
        
        for(int i = 0; i < length; i++)
        {
            xArray[i] = xDecodedList[idx];
            yArray[i] = yDecodedList[idx];
            zArray[i] = zDecodedList[idx];
            idx ++;
        }
        xDecoded.append(xArray);
        yDecoded.append(yArray);
        zDecoded.append(zArray);
    }
    
    return true;
}


/********************************************//**
\brief Compute symbols and probabilities useful to Huffman dictionnary
\param array : array to be encoded
\param symbols : array of symbols that will be filled
\param probabilities : array of probabilities that will be filled
***********************************************/
template<typename Type, typename DictionnaryType>
bool fillHuffmanSymbolsAndProbabilities(const QList< QVector<Type> >& xArray,
                                        const QList< QVector<Type> >& yArray,
                                        const QList< QVector<Type> >& zArray,
                                        DictionnaryType& dict,
                                        bool* stop)
{
    assert(xArray.size() == yArray.size() && xArray.size() == zArray.size());
    int totalSize = 0;
    
    for(int i = 0; i < xArray.size(); i++)
    {
        assert(xArray[i].size() == yArray[i].size() && xArray[i].size() == zArray[i].size());
        
        if(*stop)
        {
            return false;
        }
        QVector<Type> symbols = xArray[i];
        symbols += yArray[i];
        symbols += zArray[i];
        
        // Sort symbols in ascending order
        std::sort(symbols.begin(), symbols.end());
        
        // Remove non-unique elements in symbols
        typename QVector<Type>::iterator it = std::unique(symbols.begin(), symbols.end());
        symbols.resize(it - symbols.begin());
        
        // Find symbols
        Type size = (Type)xArray[i].size();
        totalSize += 3 * size;
        for(int k = 0; k < symbols.size(); k++)
        {
            float prob = xArray[i].count(symbols[k]) +
                         yArray[i].count(symbols[k]) +
                         zArray[i].count(symbols[k]);
            
            // Symbol found
            //std::map<int, float>::iterator it = std::find(dict.symbols.begin(), dict.symbols.end(), symbols[k]);
            int symIdx = dict.symbols.indexOf(symbols[k]);
            if(symIdx != -1 )
            {
                std::map<int, float>::iterator it = dict.frequencies.find(symIdx);
                assert(it != dict.frequencies.end());
                it->second += prob;
            }
            else
            {
               dict.symbols.push_back(symbols[k]);
               int idx = dict.symbols.size() - 1;
               dict.frequencies.insert(std::pair<int, float>(idx, prob));
            }
        }

    }
    for(std::map<int, float>::iterator it = dict.frequencies.begin(); it != dict.frequencies.end(); it++)
    {
        it->second /= totalSize;
    }

    return true;
}


/********************************************//**
\brief Encode array with Arithmetic algorithm
\param array : input array to encode
\param dict : dictionnary that will be filled
\return encoded array with arithmetic
***********************************************/
template<typename Type>
bool arithmeticEncoding(const QList< QVector<Type> >& xArray,
                        const QList< QVector<Type> >& yArray,
                        const QList< QVector<Type> >& zArray,
                        ADictionnary<Type>& dict,
                        std::string& encoded,
                        bool* stop)
{
    ArithmeticCoder<Type> arithmetic;
    
    // Fill dictionnary
    if(!arithmetic.fillDictionnary(xArray, yArray, zArray, dict, stop))
    {
        return false;
    }
    
    // Encode array
    if(!arithmetic.encode(xArray, yArray, zArray, encoded, stop))
    {
        
    }
    return true;
}


/********************************************//**
\brief Decode array with Arithmetic algorithm
\param dict : dictionnary used to decode
\param encoded : input array to decode
\return decoded array with arithmetic
***********************************************/
template<typename Type>
bool arithmeticDecoding(ADictionnary<Type>& dict,
                        const std::string& encoded,
                        QList< QVector<Type> >& xDecoded,
                        QList< QVector<Type> >& yDecoded,
                        QList< QVector<Type> >& zDecoded,
                        int countLines,
                        QVector<int> lineLength,
                        bool* stop)
{
    ArithmeticCoder<Type> arithmetic;
    
    // Set dictionnary
    arithmetic.setDictionnary(dict);
    
    // Decode symbols
    QVector<Type> decodedList;
    if(!arithmetic.decode(encoded, decodedList, stop))
    {
        return false;
    }
    
    assert(decodedList.size() > 0);
    xDecoded.clear();
    yDecoded.clear();
    zDecoded.clear();
    
    // Fit 1D vector of each dimension in m_xPointArray, m_yPointArray and m_zPointArray
    int idx = 0;
    for(int l = 0; l < countLines; l++)
    {
        // Check flag to know if we need to stop the process or not
        if(*stop)
        {
            return false;
        }
        
        int length = lineLength[l];
        QVector<Type> xArray(length), yArray(length), zArray(length);
        
        for(int i = 0; i < length; i++)
        {
            xArray[i] = decodedList[idx];
            yArray[i] = decodedList[idx + length];
            zArray[i] = decodedList[idx + (2 * length)];
            idx ++;
        }
        idx += 2 * length;
        xDecoded.append(xArray);
        yDecoded.append(yArray);
        zDecoded.append(zArray);
    }

    return true;
}


/********************************************//**
\brief Write Huffman encoded signal and its dictionnary
\param myfile : file to write into
\param encoded : input array to write
\param dict : dictionnary to write
***********************************************/
template<typename Type, typename DictionnaryType>
bool writeHuffmanEncodedSignal(std::ofstream& myfile,
                               const std::string& xEncoded,
                               const std::string& yEncoded,
                               const std::string& zEncoded,
                               DictionnaryType& dict,
                               bool* stop)
{
    int xEncodedSize = xEncoded.size();
    int yEncodedSize = yEncoded.size();
    int zEncodedSize = zEncoded.size();
    int dictSize = dict.frequencies.size();
    
    // Write size of encoded signals and size of dictionnary
    myfile.write((const char*) &(xEncodedSize), sizeof(xEncodedSize));
    myfile.write((const char*) &(yEncodedSize), sizeof(yEncodedSize));
    myfile.write((const char*) &(zEncodedSize), sizeof(zEncodedSize));
    myfile.write((const char*) &(dictSize), sizeof(dictSize));
    
    // Save symbols and frequencies of dictionary
    for(int i = 0; i < dictSize; i++)
    {
        if(*stop)
        {
            return false;
        }
        Type val = dict.symbols[i];
        myfile.write((const char*) &(val), sizeof(val));
        float valFreq = dict.frequencies[i];
        myfile.write((const char*) &(valFreq), sizeof(valFreq));
    }
    
    // Save x encoded signal
    int iter = 0;
    std::bitset<8>bits;
    bits.reset();
    unsigned char val;
    for(int i = 0; i < xEncodedSize; i++)
    {
        if(*stop)
        {
            return false;
        }
        bits.set(iter, (int)(xEncoded[i]));
        iter=(iter+1)%8;
        
        if(iter==0)
        {
            val = (unsigned char)(bits.to_ulong());
            myfile.write((const char*) &(val), sizeof(val));
            bits.reset();
        }
    }
    if(iter!=0)
    {
        val = (unsigned char)(bits.to_ulong());
        myfile.write((const char*) &(val), sizeof(val));
    }
    
    // Save y encoded signal
    iter = 0;
    bits.reset();
    for(int i = 0; i < yEncodedSize; i++)
    {
        if(*stop)
        {
            return false;
        }
        bits.set(iter, (int)(yEncoded[i]));
        iter=(iter+1)%8;
        
        if(iter==0)
        {
            val = (unsigned char)(bits.to_ulong());
            myfile.write((const char*) &(val), sizeof(val));
            bits.reset();
        }
    }
    if(iter!=0)
    {
        val = (unsigned char)(bits.to_ulong());
        myfile.write((const char*) &(val), sizeof(val));
    }
    
    // Save z encoded signal
    iter = 0;
    bits.reset();
    for(int i = 0; i < zEncodedSize; i++)
    {
        if(*stop)
        {
            return false;
        }
        bits.set(iter, (int)(zEncoded[i]));
        iter=(iter+1)%8;
        
        if(iter==0)
        {
            val = (unsigned char)(bits.to_ulong());
            myfile.write((const char*) &(val), sizeof(val));
            bits.reset();
        }
    }
    if(iter!=0)
    {
        val = (unsigned char)(bits.to_ulong());
        myfile.write((const char*) &(val), sizeof(val));
    }

    return true;
}


/********************************************//**
\brief Read Huffman encoded signal and fill its dictionnary
\param myfile : file to read from
\param encoded : encoded array to fill during read
\param dict : dictionnary to fill during read
\return true if everything went well, false otherwise
***********************************************/
template<typename Type, typename DictionnaryType>
bool readHuffmanEncodedSignal(std::ifstream& myfile,
                              std::string& xEncoded,
                              std::string& yEncoded,
                              std::string& zEncoded,
                              DictionnaryType& dict,
                              bool* stop)
{
    if(myfile.eof() || myfile.fail())
    {
        return false;
    }
    
    int xEncodedSize = -1;
    int yEncodedSize = -1;
    int zEncodedSize = -1;
    int dictSize = -1 ;
    
    // Read size of encoded signal and size of dictionnary
    myfile.read((char*)&xEncodedSize, sizeof(xEncodedSize));
    myfile.read((char*)&yEncodedSize, sizeof(yEncodedSize));
    myfile.read((char*)&zEncodedSize, sizeof(zEncodedSize));
    myfile.read((char*)&dictSize, sizeof(dictSize));
    
    if(xEncodedSize < 0 || yEncodedSize < 0 || zEncodedSize < 0 || dictSize < 0)
    {
        return false;
    }
    
    // Load symbols and frequencies
    for(int i = 0; i < dictSize; i++)
    {
        if(*stop)
        {
            return false;
        }
        
        Type val;
        myfile.read((char*)&val, sizeof(val));
        dict.symbols.push_back(val);
        
        float valFreq;
        myfile.read((char*)&valFreq, sizeof(valFreq));
        dict.frequencies[i] = valFreq;
    }
    
    // Load x encoded signal
    int iter = 0;
    xEncoded = yEncoded = zEncoded = "";
    while(iter < xEncodedSize)
    {
        if(*stop)
        {
            return false;
        }
        
        char valRaw;
        unsigned char val;
        myfile.read(&valRaw, 1);
        val = (unsigned char)valRaw;
        std::bitset<8>bits(val);
        for(int i = 0; i < 8; i++)
        {
            iter++;
            xEncoded.push_back(bits[i]);
            if(iter >= xEncodedSize)
                break;
        }
    }
    
    // Load y encoded signal
    iter = 0;
    while(iter < yEncodedSize)
    {
        if(*stop)
        {
            return false;
        }
        
        char valRaw;
        unsigned char val;
        myfile.read(&valRaw, 1);
        val = (unsigned char)valRaw;
        std::bitset<8>bits(val);
        for(int i = 0; i < 8; i++)
        {
            iter++;
            yEncoded.push_back(bits[i]);
            if(iter >= yEncodedSize)
                break;
        }
    }
    
    // Load z encoded signal
    iter = 0;
    while(iter < zEncodedSize)
    {
        if(*stop)
        {
            return false;
        }
        
        char valRaw;
        unsigned char val;
        myfile.read(&valRaw, 1);
        val = (unsigned char)valRaw;
        std::bitset<8>bits(val);
        for(int i = 0; i < 8; i++)
        {
            iter++;
            zEncoded.push_back(bits[i]);
            if(iter >= zEncodedSize)
                break;
        }
    }

    return true;
}


/********************************************//**
\brief Write arithmetic dictionnary
\param myfile : file to write into
\param dict : dictionnary to write
***********************************************/
template<typename Type>
inline bool writeArithmeticDictionnary(std::ofstream& myfile,
                                       ADictionnary<Type>& dict,
                                       bool* stop)
{
    // Save size of dictionnary
    int dictSize = dict.cum_counts.size();
    myfile.write((const char*) &(dictSize), sizeof(dictSize));
    
    // Write symbols and cum_counts
    for (int i = 0; i < dictSize; i++)
    {
        if(*stop)
        {
            return false;
        }
        Type symbol = (dict.cum_counts[i].first);
        myfile.write((const char*) &(symbol), sizeof(symbol));
        int cum_count = (dict.cum_counts[i].second);
        myfile.write((const char*) &(cum_count), sizeof(cum_count));;
    }
    return true;
}


/********************************************//**
\brief Write arithmetic encoded signal
\param myfile : file to write into
\param encoded : input array to write
\param dict : dictionnary to write
***********************************************/
template<typename Type>
bool writeArithmeticEncodedSignal(std::ofstream& myfile,
                                  const std::string& encoded,
                                  ADictionnary<Type>& dict,
                                  bool* stop)
{
    int encodedSize = encoded.size();
    int m = dict.m;
    int nbSymbolsEncoded = dict.nbSymbolsEncoded;
    int totalCount = dict.totalCount;
    
    // Write size of encoded signal, m, the number of symbols encoded and the total count
    myfile.write((const char*) &(encodedSize), sizeof(encodedSize));
    myfile.write((const char*) &(m), sizeof(m));
    myfile.write((const char*) &(nbSymbolsEncoded), sizeof(nbSymbolsEncoded));
    myfile.write((const char*) &(totalCount), sizeof(totalCount));
    
    // Save dictionary
    writeArithmeticDictionnary(myfile, dict, stop);
    
    // Save the encoded signal
    int iter = 0;
    std::bitset<8>bits;
    bits.reset();
    unsigned char val;
    for(int i = 0; i < encodedSize; i++)
    {
        if(*stop)
        {
            return false;
        }
        
        bits.set(iter, (int)(encoded[i]));
        iter=(iter+1)%8;
        
        if(iter==0)
        {
            val = (unsigned char)(bits.to_ulong());
            myfile.write((const char*) &(val), sizeof(val));
            bits.reset();
        }
    }
    if(iter!=0)
    {
        val = (unsigned char)(bits.to_ulong());
        myfile.write((const char*) &(val), sizeof(val));
    }
    
    return true;
}


/********************************************//**
\brief Read arithmetic dictionnary
\param myfile : file to read from
\param dict : dictionnary to fill
\return true if everything went well, false otherwise
***********************************************/
template<typename Type>
bool readArithmeticDictionnary(std::ifstream& myfile,
                               ADictionnary<Type>& dict,
                               bool* stop)
{
    int dictSize = -1;
    
    // Read dictionnary size
    myfile.read((char*)&dictSize, sizeof(dictSize));
    if(dictSize < 0)
    {
        return false;
    }
    
    // Load symbols and cum_counts
    for (int i = 0; i < dictSize; i++)
    {
        if(*stop)
        {
            return false;
        }
        
        Type symbol;
        myfile.read((char*) &(symbol), sizeof(symbol));
        int cum_count;
        myfile.read((char*) &(cum_count), sizeof(cum_count));
        dict.cum_counts.push_back(QPair<Type, int>(symbol, cum_count));
    }
    return true;
}


/********************************************//**
\brief Read arithmetic encoded signal and fill its dictionnary
\param myfile : file to read from
\param encoded : encoded array to fill during read
\param dict : dictionnary to fill during read
\return true if everything went well, false otherwise
***********************************************/
template<typename Type>
bool readArithmeticEncodedSignal(std::ifstream& myfile,
                                 std::string& encoded,
                                 ADictionnary<Type>& dict,
                                 bool* stop)
{
    if(myfile.eof() || myfile.fail())
    {
        return false;
    }
    
    int encodedSize = -1;
    int m = -1;
    int nbSymEncoded = -1;
    int totalCount = -1;
    
    // Read size of encoded signal, m, the number of symbols encoded and the total count
    myfile.read((char*)&encodedSize, sizeof(encodedSize));
    myfile.read((char*)&(m), sizeof(m));
    myfile.read((char*)&(nbSymEncoded), sizeof(nbSymEncoded));
    myfile.read((char*)&(totalCount), sizeof(totalCount));
    
    if(encodedSize < 0 || m < 0 || nbSymEncoded < 0 || totalCount < 0)
    {
        return false;
    }
    
    // Load dictionnary
    if(!readArithmeticDictionnary(myfile, dict, stop))
    {
        return false;
    }
    
    dict.m = m;
    dict.nbSymbolsEncoded = nbSymEncoded;
    dict.totalCount = totalCount;
    
    // Load the encoded signal
    int iter = 0;
    encoded.clear();
    
    while(iter < encodedSize)
    {
        if(*stop)
        {
            return false;
        }
        
        char valRaw;
        unsigned char val;
        myfile.read(&valRaw, 1);
        val = (unsigned char)valRaw;
        std::bitset<8>bits(val);
        for(int i = 0; i < 8; i++)
        {
            iter++;
            encoded.push_back(bits[i]);
            if(iter >= encodedSize)
                break;
        }
    }
    
    return true;
}

#endif