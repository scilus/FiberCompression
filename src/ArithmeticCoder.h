#ifndef ARITHMETICCODER_H
#define ARITHMETICCODER_H

#include <algorithm>
#include <cstdlib>
#include <math.h>
#include <set>
#include <vector>

#include <QPair>
#include <QSet>
#include <QString>
#include <QVector>

#include "Utils.h"


union converterByteChar
{
    unsigned char b[4];
    unsigned int i;
};


/********************************************//**
\brief Operator < for QPair<Type, float>
\param left : First pair to compare
\param right : Second pair to compare
\return true if < or false otherwise
***********************************************/
template<typename Type>
bool operator<(const QPair<Type, float>& left,
               const QPair<Type, float>& right)
{
    return left.second < right.second;
}
    
    
/********************************************//**
\brief Method that indicates if two pairs are equals
\return true if pairs are equal or false otherwise
***********************************************/
template <class T,class S> struct pair_equal_to : std::binary_function <T,QPair<T,S>,bool>
{
    bool operator() (const T& y, const QPair<T,S>& x) const
    {
        return x.first == y;
    }
};
    
    
// Arithmetic Dictionnary structure
template<typename Type>
struct ADictionnary
{
    QVector< QPair<Type, int> > cum_counts;     // cumulative counts
    int totalCount;                             // total count
    int m;                                      // number of bytes
    int nbSymbolsEncoded;                       // number of symbols to encode
        
    void cleanUp()
    {
        cum_counts.clear();
        totalCount = 0;
        m = 0;
        nbSymbolsEncoded = 0;
    }
};
    
    
template<typename Type>
class ArithmeticCoder
{
public:
    // Default constructor
    ArithmeticCoder(){};
    
    // Default destructor
    ~ArithmeticCoder(){};
    
    // Set dictionnary
    bool fillDictionnary(const QList< QVector<Type> >& xArray,
                         const QList< QVector<Type> >& yArray,
                         const QList< QVector<Type> >& zArray,
                         ADictionnary<Type>& dict,
                         bool* stop);
    void setDictionnary(const ADictionnary<Type>& dict);
    
    // Encoding
    bool encode(const QList< QVector<Type> >& xArray,
                const QList< QVector<Type> >& yArray,
                const QList< QVector<Type> >& zArray,
                std::string& encoded,
                bool* stop);
    
    bool encodeSymbol(const QVector<Type>& array,
                      int& l,
                      int& u,
                      int& scaleE3,
                      QString& l_bin,
                      std::string& encoded,
                      bool* stop);
    
    // Decoding
    bool decode(const std::string& encoded, QVector<Type>& decoded, bool* stop);
    
private:
    QVector<Type> removeDuplicatesFromVector(const QVector<Type>& input);
    
    // Private variables
    ADictionnary<Type> m_dictionnary;       // dictionnary
};
    
    
/********************************************//**
\brief Set dictionnary with an existing one
\param dict : dictionnary to copy from
***********************************************/
template<typename Type>
void ArithmeticCoder<Type>::setDictionnary(const ADictionnary<Type>& dict)
{
    m_dictionnary = dict;
}


/********************************************//**
\brief Remove duplicates entry in a QVector without sorting so each entry is unique and in a particular
order
\param input : input vector
\return vector unsorted with unique elements
***********************************************/
template<typename Type>
QVector<Type> ArithmeticCoder<Type>::removeDuplicatesFromVector(const QVector<Type>& input)
{
    std::vector<Type> vec = input.toStdVector();
    typename std::vector<Type>::iterator r, w;
    std::set<Type> tmpset ;
    
    for(r = vec.begin(), w = vec.begin(); r != vec.end(); ++r)
    {
        if(tmpset.insert(*r).second)
        {
            *w++ = *r ;
        }
    }
    // Erase duplicates
    vec.erase(w, vec.end());
    QVector<Type> result;
    return result.fromStdVector(vec);
}

#include <iostream>

/********************************************//**
\brief Fill dictionnary with input array
\param array : input array which contains all the data to encode later
\return arithmetic dictionnary
***********************************************/
template<typename Type>
bool ArithmeticCoder<Type>::fillDictionnary(const QList< QVector<Type> >& xArray,
                                            const QList< QVector<Type> >& yArray,
                                            const QList< QVector<Type> >& zArray,
                                            ADictionnary<Type>& dict,
                                            bool* stop)
{
    assert(xArray.size() == yArray.size() && xArray.size() == zArray.size());
    m_dictionnary.nbSymbolsEncoded = 0;
    
    for(int i = 0; i < xArray.size(); i++)
    {
        m_dictionnary.nbSymbolsEncoded += xArray[i].size() + yArray[i].size() + zArray[i].size();
        
        if(*stop)
        {
            return false;
        }
        // Remove non-unique elements to get symbols
        QVector<Type> symbols = xArray[i];
        symbols += yArray[i];
        symbols += zArray[i];
        symbols = removeDuplicatesFromVector(symbols);
        
        // Set cum_counts
        int cumcount = 0;
        for (int k = 0; k < symbols.size(); k++)
        {
            if(*stop)
            {
                return false;
            }
            int count = xArray[i].count(symbols[k]) +
                        yArray[i].count(symbols[k]) +
                        zArray[i].count(symbols[k]);
            cumcount += count;
   
            // Try to find if symbol already exist
            bool found = false;
            for(int f = 0; f < m_dictionnary.cum_counts.size(); f++)
            {
                found = false;
                if(*stop)
                {
                    return false;
                }
                // If found, update count only
                if(m_dictionnary.cum_counts[f].first == symbols[k])
                {
                    m_dictionnary.cum_counts[f].second += count;
                    found = true;
                    break;
                }
            }
            // If symbol not found, add symbol and count
            if(!found)
            {
                m_dictionnary.cum_counts.push_back(QPair<Type, int>(symbols[k], count));
            }
        }
    }
    
    // Adjust cumcounts
    int cumCount = 0;
    for (int i = 0; i < m_dictionnary.cum_counts.size(); i++)
    {
        if(*stop)
        {
            return false;
        }
        cumCount += m_dictionnary.cum_counts[i].second;
        m_dictionnary.cum_counts[i].second = cumCount;
    }
    
    // Set totalCount
    m_dictionnary.totalCount = cumCount;
    
    // Set m
    m_dictionnary.m = 0;
    while(pow(2.0, m_dictionnary.m) / 4 < m_dictionnary.totalCount)
    {
        m_dictionnary.m++;
    }
    
    dict = m_dictionnary;
    return true;
}


/********************************************//**
\brief Encode input array according to dictionnary to get 0 and 1 bytes
\param array : input array which contains all the data to encode
\return vector of bytes corresponding to vector array encoded
***********************************************/
template<typename Type>
bool ArithmeticCoder<Type>::encode(const QList< QVector<Type> >& xArray,
                                   const QList< QVector<Type> >& yArray,
                                   const QList< QVector<Type> >& zArray,
                                   std::string& encoded,
                                   bool* stop)
{
    assert(xArray.size() > 0 && yArray.size() > 0 && zArray.size());
    assert(xArray.size() == yArray.size() && xArray.size() == zArray.size());
    encoded = "";
    
    int scaleE3 = 0;
    QString l_bin;
    
    int l = 0; // lower bound of the interval
    int u = pow(2.0, m_dictionnary.m) - 1;  // upper bound of the interval
    
    for (int k = 0; k < xArray.size(); k++)
    {
        if(*stop)
        {
            return false;
        }
        if(!encodeSymbol(xArray[k], l, u, scaleE3, l_bin, encoded, stop))
        {
            return false;
        }

        if(!encodeSymbol(yArray[k], l, u, scaleE3, l_bin, encoded, stop))
        {
            return false;
        }
    
        if(!encodeSymbol(zArray[k], l, u, scaleE3, l_bin, encoded, stop))
        {
            return false;
        }
    }
    
    // Encode first byte of last status
    l_bin = int2binarystring(l, m_dictionnary.m);
    encoded.push_back((char)(l_bin[0] == '1'));
    
    // If scale of E3 is not null, then we need to add complement bits
    while(scaleE3 > 0)
    {
        encoded.push_back((char)!(l_bin[0] == '1'));
        scaleE3--;
    }
    // Add remaining bits
    for(int i = 1; i < l_bin.length(); i++)
    {
        encoded.push_back((char)(l_bin[i] == '1'));
    }
    
    // Pad if necessary
    std::div_t div;
    div = std::div(encoded.size(), m_dictionnary.m);
    int remainder = m_dictionnary.m - div.rem;
    
    bool invBit = !(bool)encoded[encoded.length() - 1];
    while(remainder > 0)
    {
        encoded.push_back(invBit);
        remainder--;
    }
    
    return true;
}
    

template<typename Type>
bool ArithmeticCoder<Type>::encodeSymbol(const QVector<Type>& array,
                                         int& l,
                                         int& u,
                                         int& scaleE3,
                                         QString& l_bin,
                                         std::string& encoded,
                                         bool* stop)
{
    QString u_bin;
    
    QChar l_msb;
    QChar u_msb;
    QChar l_2msb;
    QChar u_2msb;
    
    for (int k = 0; k < array.size(); k++)
    {
        if(*stop)
        {
            return false;
        }
        // Find cum_count of the current symbol
        typename QVector< QPair<Type, int> >::iterator it = find_if(m_dictionnary.cum_counts.begin(), m_dictionnary.cum_counts.end(), bind1st(pair_equal_to<Type, int>(), array[k]));
        int cumCount = it->second;
        int cumCountPrev = 0;
        
        it--;
        if(array[k] != m_dictionnary.cum_counts[0].first)
        {
            cumCountPrev = it->second;
        }
        it++;
        
        // update lower and upper bounds
        double u_l_one = (double)(u - l + 1);
        int lNew = l + floor(u_l_one * cumCountPrev / (double)m_dictionnary.totalCount);
        int uNew = l + floor(u_l_one * cumCount / (double)m_dictionnary.totalCount) - 1;
        l = lNew;
        u = uNew;
        
        // Convert bounds to binary code
        l_bin = int2binarystring(l, m_dictionnary.m);
        u_bin = int2binarystring(u, m_dictionnary.m);
        
        // Get most significant bit
        l_msb = l_bin[0];
        u_msb = u_bin[0];
        
        // Get second most significant bit
        l_2msb = l_bin[1];
        u_2msb = u_bin[1];
        
        // Check if bounds are in the same half or in the middle of the interval(E1, E2 or E3 mapping)
        while((l_msb == u_msb) || ((u_2msb == '0') && (l_2msb == '1')))
        {
            if(*stop)
            {
                return false;
            }
            // Check if bounds are in the same half of the interval(E1/E2 mapping)
            if(l_msb == u_msb)
            {
                // Remove most significant bit (shift to the left by 1 bit)
                l_bin = l_bin.mid(1, l_bin.size()-1);
                u_bin = u_bin.mid(1, u_bin.size()-1);
                
                // Store most significant bit
                encoded.push_back(l_msb == '1');
                
                // Shift bit
                l_bin.push_back('0');
                u_bin.push_back('1');
                
                while(scaleE3 > 0)
                {
                    encoded.push_back(!(l_msb == '1'));
                    scaleE3--;
                }
                
                // Update most significant bit
                l_msb = l_bin[0];
                u_msb = u_bin[0];
                
                // Update second most significant bit
                l_2msb = l_bin[1];
                u_2msb = u_bin[1];
            }
            // Check if bounds are in the middle of the interval (E3 mapping)
            if((u_2msb == '0') && (l_2msb == '1'))
            {
                // Remove most significant bit (shift to the left by 1 bit)
                l_bin = l_bin.mid(2, l_bin.size()-1);
                u_bin = u_bin.mid(2, u_bin.size()-1);
                
                // Shift bit
                l_bin = u_2msb + l_bin + '0';
                u_bin = l_2msb + u_bin + '1';
                
                // Update most significant bit
                l_msb = l_bin[0];
                u_msb = u_bin[0];
                
                // Update second most significant bit
                l_2msb = l_bin[1];
                u_2msb = u_bin[1];
                
                scaleE3++;
            }
        }
        // Convert bounds back to integers
        l = binary2int(l_bin);
        u = binary2int(u_bin);
    }
    return true;
}


/********************************************//**
\brief Decode input array according to dictionnary to retrieve symbols
\param array : input array which contains 0 and 1 bytes to be decoded
\return vector of symbols corresponding to vector array decoded
***********************************************/
template<typename Type>
bool ArithmeticCoder<Type>::decode(const std::string& array, QVector<Type>& decoded, bool* stop)
{
    decoded.clear();
    QString l_bin;
    QString u_bin;
    QString t_bin;
    
    int nbSymDecoded = 0;
    int t = 0;
    int l = 0;
    int u = pow(2.0, m_dictionnary.m) - 1;
    
    // Read first m bits
    for (int j = 0; j < m_dictionnary.m; j++)
    {
        t_bin += bool2char(array[j]);
    }
    
    bool print = false;
    unsigned int currIdxArray = m_dictionnary.m - 1;
    while(nbSymDecoded < m_dictionnary.nbSymbolsEncoded && currIdxArray < array.size())
    {
        if(*stop)
        {
            return false;
        }
        // Get tag
        t = binary2int(t_bin);
        double num = t - l + 1;
        num = (num * m_dictionnary.totalCount) - 1;
        double u_l_one = (double)(u - l + 1);
        int tag = floor(num / u_l_one);
        
        // Try to decode the symbol
        for(int c = 0; c < m_dictionnary.cum_counts.size(); c++)
        {
            if(*stop)
            {
                return false;
            }
            bool tagFound = false;
            int lNew, uNew;
            if(c == 0)
            {
                tagFound = (tag >= 0 && tag < m_dictionnary.cum_counts[c].second);
                lNew = l;
                uNew = l + floor(u_l_one * m_dictionnary.cum_counts[c].second / (double)m_dictionnary.totalCount) - 1;
            }
            else if(c > 0 && c < m_dictionnary.cum_counts.size() - 1)
            {
                tagFound = (tag >= m_dictionnary.cum_counts[c-1].second && tag < m_dictionnary.cum_counts[c].second);
                lNew = l + floor(u_l_one * m_dictionnary.cum_counts[c-1].second / (double)m_dictionnary.totalCount);
                uNew = l + floor(u_l_one * m_dictionnary.cum_counts[c].second / (double)m_dictionnary.totalCount) - 1;
            }
            else
            {
                tagFound = (tag >= m_dictionnary.cum_counts[c-1].second && tag <= m_dictionnary.cum_counts[c].second);
                lNew = l + floor(u_l_one * m_dictionnary.cum_counts[c-1].second / (double)m_dictionnary.totalCount);
                uNew = l + floor(u_l_one * m_dictionnary.cum_counts[c].second / (double)m_dictionnary.totalCount) - 1;
            }
            // If tag is found, decode symbol associated to it
            if(tagFound)
            {
                decoded.push_back(m_dictionnary.cum_counts[c].first);
                nbSymDecoded++;
                if(nbSymDecoded == m_dictionnary.nbSymbolsEncoded)
                {
                    return true;
                }
                // Update bounds
                l = lNew;
                u = uNew;
                
                if((float)m_dictionnary.cum_counts[c].first == 129.2f && l == 214687)
                {
                    print = true;
                }
                break;
            }
        }
        // Convert bounds to binary code
        l_bin = int2binarystring(l, m_dictionnary.m);
        u_bin = int2binarystring(u, m_dictionnary.m);
        
        // Get most significant bit
        QChar l_msb = l_bin[0];
        QChar u_msb = u_bin[0];
        
        // Get second most significant bit
        QChar l_2msb = l_bin[1];
        QChar u_2msb = u_bin[1];
        
        // Check if bounds are in the same half or in the middle of the interval(E1, E2 or E3 mapping)
        while((l_msb == u_msb) || ((u_2msb == '0') && (l_2msb == '1')))
        {
            if(*stop)
            {
                return false;
            }
            // Check if bounds are in the same half of the interval(E1/E2 mapping)
            if(l_msb == u_msb)
            {
                // Remove most significant bit (shift to the left by 1 bit)
                l_bin = l_bin.mid(1, l_bin.size()-1);
                u_bin = u_bin.mid(1, u_bin.size()-1);
                t_bin = t_bin.mid(1, t_bin.size()-1);
                
                // Shift bit
                l_bin.push_back('0');
                u_bin.push_back('1');
                currIdxArray++;
                t_bin.push_back(bool2char(array[currIdxArray]));
                
                // Update most significant bit
                l_msb = l_bin[0];
                u_msb = u_bin[0];
                
                // Update second most significant bit
                l_2msb = l_bin[1];
                u_2msb = u_bin[1];
            }
            // Check if bounds are in the middle of the interval (E3 mapping)
            if((u_2msb == '0') && (l_2msb == '1'))
            {
                char compMSB = '0';
                if( t_bin[1] == '0')
                {
                    compMSB = '1';
                }
                
                // Remove most significant bit (shift to the left by 1 bit)
                l_bin = l_bin.mid(2, l_bin.size()-1);
                u_bin = u_bin.mid(2, u_bin.size()-1);
                t_bin = t_bin.mid(2, t_bin.size()-1);
                
                // Shift bit
                l_bin = u_2msb + l_bin + '0';
                u_bin = l_2msb + u_bin + '1';
                
                currIdxArray++;
                t_bin = compMSB + t_bin + bool2char(array[currIdxArray]);
                
                // Update most significant bit
                l_msb = l_bin[0];
                u_msb = u_bin[0];
                
                // Update second most significant bit
                l_2msb = l_bin[1];
                u_2msb = u_bin[1];
            }
        }
        // Convert bounds back to integers
        l = binary2int(l_bin);
        u = binary2int(u_bin);
    }
    
    return true;
}

#endif