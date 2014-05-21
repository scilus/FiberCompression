#ifndef HUFFMAN_UTILS_IO
#define HUFFMAN_UTILS_IO

#include "Huffman.h"

#include <algorithm>
#include <assert.h>
#include <bitset>
#include <fstream>
#include <iostream>
#include <iterator>
#include <ostream>
#include <map>
#include <math.h>

#include <QString>
#include <QVector>


template<typename Type, typename DictionnaryType>
bool symbolListToIndexList(DictionnaryType& dict,
                           const QVector<Type> symbols,
                           QVector<int>& indexList,
                           bool* stop)
{
	indexList.clear();
    assert(symbols.size() > 0);
    
	for (int i = 0; i < (int)(symbols.size()); i++)
    {
        if(*stop)
        {
            return false;
        }
        int index = dict.symbols.indexOf(symbols[i]);

        if(index != -1)
        {
            indexList.push_back(index);
        }
	}
	return true;
}

#include "limits.h"

template<typename Type, typename DictionnaryType>
bool indexListToSymbolList(DictionnaryType& dict,
                           const QVector<int> list,
                           QVector<Type>& symbolList,
                           bool* stop)
{
	symbolList.clear();
	for (int i = 0; i < (int)(list.size()); i++)
    {
        if(*stop)
        {
            return false;
        }
		symbolList.push_back(dict.symbols[list[i]]);
        //std::cout << dict.symbols[list[i]] << "(" << list[i] << ") ";
	}
	return true;
}

#endif