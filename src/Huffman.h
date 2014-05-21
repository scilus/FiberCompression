/* Copyright (c) 2012 the authors listed at the following URL, and/or
 the authors of referenced articles or incorporated external code:
 http://en.literateprograms.org/Huffman_coding_(C_Plus_Plus)?action=history&offset=20090129100015
 
 Permission is hereby granted, free of charge, to any person obtaining
 a copy of this software and associated documentation files (the
 "Software"), to deal in the Software without restriction, including
 without limitation the rights to use, copy, modify, merge, publish,
 distribute, sublicense, and/or sell copies of the Software, and to
 permit persons to whom the Software is furnished to do so, subject to
 the following conditions:
 
 The above copyright notice and this permission notice shall be
 included in all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 
 Retrieved from: http://en.literateprograms.org/Huffman_coding_(C_Plus_Plus)?oldid=16057
 and modified to suit our needs
 */

#ifndef HUFFMAN_H
#define HUFFMAN_H

#include <algorithm>
#include <bitset>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <math.h>
#include <ostream>
#include <queue>

#include <QSet>
#include <QString>
#include <QVector>


struct ColorDictionnary
{
	QVector<unsigned char> symbols;
    std::map<int, float> frequencies;
    
    void cleanUp()
    {
        symbols.clear();
        frequencies.clear();
    }
};

struct Dictionnary
{
    QVector<float> symbols;
    std::map<int, float> frequencies;
    
    void cleanUp()
    {
        symbols.clear();
        frequencies.clear();
    }
};


template<typename DataType, typename Frequency> class Hufftree
{
public:
    template<typename InputIterator>
    Hufftree(InputIterator begin, InputIterator end, bool* stop);
    
    ~Hufftree() { delete tree; }
    
    template<typename InputIterator>
    bool encode(InputIterator begin, InputIterator end, std::string& encoded, bool* stop);
    
    QVector<bool> encode(DataType const& value)
    {
        return encoding[value];
    }
    
    template<typename OutputIterator>
    bool decode(std::string const& v, OutputIterator iter, bool* stop);
    
    bool isStopped()
    {
        return stopped;
    }
    
private:
    struct Node;
    Node* tree;
    
    typedef std::map<DataType, QVector<bool> > encodemap;
    encodemap encoding;
    
    struct NodeOrder;
    bool stopped;
};


template<typename DataType, typename Frequency>
struct Hufftree<DataType, Frequency>::Node
{
    Frequency frequency;
    Node* leftChild;
    union
    {
        Node* rightChild; // if leftChild != 0
        DataType* data;   // if leftChild == 0
    };
    
    Node(Frequency f, DataType d):
    frequency(f),
    leftChild(0),
    data(new DataType(d))
    {}
    
    Node(Node* left, Node* right):
    frequency(left->frequency + right->frequency),
    leftChild(left),
    rightChild(right)
    {}
    
    ~Node()
    {
        if (leftChild)
        {
            delete leftChild;
            delete rightChild;
        }
        else
        {
            delete data;
        }
    }
    
    void fill(std::map<DataType, QVector<bool> >& encoding, QVector<bool>& prefix)
    {
        if (leftChild)
        {
            prefix.push_back(0);
            leftChild->fill(encoding, prefix);
            prefix.back() = 1;
            rightChild->fill(encoding, prefix);
            prefix.pop_back();
        }
        else
        {
            encoding[*data] = prefix;
        }
    }
};


template<typename DataType, typename Frequency>
template<typename InputIterator> Hufftree<DataType, Frequency>::Hufftree(InputIterator begin, InputIterator end, bool* stop):
tree(0),
stopped(false)
{
    std::priority_queue<Node*, QVector<Node*>, NodeOrder> pqueue;
    
    while (begin != end)
    {
        if(*stop)
        {
            stopped = true;
            return;
        }
        Node* dataNode = new Node(begin->second, begin->first);
        pqueue.push(dataNode);
        ++begin;
    }
    while (!pqueue.empty())
    {
        if(*stop)
        {
            stopped = true;
            return;
        }
        Node* top = pqueue.top();
        pqueue.pop();
        if (pqueue.empty())
        {
            tree = top;
        }
        else
        {
            Node* top2 = pqueue.top();
            pqueue.pop();
            pqueue.push(new Node(top, top2));
        }
    }
    QVector<bool> bitvec;
    tree->fill(encoding, bitvec);
}

template<typename DataType, typename Frequency>
struct Hufftree<DataType, Frequency>::NodeOrder
{
    bool operator()(Node* a, Node* b)
    {
        if (b->frequency < a->frequency)
            return true;
        if (a->frequency < b->frequency)
            return false;
        
        if (!a->leftChild && b->leftChild)
            return true;
        if (a->leftChild && !b->leftChild)
            return false;
        
        if (a->leftChild && b->leftChild)
        {
            if ((*this)(a->leftChild, b->leftChild))
                return true;
            if ((*this)(b->leftChild, a->leftChild))
                return false;
            return (*this)(a->rightChild, b->rightChild);
        }
        return *(a->data) < *(b->data);
    }
};

template<typename DataType, typename Frequency>
template<typename InputIterator>
bool Hufftree<DataType, Frequency>::encode(InputIterator begin,
                                           InputIterator end,
                                           std::string& encoded,
                                           bool* stop)
{
    encoded.clear();
    while (begin != end)
    {
        if(*stop)
        {
            return false;
        }
        typename encodemap::iterator i = encoding.find(*begin);
        for (int j = 0; j < i->second.size(); j++)
            //QVector<bool>::const_iterator it = i->second.begin(); it != i->second.end(); ++it)
        {
            //encoded.push_back((char)*it);
            encoded += i->second[j] ? 1 : 0;
        }
        ++begin;
    }
    return true;
}

template<typename DataType, typename Frequency>
template<typename OutputIterator>
bool Hufftree<DataType, Frequency>::decode(std::string const& v,
                                           OutputIterator iter,
                                           bool* stop)
{
    Node* node = tree;
    for (std::string::const_iterator i = v.begin(); i != v.end(); ++i)
    {
        if(*stop)
        {
            return false;
        }
        node = *i? node->rightChild : node->leftChild;
        if (!node -> leftChild)
        {
            *iter++ = *(node->data);
            node = tree;
        }
    }
    return true;
}


#endif