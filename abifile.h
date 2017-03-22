/*
 * abi_file.h
 *
 * The header file for the class implementation to access the ABI tracefile
 * Written by Conrad Shyu, July 30, 2004
 *
 * Initiative for Bioinformatics and Evolutionary Studies (IBEST)
 * Department of Bioinformatics and Computational Biology (BCB)
 * Department of Biological Sciences
 * University of Idhao, Moscow, ID 83844
 *
 * last updated on August 1, 2004
*/
#ifndef _ABI_FILE_H
#define _ABI_FILE_H

#include <sys/stat.h>

// C++ header files
#include <list>
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <cstring>
#include <iostream>

//#define _DEBUG

using namespace std;

const int abifTAGSIZE   = 28;
const int sizeFLOAT     = sizeof( float );
const int sizeDOUBLE    = sizeof( double );
const int sizeLDOUBLE   = sizeof( long double );
const int sizeDATE      = 20;       // format: mm/dd/yyyy
const int sizeTIME      = 20;       // format: hh:ss:mm.tt
const int abiSHORT      = 2;
const int abiLONG       = 4;
const int abiFLOAT      = 4;
const int abiBOOL       = 2;

// a total of 96 bytes
struct PEAKDATA
{
    int nPoint;     // data point; 4 bytes
    int nHeight;    // peak height; 2 bytes
    int nBegin;     // peak begin position; 4 bytes
    int nEnd;       // peak end position; 4 bytes
    int nBeginHi;   // peak begin height; 2 bytes
    int nEndHi;     // peak end height; 2 bytes
    int nArea;      // peak area; 4 bytes
    int nVolume;    // peak volume; 4 bytes
    double dSize;   // fragment size of peak (basepairs); 4 bytes
    bool bEdit;     // has this peak been edited; 2 bytes
    string szLabel; // peak label; 64 bytes
};

struct SIGNAL
{
    string szCaption;       // signal caption
    vector<int> vSignal;    // relative fluorescent intensity
};

struct PEAK
{
    string szCaption;
    list<PEAKDATA> lpPeak;
};

/*
 * class implementation to access the ABI tracefile
*/
class AbiFile
{
public:
    AbiFile();
    AbiFile( const char* );
    AbiFile( string&  );
    ~AbiFile()  { delete szAbifBuffer; }

    bool LoadFile( const char* );
    list<SIGNAL>&   GetCCDData( list<SIGNAL>& );
    list<SIGNAL>&   GetGSData( list<SIGNAL>& );
    list<SIGNAL>&   GetEPData( list<SIGNAL>& );
    list<PEAK>&     GetPeakData( list<PEAK>& );
    list<AbiTagRecord>& GetTagRecord( list<AbiTagRecord>& ) const;

private:
    list<AbiTagRecord>  abiTagList;
    unsigned char*      szAbifBuffer;

    bool    GetBool( int );
    int     GetShort( int );
    int     GetLong( int );
    char    GetChar( int );
    double  GetFloat( int );
    double  ReadFloat( unsigned int );
    int     GetShort( list<AbiTagRecord>::iterator );
    int     GetLong( list<AbiTagRecord>::iterator );
    char    GetChar( list<AbiTagRecord>::iterator );
    double  GetFloat( list<AbiTagRecord>::iterator );
    vector<int>&    GetShort( list<AbiTagRecord>::iterator, vector<int>& );
    vector<int>&    GetLong( list<AbiTagRecord>::iterator, vector<int>& );
    vector<char>&   GetChar( list<AbiTagRecord>::iterator, vector<char>& );
    vector<double>& GetFloat( list<AbiTagRecord>::iterator, vector<double>& );

    double  GetTime( list<AbiTagRecord>::iterator );
    string& GetTime( list<AbiTagRecord>::iterator, string& );
    string& GetDate( list<AbiTagRecord>::iterator, string& );
    string& GetString( list<AbiTagRecord>::iterator, string& );
    string& GetString( int, int, string& );
    list<AbiTagRecord>::iterator FindFlag( const string&, const int );
    list<PEAKDATA>& GetPeakRecord( list<AbiTagRecord>::iterator, list<PEAKDATA>&, int );
};

#endif  // _ABI_FILE_H
