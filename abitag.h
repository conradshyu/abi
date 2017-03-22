/*
 * abitag.h
 *
 * ABI tracefile data structure
 *
 * Initiative for Bioinformatics and Evolutionary Studies (IBEST)
 * Department of Bioinformatics and Computational Biology (BCB)
 * Department of Biological Sciences
 * University of Idhao, Moscow, ID 83844
 *
 * last updated on August 1, 2004
*/
#ifndef _ABI_TAG_H
#define _ABI_TAG_H

// C++ header files
#include <vector>
#include <string>
#include <cstdlib>
#include <iostream>

using namespace std;

/*
 * all ABI FLAG records are 28 byes in length and exhibit the following structure:
*/
class AbiTagRecord
{
public:
    AbiTagRecord( unsigned char*, const unsigned int );
    ~AbiTagRecord() {}

    const string& GetFlagName() const   { return( szFlagName ); }
    const string& GetTypeName() const   { return( szTypeName ); }

    int GetFlagID() const       { return( nFlagID ); }
    int GetDataType() const     { return( nDataType ); }
    int GetRecordSize() const   { return( nRecordSize ); }
    int GetRecordCount() const  { return( nRecordCount ); }
    int GetRecordLength() const { return( nRecordLength ); }
    int GetDataValue() const    { return( nDataValue ); }
    int GetDataPadding() const  { return( nDataPadding ); }

private:
    string  szTypeName;     // ABI designated data type name
    string  szFlagName;     // ASCII flag name
    int nFlagID;            // id of the flag
    int nDataType;          // ABI designated data type
    int nRecordSize;        // length of referenced data records (bytes)
    int nRecordCount;       // number of referenced data records
    int nRecordLength;      // length of referenced data array (bytes)
    int nDataValue;         // either (1) the data itself or (2) a pointer
    int nDataPadding;       // purpose is unknown

    unsigned char*  szBuffer;
    unsigned int    nEntry;

    string& GetFlag( string& );
    int GetShort();
    int GetLong();
    int GetDataType( int );
};

#endif  // _ABI_TAG_H
