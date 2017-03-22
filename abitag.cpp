/*
 * abitag.cpp
 *
 * Written by Conrad Shyu, July 31, 2004
 * Initiative for Bioinformatics and Evolutionary Studies (IBEST)
 * Department of Bioinformatics and Computational Biology (BCB)
 * Department of Biological Sciences
 * University of Idaho, Moscow, ID 83844
*/
#include <abitag.h>

// ABI designated data type
string abiTYPENAME[ 26 ] =
{
    "IllegalType", "Byte", "Char", "Word", "Short", "Long", "Rational", "Float", "Double", "BCD",
    "Date", "Time", "Thumb", "Boolean", "Point", "Rect", "VPoint", "VRect", "PString", "CString",
    "Tag", "DeltaLZWcompression", "LZWcompression", "Directory", "UserType", "CustomUserType"
};

/*
 * all ABI processed tracefiles (chromatograms) begin with a single 128 byte
 * header:
 *
 * 41 42 49 46 00 65 74 64 69 72 00 00 00 01 03 FF
 * 00 1C 00 00 00 3C 00 00 07 00 00 00 0C 5B 00 00
 * 00 00 FF FF FF FF FF FF FF FF FF FF FF FF FF FF
 * FF FF FF FF FF FF FF FF FF FF FF FF FF FF FF FF
 * FF FF FF FF FF FF FF FF FF FF FF FF FF FF FF FF
 * FF FF FF FF FF FF FF FF FF FF FF FF FF FF FF FF
 * FF FF FF FF FF FF FF FF FF FF FF FF FF FF FF FF
 * FF FF FF FF FF FF FF FF FF FF FF FF FF FF FF FF
 *
 * contained within this header are four bytes (41 42 49 46) that designate the
 * file as an ABI processed tracefile (ABIF); two bytes of an unknown purpose
 * (00 65) and a single 28 byte data record (74 .. 6C) offset 6 bytes from the
 * beginning of the tracefile. This statically located record provides information
 * about the structure and location of the dynamically located multiple FLAG
 * region. The FLAG records function as the data-organizing unit of each ABI
 * processed tracefile. The reminder of the 128 byte header is padded with FF
 * bytes.
*/
AbiTagRecord::AbiTagRecord(
    unsigned char* _s, const unsigned int _i ) :
    szBuffer( _s ), nEntry( _i )
{
    // tag record can't start from 0
    if ( !( nEntry > 0 ) )
    {
        exit( 1 );
    }

    // parse the record; note: the order is very important!
    GetFlag( szFlagName );      // ascii flag name; 4 bytes
    nFlagID = GetLong();        // id of the flag; 4 bytes
    nDataType = GetDataType( GetShort() );  // ABI designated data type; 2 bytes
    nRecordSize = GetShort();   // size of referenced record (bytes); 2 bytes
    nRecordCount = GetLong();   // number of referenced data records; 4 bytes
    nRecordLength = GetLong();  // length of referenced data array (bytes); 4 bytes
    nDataValue = GetLong();     // either the data itself or a pointer; 4 bytes
    nDataPadding = GetLong();   // purpose is unknown; 4 bytes
}

/*
 * get a string from the binary tracefile
*/
string& AbiTagRecord::GetFlag(
    string& _s )
{
    _s.resize( 4 );

    _s[ 0 ] = szBuffer[ nEntry++ ];
    _s[ 1 ] = szBuffer[ nEntry++ ];
    _s[ 2 ] = szBuffer[ nEntry++ ];
    _s[ 3 ] = szBuffer[ nEntry++ ];

    return( _s );
}

/*
 * get an integer from the binary tracefile
*/
int AbiTagRecord::GetLong()
{
    int value;

    value  = szBuffer[ nEntry++ ] << 0x18;
    value += szBuffer[ nEntry++ ] << 0x10;
    value += szBuffer[ nEntry++ ] << 0x8;
    value += szBuffer[ nEntry++ ];

    return( value );
}

int AbiTagRecord::GetShort()
{
    int value;

   value  = szBuffer[ nEntry++ ] << 0x8;
    value += szBuffer[ nEntry++ ];

    return( value );
}

/*
 * translate the ABI designated data type
*/
int AbiTagRecord::GetDataType(
    int _type )
{
    switch ( _type )
    {
    case 0:     _type = 0; break;   // illegal type
    case 1:     _type = 1; break;   // byte
    case 2:     _type = 2; break;   // char
    case 3:     _type = 3; break;   // word
    case 4:     _type = 4; break;   // short
    case 5:     _type = 5; break;   // long
    case 6:     _type = 6; break;   // rational
    case 7:     _type = 7; break;   // float
    case 8:     _type = 8; break;   // double
    case 9:     _type = 9; break;   // BCD
    case 10:    _type = 10; break;  // date
    case 11:    _type = 11; break;  // time
    case 12:    _type = 12; break;  // thumb
    case 13:    _type = 13; break;  // boolean
    case 14:    _type = 14; break;  // point
    case 15:    _type = 15; break;  // rect
    case 16:    _type = 16; break;  // vpoint
    case 17:    _type = 17; break;  // vrect
    case 18:    _type = 18; break;  // pstring
    case 19:    _type = 19; break;  // cstring
    case 20:    _type = 20; break;  // tag
    case 128:   _type = 21; break;  // delta lzw compression
    case 256:   _type = 22; break;  // lzw compression
    case 1023:  _type = 23; break;  // directory
    case 1024:  _type = 24; break;  // user type
    default:    _type = 25;         // custom user type
    }

    szTypeName = abiTYPENAME[ _type ];
    return( _type );
}
