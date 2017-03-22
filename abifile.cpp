/*
 * abi_file.cpp
 *
 * the class implementation to access the ABI trace files
 * Written by Conrad Shyu, July 30, 2004
 *
 * Initiative for Bioinformatics and Evolutionary Studies (IBEST)
 * Department of Bioinformatics and Computational Biology (BCB)
 * Department of Biological Sciences
 * University of Idhao, Moscow, ID 83844
 *
 * last updated on August 1, 2004
*/
#include <abitag.h>
#include <abifile.h>

/*
 * open the ABI trace file
*/
AbiFile::AbiFile(
    const char* _szFile )
{
    if ( !LoadFile( _szFile ) )
    {
#ifdef _DEBUG
        cout << "file failed to load" << endl;
#endif
        exit( 1 );
    }
}

AbiFile::AbiFile(
    string& _szFile )
{
    if ( !LoadFile( _szFile.c_str() ) )
    {
#ifdef _DEBUG
        cout << "file failed to load" << endl;
#endif
        exit( 1 );
    }
}

/*
 * load the entire tracefile into memory
*/
bool AbiFile::LoadFile(
    const char* _szFilename )
{
    struct stat fs;

    // try to get the status of the file
    if ( stat( _szFilename, &fs ) )
    {
        return( false );
    }

    szAbifBuffer = new unsigned char [ fs.st_size ];
    ifstream ifTraceFile( _szFilename, ios::in | ios::binary );

    if ( !ifTraceFile )
    {
        return( false );
    }

    ifTraceFile.read( reinterpret_cast<char*>( szAbifBuffer ), fs.st_size );
    ifTraceFile.close();

    // make sure the file contains the ABI signature "ABIF"
    if ( strncmp( reinterpret_cast<const char*>( szAbifBuffer ), "ABIF", 4 ) )
    {
        return( false );
    }

    // now parse the file, begin with the main tag list
    AbiTagRecord abiMainTag( szAbifBuffer, 6 );
    unsigned int entry = abiMainTag.GetDataValue();

    for ( int i = 0; i < abiMainTag.GetRecordCount(); ++i )
    {
        AbiTagRecord data( szAbifBuffer, entry );
        abiTagList.push_back( data );
        entry += abifTAGSIZE;
    }

#ifdef _DEBUG
    // print out the tag records
    for ( list<AbiTagRecord>::iterator tag = abiTagList.begin();
        !( tag == abiTagList.end() ); ++tag )
    {
        cout << ( *tag ).GetFlagName() << ": ";
        cout << ( *tag ).GetTypeName() << ", data value: ";
        cout << ( *tag ).GetDataValue() << ", record size: ";
        cout << ( *tag ).GetRecordSize() << ", record count: ";
        cout << ( *tag ).GetRecordCount() << endl;
    }
#endif

    return( true );
}

/*
 * export GeneScan analyzed data
*/
list<SIGNAL>& AbiFile::GetGSData(
    list<SIGNAL>& _data )
{
    string szCAPTION[] = { "Filter 1", "Filter 2", "Filter 3", "Filter 4" };
    list<AbiTagRecord>::iterator tag;
    SIGNAL stData;

    // loop through index 9, 10, 11, 12
    for ( int i = 0; i < 4; ++i )
    {
        tag = FindFlag( "DATA", ( i + 9 ) );

        if ( tag == abiTagList.end() )
        {
            continue;
        }

        stData.szCaption = szCAPTION[ i ];
        GetShort( tag, stData.vSignal );

#ifdef _DEBUG
        cout << stData.szCaption << ": " << stData.vSignal.size() << endl;
#endif

        _data.push_back( stData );
    }

#ifdef _DEBUG
    cout << "record: " << _data.size() << endl;
#endif

    return( _data );
}

/*
 * export the CCD raw data
*/
list<SIGNAL>& AbiFile::GetCCDData(
    list<SIGNAL>& _data )
{
    string szCAPTION[] = { "Filter 1", "Filter 2", "Filter 3", "Filter 4" };
    list<AbiTagRecord>::iterator tag;
    SIGNAL stData;

    // loop through index 1, 2, 3, 4
    for ( int i = 0; i < 4; ++i )
    {
        tag = FindFlag( "DATA", ( i + 1 ) );

        if ( tag == abiTagList.end() )
        {
            continue;
        }

        stData.szCaption = szCAPTION[ i ];
        stData.vSignal.clear();
        GetShort( tag, stData.vSignal );

#ifdef _DEBUG
        cout << stData.szCaption << ": " << stData.vSignal.size() << endl;
#endif

        _data.push_back( stData );
    }

#ifdef _DEBUG
    cout << "record: " << _data.size() << endl;
#endif

    return( _data );
}

/*
 * export electrophoresis status
*/
list<SIGNAL>& AbiFile::GetEPData(
    list<SIGNAL>& _data )
{
    string szCAPTION[] = { "Voltage", "mAmps", "Watts", "Temperature" };
    list<AbiTagRecord>::iterator tag;
    SIGNAL stData;

    // loop through index 5, 6, 7, 8
    for ( int i = 0; i < 4; ++i )
    {
        tag = FindFlag( "DATA", ( i + 5 ) );

        if ( tag == abiTagList.end() )
        {
            continue;
        }

        stData.szCaption = szCAPTION[ i ];
        stData.vSignal.clear();
        GetShort( tag, stData.vSignal );

#ifdef _DEBUG
        cout << stData.szCaption << ": " << stData.vSignal.size() << endl;
#endif

        _data.push_back( stData );
    }

#ifdef _DEBUG
    cout << "record: " << _data.size() << endl;
#endif

    return( _data );
}

/*
 * export peak record
*/
list<PEAK>& AbiFile::GetPeakData(
    list<PEAK>& _data )
{
    string szCAPTION[] = { "Filter 1", "Filter 2", "Filter 3", "Filter 4" };
    list<AbiTagRecord>::iterator tag;
    PEAK stPeak;
    int count;

    for ( int i = 0; i < 4; ++i )
    {
        // locate the peak record
        tag = FindFlag( "PK_#", ( i + 1 ) );

        if ( tag == abiTagList.end() )
        {
            continue;
        }

        // get the number of records
        count = GetShort( tag );

        if ( !( count > 0 ) )
        {
            continue;
        }

        tag = FindFlag( "PEAK", ( i + 1 ) );
        stPeak.szCaption = szCAPTION[ i ];
        stPeak.lpPeak.clear();
        GetPeakRecord( tag, stPeak.lpPeak, count );

#ifdef _DEBUG
        cout << stPeak.szCaption << ": " << stPeak.lpPeak.size() << endl;
#endif

        _data.push_back( stPeak );
    }

#ifdef _DEBUG
    cout << "record: " << _data.size() << endl;
#endif

    return( _data );
}   // end of GetPeakData()

/*
 * extract the peak data from the file
*/
list<PEAKDATA>& AbiFile::GetPeakRecord(
    list<AbiTagRecord>::iterator _tag,
    list<PEAKDATA>& _data,
    int _count )
{
    int entry = ( *_tag ).GetDataValue();
    PEAKDATA stPeakData;

    for ( int i = 0; i < _count; ++i )
    {
        stPeakData.nPoint = GetLong( entry );       entry += abiLONG;
        stPeakData.nHeight = GetShort( entry );     entry += abiSHORT;
        stPeakData.nBegin = GetLong( entry );       entry += abiLONG;
        stPeakData.nEnd = GetLong( entry );         entry += abiLONG;
        stPeakData.nBeginHi = GetShort( entry );    entry += abiSHORT;
        stPeakData.nEndHi = GetShort( entry );      entry += abiSHORT;
        stPeakData.nArea = GetLong( entry );        entry += abiLONG;
        stPeakData.nVolume = GetLong( entry );      entry += abiLONG;
        stPeakData.dSize = GetFloat( entry );       entry += abiFLOAT;
        stPeakData.bEdit = GetBool( entry );        entry += abiBOOL;
        GetString( entry, 64, stPeakData.szLabel ); entry += 64;
        _data.push_back( stPeakData );
    }

    return( _data );
}   // end of GetPeakRecord()

/*
 * find the tag record with a specified flag name and id
*/
list<AbiTagRecord>::iterator AbiFile::FindFlag(
   const string& _flag, const int _fid )
{
    list<AbiTagRecord>::iterator i;

    for ( i = abiTagList.begin(); !( i == abiTagList.end() ); ++i )
    {
        if ( ( ( *i ).GetFlagName() == _flag ) && ( ( *i ).GetFlagID() == _fid ) )
        {
            return( i );
        }
    }   // search the entire list for the given flag and file identification

    return( abiTagList.end() );
}   // end of FindFlag()

/*
 * get a character from the file
*/
char AbiFile::GetChar(
    int _entry )
{
    return( szAbifBuffer[ _entry ] );
}

char AbiFile::GetChar(
    list<AbiTagRecord>::iterator _i )
{
    return( static_cast<char>( ( ( *_i ).GetDataValue() >> 0x18 ) & 0xFF ) );
}

vector<char>& AbiFile::GetChar(
    list<AbiTagRecord>::iterator _i, vector<char>& _v )
{
    int entry = ( *_i ).GetDataValue();
    _v.clear();

    for ( int i = 0; i < ( *_i ).GetRecordCount(); ++i, ++entry )
    {
        _v.push_back( szAbifBuffer[ entry ] );
    }

    return( _v );
}

bool AbiFile::GetBool(
    int _entry )
{
    return( static_cast<bool>( GetShort( _entry ) ) );
}

/*
 * get a two-byte integer from the file
*/
int AbiFile::GetShort(
    int _entry )
{
    int value;

    value  = szAbifBuffer[ _entry++ ] << 0x8;
    value += szAbifBuffer[ _entry ];

    return( value );
}

int AbiFile::GetShort(
    list<AbiTagRecord>::iterator _i )
{
    // only return the first two bytes
    return( ( *_i ).GetDataValue() >> 0x10 );
}

/*
 * extract all the data from the file
*/
vector<int>& AbiFile::GetShort(
    list<AbiTagRecord>::iterator _i, vector<int>& _v )
{
    int entry = ( *_i ).GetDataValue();
    int value;
    _v.clear();     // clear all elements

    for ( int i = 0; i < ( *_i ).GetRecordCount(); ++i )
    {
        value  = szAbifBuffer[ entry++ ] << 0x8;
        value += szAbifBuffer[ entry++ ];
        _v.push_back( value );
    }

    return( _v );
}

/*
 * extaact a long integer from the file
*/
int AbiFile::GetLong(
    int _entry )
{
    int value;

    value  = szAbifBuffer[ _entry++ ] << 0x18;
    value += szAbifBuffer[ _entry++ ] << 0x10;
    value += szAbifBuffer[ _entry++ ] << 0x8;
    value += szAbifBuffer[ _entry ];

    return( value );
}

int AbiFile::GetLong(
    list<AbiTagRecord>::iterator _i )
{
    return( ( *_i ).GetDataValue() );
}

vector<int>& AbiFile::GetLong(
    list<AbiTagRecord>::iterator _i, vector<int>& _v )
{
    int entry = ( *_i ).GetDataValue();
    int value;
    _v.clear();     // clear all elements

    for ( int i = 0; i < ( *_i ).GetRecordCount(); ++i )
    {
        value  = szAbifBuffer[ entry++ ] << 0x18;
        value += szAbifBuffer[ entry++ ] << 0x10;
        value += szAbifBuffer[ entry++ ] << 0x8;
        value += szAbifBuffer[ entry++ ];
        _v.push_back( value );
    }

    return( _v );
}

double AbiFile::ReadFloat(
    unsigned int _data )
{
/*
    int bias = ( ( _data >> 0x17 ) & 0xFF ) - 127;
    double sign = ( ( _data >> 0x1F ) & 0x1 ) ? -1.0 : 1.0;
    double mantissa = ( _data & 0x7FFFFF ) / 8388607.0 + 1.0;

    return( sign * mantissa * pow( 2.0, bias ) );
*/
    union INT2FLOAT
    {
        unsigned int i;     // share the 32-bit space
        float f;
    };

    INT2FLOAT value; value.i = _data;

    return( static_cast<double>( value.f ) );
}

double AbiFile::GetFloat(
    int _entry )
{
    unsigned int value;

    value  = szAbifBuffer[ _entry++ ] << 0x18;
    value += szAbifBuffer[ _entry++ ] << 0x10;
    value += szAbifBuffer[ _entry++ ] << 0x8;
    value += szAbifBuffer[ _entry ];

    return( ReadFloat( value ) );
}

double AbiFile::GetFloat(
    list<AbiTagRecord>::iterator _i )
{
    return( ReadFloat( static_cast<unsigned int>( ( *_i ).GetDataValue() ) ) );
}

vector<double>& AbiFile::GetFloat(
    list<AbiTagRecord>::iterator _i, vector<double>& _v )
{
    int entry = ( *_i ).GetDataValue();
    unsigned int value;

    for ( int i = 0; i < ( *_i ).GetRecordCount(); ++i )
    {
        value  = szAbifBuffer[ entry++ ] << 0x18;
        value += szAbifBuffer[ entry++ ] << 0x10;
        value += szAbifBuffer[ entry++ ] << 0x8;
        value += szAbifBuffer[ entry++ ];
        _v.push_back( ReadFloat( value ) );
    }

    return( _v );
}

/*
 * get the date information from the file in string format
*/
string& AbiFile::GetDate(
    list<AbiTagRecord>::iterator _i, string& _s )
{
    int value = ( *_i ).GetDataValue();
    int year = ( value >> 0x10 ) & 0xFFFF;  // byte[1][2]: year
    int month = ( value >> 0x8 ) & 0xFF;    // byte[3]: month
    int day = value & 0xFF;         // byte[4]: day
    char sz_date[ sizeDATE ];

    sprintf( sz_date, "%02d/%02d/%4d", month, day, year );
    _s.assign( sz_date );

    return( _s );
}

/*
 * time: byte[1][2][3][4]=hh:mm:ss:tt
*/
double AbiFile::GetTime(
    list<AbiTagRecord>::iterator _i )
{
    int value = ( *_i ).GetDataValue();
    double tt = ( value & 0xFF ) / 1000.0;          // byte[4]: one thousandth
    int ss = ( ( value >> 0x8 ) & 0xFF );           // byte[3]: seconds
    int mm = ( ( value >> 0x10 ) & 0xFF ) * 60;     // byte[2]: minutes
    int hh = ( ( value >> 0x18 ) & 0xFF ) * 3600;   // byte[1]: hours
    tt += ( ss + mm + hh );

    return( tt );
}

/*
 * get the time information from the file in string format
*/
string& AbiFile::GetTime(
    list<AbiTagRecord>::iterator _i, string& _s )
{
    int value = ( *_i ).GetDataValue();
    int tt = ( value & 0xFF );              // byte[0]: one thousandth
    int ss = ( ( value >> 0x8 ) & 0xFF );   // byte[3]: seconds
    int mm = ( ( value >> 0x10 ) & 0xFF );  // byte[2]: minutes
    int hh = ( ( value >> 0x18 ) & 0xFF );  // byte[1]: hours

    char sz_time[ sizeTIME ];

    sprintf( sz_time, "%2d:%2d:%2d.%2d", hh, mm, ss, tt );
    _s.assign( sz_time );

    return( _s );
}

/*
 * get a string of given size from the file
*/
string& AbiFile::GetString(
    int _entry, int _size, string& _s )
{
    _s.resize( _size );

    for ( int i = 0; i < _size; ++i, ++_entry )
    {
        _s[ 0 ] = szAbifBuffer[ _entry ];
    }

    return( _s );
}

string& AbiFile::GetString(
    list<AbiTagRecord>::iterator _i, string& _s )
{
   int entry = ( *_i ).GetDataValue();
    int length = ( *_i ).GetRecordSize();
    _s.resize( length );

    if ( length > 4 )
    {
        for ( int i = 0; i < length; ++i )
        {
            _s[ i ] = szAbifBuffer[ ++entry ];  // first byte is the length
        }
    }
    else
    {
        _s[ 0 ] = static_cast<char>( ( entry >> 0x18 ) & 0xFF );
        _s[ 1 ] = static_cast<char>( ( entry >> 0x10 ) & 0xFF );
        _s[ 2 ] = static_cast<char>( ( entry >> 0x8 ) & 0xFF );
        _s[ 3 ] = static_cast<char>( entry & 0xFF );
    }

    return( _s );
}

/*
 * test driver program
*/
/*
int main( int argc, char** argv )
{
    if ( argc < 2 )
    {
        cout << "usage: " << argv[ 0 ] << " abi_tracefile" << endl;
        return( 1 );
    }

    AbiFile q( argv[ 1 ] );
    list<SIGNAL> signal;
    list<PEAK> peak;

    q.GetGSData( signal );
    q.GetCCDData( signal );
    q.GetEPData( signal );
    q.GetPeakData( peak );
}
*/
