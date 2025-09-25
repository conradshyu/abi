/*
 * abi2csv.cpp
 *
 * convert ABI and AB1 files into CSV format
 *
 * written by Conrad Shyu, November 8, 2004
 *
 * Initiative for Bioinformatics and Evolutionary Studies (IBEST)
 * Program of Bioinformatics and Computational Biology (BCB)
 * Department of Biological Sciences
 * University of Idaho, Moscow, ID 83844
*/

// for standard c libraries
#include <errno.h>
#include <dirent.h>
#include <unistd.h>
#include <fnmatch.h>

#include <abitag.h>
#include <abifile.h>

// for c++ standard template library
#include <list>
#include <string>
#include <fstream>
#include <iostream>

using namespace std;

//#define _DEBUG

/*
 * platform dependent implementation
*/
#ifdef _WIN32
    #include <windows.h>
#else
    #include <dirent.h>
    #include <sys/types.h>
#endif

const unsigned int BUFFER_SIZE = 4192;

/*
 * get the list of filenames that match the pattern
*/
bool GetFilename(
    list<string>& _list,        // list of filename in given directory
    const string& _pattern )    // filename patterns
{
#ifdef _WIN32
    // for Win32 implementation
    WIN32_FIND_DATA stFindData;
    HANDLE          hFindFile;

    // find the first file
    hFindFile = FindFirstFile( _pattern, &stFindData );

    // search all files that match the pattern
    do {
        // make sure it is not a directory
        if ( !( stFindData.cFileName[ 0 ] == '.' ) )
        {
            _list.push_back( stFindData.cFileName );
        }

        FindNextFile( hFindFile, &stFindData );     // search for all matched files
    } while ( !( GetLastError() == ERROR_NO_MORE_FILES ) );

    FindClose( hFindFile );

#else
    // for other implementations
    struct dirent *st_dir;
    char buffer[ BUFFER_SIZE ]; string file;
    DIR *p_dir = opendir( "." );

    if ( !p_dir )
    {
        cout << "directory cannot be opened" << endl; return( false );
    }   // make sure the directory can be opened successfully

    while ( st_dir = readdir( p_dir ) )
    {
        if ( st_dir->d_name[ 0 ] == '.' )
        {
            continue;
        }   // skip the root directories and hidden files

        if ( st_dir->d_type == DT_DIR )
        {
            chdir( st_dir->d_name );    // change to the next directory
            GetFilename( _list, _pattern );
            chdir( ".." );              // revert back to the previous directory
            continue;
        }   // recursively search each directory

        if ( fnmatch( _pattern.c_str(), st_dir->d_name, FNM_FILE_NAME | FNM_PERIOD ) == FNM_NOMATCH )
        {
            continue;
        }   // search for files with specific patterns

        // put together the path and filename
        file.assign( getcwd( buffer, BUFFER_SIZE ) );
        file.append( "/" ); file.append( st_dir->d_name );

        _list.push_back( file );
    }   // scan through the entire directory

    closedir( p_dir );
#endif

    return( true );
}   // end of GetFilename()

bool WriteCSV(
    string& _filename,
    list<SIGNAL>& _signal )
{
    ofstream csv( _filename.c_str(), ios::out | ios::trunc );

    if ( !csv )
    {
        return( false );
    }

    list<SIGNAL>::iterator filter;

    // first write all the captions
    filter = _signal.begin();
    csv << "\"" << ( *filter ).szCaption.c_str() << "\"";

    for ( ++filter; !( filter == _signal.end() ); ++filter )
    {
        csv << ",\"" << ( *filter ).szCaption.c_str() << "\"";
    }

    csv << endl;

    // now write all the data
    for ( unsigned int i = 0; i < ( _signal.front() ).vSignal.size(); ++i )
    {
        filter = _signal.begin();
        csv << static_cast<int>( ( *filter ).vSignal[ i ] );

        for ( ++filter; !( filter == _signal.end() ); ++filter )
        {
            csv << "," << static_cast<int>( ( *filter ).vSignal[ i ] );
        }

        csv << endl;
    }

    csv.close(); return( true );
}

bool WriteCSV(
    string& _filename,
    list<PEAK>& _peak )
{
    ofstream csv( _filename.c_str(), ios::out | ios::trunc );

    if ( !csv )
    {
        return( false );
    }

    list<PEAK>::iterator peak;

    for ( peak = _peak.begin(); !( peak == _peak.end() ); ++peak )
    {
        // first write the filter name
        csv << "\"" << ( *peak ).szCaption.c_str() << "\"" << endl;
        csv << "\"Position\",\"Height\",\"BeginPeak\",\"EndPeak\",";
        csv << "\"BeginHeight\",\"EndHeight\",\"Area\",\"Size\"" << endl;

        list<PEAKDATA>::iterator p;

        // now, write all the data
        for ( p = ( *peak ).lpPeak.begin(); !( p == ( *peak ).lpPeak.end() ); ++p )
        {
            csv << ( *p ).nPoint << "," << ( *p ).nHeight << ",";
            csv << ( *p ).nBegin << "," << ( *p ).nEnd << ",";
            csv << ( *p ).nBeginHi << "," << ( *p ).nEndHi << ",";
            csv << ( *p ).nArea << "," << ( *p ).dSize << endl;
        }
    }

    csv.close(); return( true );
}

/*
 * main procedure
*/
int main( int argc, char** argv )
{
    // make sure we have enough parameters
    if ( argc < 2 )
    {
        cout << "usage: " << argv[ 0 ] << " extension" << endl;
        cout << "convert the ABI and AB1 files into CSV format" << endl;
        exit( 1 );
    }

    list<string> lpFile; lpFile.clear();
    list<string>::iterator i;
    string szFilename;
    string szExtension( "*." ); szExtension.append( argv[ 1 ] );

    GetFilename( lpFile, szExtension ); lpFile.sort();

    list<SIGNAL> signal; list<PEAK> peak;

#ifdef _DEBUG
    cout << "number of file(s): " << lpFile.size() << endl;
#endif

    for ( i = lpFile.begin(); !( i == lpFile.end() ); ++i )
    {
        szFilename = ( *i );

#ifdef _DEBUG
        cout << "filename: " << ( *i ) << endl;
#endif

        cout << "processing file " << szFilename << "...";
        AbiFile abi( szFilename );

        szFilename.resize( szFilename.length() - 4 );
        szFilename.append( "_raw.csv" );
        signal.clear();
        abi.GetGSData( signal ); abi.GetCCDData( signal );

        if ( !WriteCSV( szFilename, signal ) )
        {
            cout << "file writing error" << endl;
        }

        szFilename = ( *i );
        szFilename.resize( szFilename.length() - 4 );
        szFilename.append( "_peak.csv" );
        peak.clear();
        abi.GetPeakData( peak );
        WriteCSV( szFilename, peak );

        cout << " done" << endl;
     }

    cout << "all tasks are completed!" << endl; return( 0 );
}
