/*  Copyright (c) 2015 HPAC
 *  Written by 	Alexis Breust <alexis.breust@gmail.fr>
 * 		Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 */

#include<string>

std::string removeExtension( std::string const& filename ) {
    auto pivot = std::find( filename.rbegin(), filename.rend(), '.' );
    return pivot == filename.rend() ? filename : std::string( filename.begin(), pivot.base() - 1 );
}
struct MatchPathSeparator
{
    bool operator()( char ch ) const
    {
        return ch == '\\' || ch == '/';
    }
};
std::string basename( std::string const& pathname ) {
    std::string noext = removeExtension(pathname);
    return std::string(  
        std::find_if( pathname.rbegin(), pathname.rend(),
                      MatchPathSeparator() ).base(),
        pathname.end() );
}

