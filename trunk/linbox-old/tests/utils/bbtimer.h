/* File: src/examples/utils/bbtimer.h
 * Author: William J Turner for LinBox group
 */

#ifndef _BBTIMER_
#define _BBTIMER_

/* class: BBTimer
 * function: GetCPUTime
 *  
 * C++ class and function to provide benchmarking resource usage.
 * GetCPUTime is (slight) modification of function by Angel Diaz.
 * 
 * Based on bbtime.h by Angel Diaz.
 * 
 * Converted into class to allow for multithread usage by 
 * eliminating the static variables.
 */
   
#include <iostream>
#include <string>
#include <algorithm>
#include <iterator>
#include <ctime>

#include <sys/times.h>
#include <sys/param.h>

extern "C" int gethostname(char *name, int namelen);

/** BlackBox Timer.
  * Provides benchmarking resource usage.  Based on bbtime.h by Angel Diaz
  * and the modified function GetCPUTime.  Conversion into class allows for
  * multithread usage by eliminating static variables.
  * @see GetCPUTime
  */
class bbtimer 
{
public:

  /** Constructor.
    * Sets host name and initial time buffer.
    */
  bbtimer()
  {
    times(&ibuffer);
    char hostname[128];
    gethostname(hostname, 128);
    
    // truancate host name and place in string
    host = string(&hostname[0], find(&hostname[0],&hostname[127],'.'));
  } // bbtimer()

private:

  /** Output operator.
    * Time elapsed is obtained by the output operator<< to an
    * output stream such as cout.
    */
  friend ostream& operator<< (ostream&, const bbtimer&);
  
  /// Name of machine on which work is being done
  string host;

  /// Initial time buffer
  struct tms ibuffer;

}; // class bbtimer

ostream& operator<< (ostream& os, const bbtimer& B) 
{
  struct tms buffer;
  times( &buffer );
  long istime  =  buffer.tms_stime - B.ibuffer.tms_stime;
  long iutime  =  buffer.tms_utime - B.ibuffer.tms_utime;
  double stime = static_cast<double>(istime) / static_cast<double>(HZ);
  double utime = static_cast<double>(iutime) / static_cast<double>(HZ);
  long mms     = (istime + iutime) / HZ;
  long hms     = mms / 3600;
  long sms     = mms % 60;
  mms          = ( mms - hms * 3600) / 60;
  os /* << task << ": \t" */ << utime << " User, " << stime << " Sys, "
       << hms << ":" << mms << ":" << sms << " Total (on "
       << B.host << ")" << endl;
  return os;
} // operator<<

/** BlackBox timer function.
  * This is a very slight modification of Angel Diaz's 
  * GetCPUTime to fix a bug in his code, to use standard C and C++,
  * and to combine two functions into one.
  * The first call to the function sets static members record
  * host machine name, initial time called, and a record that it has
  * been called already.
  * Subsequent calls cause task name and elapsed time to be
  * printed to output stream.
  * @param  task  name of current task to be printed to output
  * @param  os    output stream to which putput is printed
  * @see bbtimer
  */
void GetCPUTime(const string& task, ostream& os = cout ) 
{
  struct tms  buffer;
  static struct tms ibuffer;
  static bool firstTime(true);
  static string host;
  
  if (firstTime)  
  {
    firstTime = false; // mark as having been here before
    times(&ibuffer);
    char hostname[128];
    gethostname( hostname, 128);
    
    // truancate host name and place in string
    host = string(&hostname[0], find(&hostname[0],&hostname[127],'.'));
  }
  else 
  { // print elapsed usage on second and subseq. calls
         // note HZ is in sys/param.h and defines the no. of ticks per clock
    times( &buffer );
    long istime =  buffer.tms_stime - ibuffer.tms_stime;
    long iutime =  buffer.tms_utime - ibuffer.tms_utime;
    double stime = static_cast<double>(istime) / static_cast<double>(HZ);
    double utime = static_cast<double>(iutime) / static_cast<double>(HZ);
    long mms    = (istime + iutime) / HZ;
    long hms    = mms / 3600;
    long sms    = mms % 60;
    mms         = ( mms - hms * 3600) / 60;
    os << task << ": \t" << utime << " User, " << stime << " Sys, "
       << hms << ":" << mms << ":" << sms << " Total (on "
       << host << ")" << endl;
  }
  return;
} // GetCPUTime()

#endif
