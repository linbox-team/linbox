/* File: src/examples/utils/fileutils.h
 * Author: William J Turner for LinBox group
 */

#ifndef _FILEUTILS_
#define _FILEUTILS_

#include <fstream>
#include <iostream>
#include <string>

/** Print error and exit.
 * Prints error message
 * "ERROR: err err2"
 * to cerr and exits the program with return code 1.
 * @param  err   string containing error message.
 * @param  err2  string containing second error message.
 */
void error(const string& err, const string& err2) 
{
  cerr << "ERROR: " << err << err2 << endl;
  exit(1);
} // error()

/** Open input stream.
 * Opens input stream for file with name filename.
 * If filename is cin, returns pointer to cin.
 * Logs messages to output stream log.
 * @return input stream pointer
 * @param  filename  name of input file to open
 * @param  log       output stream on which to log operations
 */
istream* openInput(const string& filename, ostream& log) 
{
  if (filename == string("cin")) 
  {
    log << "Using standard input cin" << endl;
    return &cin;
  } // (filename == "cin")

  log << "Opening input file " << filename << endl;
  ifstream* input = new ifstream(filename.c_str());

  if (!input->is_open()) error("cannot open input file ",filename);
    
  return input;
} // openInput()

/** Open output stream.
 * Opens output stream for file with name filename.
 * If filename is cout, returns pointer to cout.
 * Logs messages to output stream log.
 * @return output    stream pointer
 * @param  filename  name of output file to open
 * @param  log       output stream on which to log operations
 */
ostream* openOutput(const string& filename, ostream& log) 
{
  if (filename == string("cout")) 
  {
    log << "Using standard output cout" << endl;
    return &cout;
  } // (filename == "cin")

  log << "Opening output file " << filename << endl;
  ofstream* output = new ofstream(filename.c_str());

  if (!output->is_open()) error("cannot open output file ", filename);
    
  return output;
} // openOutput()

/** Open output stream for logging.
 * Opens output stream for file with name filename.
 * If filename is clog, returns pointer to clog.
 * Logs messages to clog.
 * @return output    stream pointer
 * @param  filename  name of output file to open
 * @param  log       output stream on which to log operations
 */
ostream* openLog(const string& filename) 
{
  if (filename == string("clog")) 
  {
    clog << "Using standard log clog" << endl;
    return &clog;
  } // (filename == "cin")

  clog << "Opening log file " << filename << endl;
  ofstream* log = new ofstream(filename.c_str());

  if (!log->is_open()) error("cannot open log file ", filename);
    
  return log;
} // openLog()

#endif
