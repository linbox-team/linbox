/* -*- mode: c; style: linux -*- */

/* linbox/util/commentator.h
 * Copyright (C) 1999 B. David Saunders,
 *                    Jean-Guillaume Dumas
 * Copyright (C) 2002 Bradford Hovinen
 *
 * Written by B. David Saunders <saunders@cis.udel.edu>,
 *            Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 *
 * This file implements the C++ interface to commentators (for 
 * providing runtime commentary to the user)
 */

#ifndef __COMMENTATOR_H
#define __COMMENTATOR_H

#include <deque>
#include <stack>
#include <map>
#include <list>
#include <string>
#include <iostream>
#include <fstream>
#include <cstring>

#include "linbox/util/timer.h"

#ifndef MAX
#  define MAX(a,b) (((a) > (b)) ? (a) : (b))
#endif

// Definitions for various common message classes
#define BRIEF_REPORT         "Brief report"
#define PROGRESS_REPORT      "Progress report"
#define TIMING_MEASURE       "Timing measure"
#define TIMING_ESTIMATE      "Timing estimate"
#define PARTIAL_RESULT       "Partial result"
#define INTERNAL_WARNING     "Internal warning"
#define INTERNAL_ERROR       "Internal error"
#define INTERNAL_DESCRIPTION "Internal description"

// Definitions for some common termination messages
#define MSG_OK               "ok"
#define MSG_DONE             "done"
#define MSG_PASSED           "passed"
#define MSG_FAILED           "FAILED"

#define MSG_STATUS(ret) (ret ? MSG_PASSED : MSG_FAILED)

namespace LinBox 
{
	// Forward declaration
	class MessageClass;

	// Utility object needed for associative containers
	struct LessThanString
	{
		bool operator () (const char *str1, const char *str2) const 
			{ return strcmp (str1, str2) < 0; }
	};

	/** Commentator object
	 * This object is used for reporting information about a computation to
	 * the user. Such information includes errors and warnings, descriptions
	 * of internal progress, performance measurements, and timing
	 * estimates. It also includes facilities for controlling the type and
	 * amount of information displayed.
	 *
	 * Typical usage follows the following pattern:
	 * @code{
	 *   void myFunction () {
	 *       commentator.start ("Doing important work", "myFunction", 100);
	 *       for (int i = 0; i < 100; i++) {
	 *           ...
	 *           commentator.progress ();
	 *       }
	 *       commentator.stop (MSG_DONE, "Task completed successfully");
	 *   }
	 * @}
	 *
	 * In the above example, the call to commentator.start () informs the
	 * commentator that some new activity has begun. This may be invoked
	 * recursively, an the commentator keeps track of nested activities. The
	 * user may elect to disable the reporting of any activities below a
	 * certain depth of nesting. The call to commentator.stop () informs the
	 * commentator that the activity started above is finished.
	 *
	 * The call to commentator.progress () indicates that one step of the
	 * activity is complete. The commentator may then output data to the
	 * console to that effect. This allows the easy implementation of
	 * progress bars and other reporting tools.
	 *
	 * In addition, commentator.report () allows reporting general messages,
	 * such as warnings, errors, and descriptions of internal progress.
	 *
	 * By default, there are two reports: a brief report that outputs to
	 * cout, and a detailed report that outputs to a file that is
	 * specified. If no file is specified, the detailed report is thrown
	 * out. The brief report is intended either for human consumption or to
	 * be piped to some other process. Therefore, there are two output
	 * formats for that report. One can further customize the report style
	 * by inheriting the commentator object.
	 *
	 * The commentator allows very precise control over what gets
	 * printed. See the Configuration section below.
	 */

	class Commentator {
	    public: 
		/** Default constructor
		 * Constructs a commentator with default settings
		 */
		Commentator ();

		/** Default destructor
		 */
		virtual ~Commentator () {}

		/** @name Reporting facilities
		 */

		//@{

		/** Start an activity
		 * Inform the commentator that some new activity has begun. This
		 * is typically called at the beginning of a library function.
		 * @param description A human-readable text description of the
		 *                    activity
		 * @param fn The fully-qualified name of the function. Use 0 to
		 *           use the same function as the surrounding activity
		 *           (default 0)
		 * @param len Number of items involved in this activity. Used
		 *            for progress reporting; use 0 if this is not
		 *            applicable to the situation in question.
		 *            (default 0)
		 */
		void start (const char *description, const char *fn = (const char *) 0, unsigned long len = 0);

		/** Stop an activity
		 * Inform the commentator that the current activity has
		 * finished.
		 * @param msg A short (one word) message describing the
		 *            termination, e.g. "passed", "FAILED", "ok", etc.
		 * @param long_msg A longer message describing the nature of
		 *                 the termination. May be 0, in which case msg
		 *                 is used.
		 * @param fn Name of the function whose activity is
		 *           complete. This is intended for checking to make
		 *           sure that each call to start () is matched with a
		 *           call to stop (). It is optional, and has no other
		 *           effect.
		 */
		void stop (const char *msg, const char *long_msg = (const char *) 0, const char *fn = (const char *) 0);

		/** Report progress in the current activity
		 * Inform the commentator that k steps of the current activity
		 * have been completed.
		 * @param k Number of steps completed; use -1 to increment the
		 *          current counter
		 * @param len Total length of the operation; use -1 to use the
		 *            existing length. This allows updating of inexact
		 *            estimates of the number of steps.
		 */
		void progress (long k = -1, long len = -1);

		/** Message level
		 * Some default settings to use for the message level
		 */
		enum MessageLevel {
			LEVEL_ALWAYS       =  1,
			LEVEL_IMPORTANT    =  2,
			LEVEL_NORMAL       =  3,
			LEVEL_UNIMPORTANT  =  4,
		};

		/** Basic reporting
		 * Send some report string to the commentator
		 * @param level Level of detail of the message
		 * @param msg_class Type of message
		 * @return A reference to the stream to which to output data
		 */
		ostream &report (long level, const char *msg_class);

		/** Indent to the correct column on the given string
		 */
		void indent (ostream &stream);

		//@} Reporting facilities

		/** @name Configuration
		 */

		//@{

		/** Output format
		 * OUTPUT_CONSOLE - output human-readable and tailor-made for the console
		 * OUTPUT_PIPE - output made for piping into other processes
		 *
		 * This affects only the brief report
		 */
		enum OutputFormat
			{ OUTPUT_CONSOLE, OUTPUT_PIPE };

		enum EstimationMethod {
			BEST_ESTIMATE,
			POLY_ESTIMATE,
			EXPO_ESTIMATE
		};

		enum TimingDelay {
			SHORT_TIMING,
			LONG_TIMING
		};

		/** Register a new message class
		 * Register some new message class with the commentator.
		 * @param msg_class Name of message class
		 * @param stream Stream to which to send data of this type
		 * @param max_depth Default maximum recursion depth at which to
		 *                  print messages of this class (default 1)
		 * @param max_level Default maximum detail level at which to
		 *                  print messages of this class (default 2, or
		 *                  LEVEL_IMPORTANT)
		 * @return A reference to the new message class object
		 */
		MessageClass &registerMessageClass (const char *msg_class, ostream &stream, unsigned long max_depth = 1, unsigned long max_level = 2);

		/** Clone an existing message class
		 * Clone the message class to construct a new message class with
		 * the same parameters.
		 * @param new_msg_class Name of the new class
		 * @param msg_class Name of class to clone
		 * @return A reference to the new message class object
		 */
		MessageClass &cloneMessageClass (const char *new_msg_class, const char *msg_class);

		/** Clone an existing message class, specifying a new output stream
		 * Clone the message class to construct a new message class with
		 * the same parameters.
		 * @param new_msg_class Name of the new class
		 * @param msg_class Name of class to clone
		 * @param stream New stream to which to direct output
		 * @return A reference to the new message class object
		 */
		MessageClass &cloneMessageClass (const char *new_msg_class, const char *msg_class, ostream &stream);

		/** Retrieve a message class by name
		 * @param msg_class Name of message class
		 * @return Reference to message class object
		 */
		MessageClass &getMessageClass (const char *msg_class);

		/** Precise control over printing
		 * Sets whether not messages between two given activity depths
		 * and up to a given detail level are printed, optionally
		 * specifying a specific function. This gives very precise
		 * control over what gets printed.
		 *
		 * This is similar to MessageClass::setPrintParameters but
		 * affects all message classes simultaneously.
		 * @param low_depth First activity depth in the range
		 * @param high_depth Last activity depth in the range; use -1
		 *                   for all depths below low_depth
		 * @param max_level Maximum detail level at which to print
		 * @param fn Fully qualified name of the function for which this
		 *           is valid; may be 0, in which case this is valid for
		 *           everything
		 */
		void setPrintParameters (long low_depth, long high_depth, long max_level, const char *fn = (const char *) 0);

		/** Set parameters for the brief report
		 * @param format Output format
		 * @param show_timing Show the CPU time for each toplevel activity
		 * @param show_progress Show a counter of the progress as each
		 *                      activity progresses
		 * @param show_est_time Show estimated time remaining
		 */
		void setBriefReportParameters (OutputFormat format, bool show_timing, bool show_progress, bool show_est_time);

		/** Determine whether a message will be printed
		 * @param depth Activity depth
		 * @param level Message level
		 * @param msg_class Type of message
		 * @param fn Fully qualified function name (0 if not applicable)
		 * @return true if the message with the given characteristics
		 *         will be printed as things are currently configured
		 */
		bool isPrinted (long depth, long level, const char *msg_class, const char *fn = (const char *) 0);

		/** Set output stream for brief report
		 * @param stream Output stream
		 */
		void setBriefReportStream (ostream &stream);

		/** Set output stream for all reports other than the brief
		 * report
		 * @param stream Output stream
		 */
		void setReportStream (ostream &stream);

		/** Set the output stream for a given message class
		 * @param msg_class Message class
		 * @param stream Output stream
		 */
		void setMessageClassStream (const char *msg_class, ostream &stream);

		/** Set default report file
		 * @param name of file
		 */
		void setDefaultReportFile (const char *filename);

		//@} Configuration

	    protected:
		struct StepsAndTime {
			StepsAndTime (long k, double t)
				: _time(t), _steps(k) {}

			double                   _time;
			long                     _steps;
		};
        
		typedef std::deque<StepsAndTime> Estimator;

		struct Activity {
			Activity (const char *desc, const char *fn, unsigned long len) 
				: _desc (desc), _fn (fn), _len (len), _progress (0) {}

			const char              *_desc;
			const char              *_fn;
			unsigned long            _len;
			unsigned long            _progress;
			Timer                    _timer;
			Estimator                _estimate;
		};

		std::stack<Activity *>           _activities;      // Stack of activity structures

		std::map<const char *, MessageClass *, LessThanString>
			                         _messageClasses;
		EstimationMethod                 _estimationMethod;     // Activity estimator

		// Brief report parameters
		OutputFormat                     _format;
		bool                             _show_timing;
		bool                             _show_progress;
		bool                             _show_est_time;

		ofstream                         _report;

		// Functions for the brief report
		virtual void printActivityReport  (Activity &activity);
		virtual void updateActivityReport (Activity &activity);
		virtual void finishActivityReport (Activity &activity, const char *msg);

	    private:
		// Null ostream prints nothing
		struct nullstreambuf : public streambuf {
			nullstreambuf() {};
			streampos seekoff(long long, ios::seek_dir, int) {return 0;}
			streampos seekpos(long long, int) {return 0;}
			streampos sys_seek(long long, ios::seek_dir) {return 0;}
			int showmanyc(void) {return 0;}
			void imbue(void *) {}
		};

		ostream cnull;
	};

	/** Message class object
	 * This object encapsulates the configuration of a given message class
	 */
	class MessageClass {
	    public:
		friend class Commentator;

		/** Constructor
		 * Constructs a new MessageClass object with the given type,
		 * outputing to the given stream. All other parameters are set
		 * to factory defaults.
		 * @param msg_class Name of message class
		 * @param stream Output stream to which to send messages of this
		 *               class
		 * @param max_depth Default maximal recursion depth at which to
		 *                  print messages of this class (default 2)
		 * @param max_level Default maximal detail level at which to
		 *                  print messages of this class (default 1)
		 */
		MessageClass (const char *msg_class, ostream &stream, unsigned long max_depth = 1, unsigned long max_level = 2);

		/** Set maximum message depth
		 * Sets the maximum activity depth, as defined by
		 * Commentator::start and Commentator::stop, at which messages
		 * of this class will be printed.
		 * @param depth Maximum depth at which to print messages; set to
		 *              -1 to print messages at any depth and 0 to
		 *              disable entirely
		 */
		void setMaxDepth (long depth);

		/** Set maximum detail level
		 * Sets the maximum detail level at which to print messages
		 * @param level Maximum detail level at which to print messages;
		 *              set to -1 to print messages at all levels and 0
		 *              to disable entirely
		 */
		void setMaxDetailLevel (long level);

		/** Precise control over printing
		 * Sets whether not messages between two given activity depths
		 * and up to a given detail level are printed, optionally
		 * specifying a specific function. This gives very precise
		 * control over what gets printed.
		 * @param low_depth First activity depth in the range
		 * @param high_depth Last activity depth in the range; use -1
		 *                   for all depths below low_depth
		 * @param max_level Maximum detail level at which to print
		 * @param fn Fully qualified name of the function for which this
		 *           is valid; may be 0, in which case this is valid for
		 *           everything
		 */
		void setPrintParameters (long low_depth, long high_depth, long max_level, const char *fn);

		/** Determine whether a given message will be printed when of
		 * this message class
		 * @param depth Activity depty of the message
		 * @param level Detail level of the message
		 * @param fn Fully qualified function name (0 if not applicable;
		 *           default 0)
		 * @return true if a message of the given characteristics is
		 *         printed
		 */
		bool isPrinted (long depth, long level, const char *fn = (const char *) 0);

	    private:
		typedef std::map <const char *, std::list<std::pair <unsigned long, unsigned long> >, LessThanString> Configuration;

		const char    *_msg_class; // Name of this message class
		ostream       &_stream;    // Stream to which to output data

		Configuration  _configuration;

		unsigned long  _max_level;
		unsigned long  _max_depth;

		// Constructor that gives existing configuration to copy
		MessageClass (const char *msg_class, ostream &stream, Configuration configuration);

		void fixDefaultConfig ();
		bool checkConfig (list <pair <unsigned long, unsigned long> > &config, long depth, long level);
		void dumpConfig () const;   // Dump the contents of configuration to stderr
	};

	// Default global commentator
	extern Commentator commentator;
}

#endif // __COMMENTATOR_H
