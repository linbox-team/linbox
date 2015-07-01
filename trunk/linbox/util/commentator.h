/* linbox/util/commentator.h
 * Copyright (C) 1999 B. David Saunders,
 *                    Jean-Guillaume Dumas
 * Copyright (C) 2002 Bradford Hovinen
 *
 * Written by B. David Saunders <saunders@cis.udel.edu>,
 *            Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
  * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *
 * This file implements the C++ interface to commentators (for
 * providing runtime commentary to the user)
 */

/*! @file util/commentator.h
 * @ingroup util
 * @brief Reports information to the user.
 */

#ifndef __LINBOX_commentator_H
#define __LINBOX_commentator_H

#include <deque>
#include <stack>
#include <map>
#include <list>
#include <string>
#include <iostream>
#include <streambuf>
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

#ifndef DISABLE_COMMENTATOR

namespace LinBox
{
	// Forward declaration
	class MessageClass;
	class Commentator;

	// \class ActivityState commentator.h linbox/util/commentator.h
	/** @internal
	 * \brief used by commentator

	 * This stores a snapshot of the state of the commentator's activity
	 * stack, so it may be restored after an exception is thrown
	 */
	class ActivityState {
	public:

		ActivityState (void *act) :
			_act (act)
		{}

	private:

		friend class Commentator;

		void *_act;
	};

	/** @brief Give information to user during runtime.
	 * \ingroup util
	 *
	 * This object is used for reporting information about a computation to
	 * the user. Such information includes errors and warnings, descriptions
	 * of internal progress, performance measurements, and timing
	 * estimates. It also includes facilities for controlling the type and
	 * amount of information displayed.
	 *
	 * Typical usage follows the following pattern:
	 \code
	 void myFunction ()
	 {
	     commentator().start ("Doing important work", "myFunction", 100);
	     for (int i = 0; i < 100; i++) {
	         ...
	         commentator().progress ();
	     }
	     commentator().stop (MSG_DONE, "Task completed successfully");
	 }
	 \endcode
	 *
	 * In the above example, the call to commentator().start () informs the
	 * commentator that some new activity has begun. This may be invoked
	 * recursively, an the commentator keeps track of nested activities. The
	 * user may elect to disable the reporting of any activities below a
	 * certain depth of nesting. The call to commentator().stop () informs the
	 * commentator that the activity started above is finished.
	 *
	 * The call to commentator().progress () indicates that one step of the
	 * activity is complete. The commentator may then output data to the
	 * console to that effect. This allows the easy implementation of
	 * progress bars and other reporting tools.
	 *
	 * In addition, commentator().report () allows reporting general messages,
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
		/** @internal
		 * Default constructor.
		 * Constructs a commentator with default settings
		 */
		Commentator ();
		Commentator (std::ostream&);


		/** @internal
		 * Default destructor.
		*/
		virtual ~Commentator ();

		/**@internal
		 * @name Reporting facilities.
		 * @{
		*/


		/** <!--@internal-->
		 * Start an activity.
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
		void start (const char *description,
			    const char *fn = (const char *) 0,
			    unsigned long len = 0);
		void start (std::string description,
			    const char *fn = (const char *) 0,
			    unsigned long len = 0);


		/** @internal
		 * Start a new iteration.
		 * This is a convenience function to indicate that an iteration
		 * has started
		 * @param iter Number of the iteration
		 * @param len Number of items involved in this activity. Used
		 *            for progress reporting; use 0 if this is not
		 *            applicable to the situation in question.
		 *            (default 0)
		 */
		void startIteration (unsigned int iter, unsigned long len = 0);

		/** <!--@internal-->
		 * Stop an activity.
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
		void stop (const char *msg = MSG_DONE,
			   const char *long_msg = (const char *) 0,
			   const char *fn = (const char *) 0);

		/** <!--@internal-->
		 * Report progress in the current activity.
		 * Inform the commentator that k steps of the current activity
		 * have been completed.
		 * @param k Number of steps completed; use -1 to increment the
		 *          current counter
		 * @param len Total length of the operation; use -1 to use the
		 *            existing length. This allows updating of inexact
		 *            estimates of the number of steps.
		 */
		void progress (long k = -1, long len = -1);

		/** @internal
		 * Message level.
		 * Some default settings to use for the message level
		 */
		enum MessageLevel {
			LEVEL_ALWAYS       =  0,
			LEVEL_IMPORTANT    =  1,
			LEVEL_NORMAL       =  2,
			LEVEL_UNIMPORTANT  =  3
		};

		/** <!--@internal-->
		 * Basic reporting.
		 * Send some report string to the commentator
		 * @param level Level of detail of the message
		 * @param msg_class Type of message
		 * @return A reference to the stream to which to output data
		 */
		std::ostream &report (long level = LEVEL_IMPORTANT, const char *msg_class = INTERNAL_DESCRIPTION);

		/** @internal
		 * Indent to the correct column on the given string.
		*/
		void indent (std::ostream &stream) const;

		//@} Reporting facilities


		/** @internal
		 * @name Activity stack restoration.
		 *
		 * The methods below facilitate restoring the activity stack
		 * after an exception has been thrown. If user code wishes to
		 * catch an exception, it should invoke the method
		 * saveActivityState before the try block, storing the result in
		 * an ActivityState object. Then it may invoke
		 * restoreActivityState, passing the ActivityState object, in
		 * each of the catch blocks.
		 */

		//@{

		/** @internal
		 * Save activity state.
		 *
		 * Saves a copy of the activity state and returns it to the
		 * caller. The caller need only pass this object back to
		 * \ref restoreActivityState to return the commentator's
		 * activity stack to the state it was when saveActivityState was
		 * called.
		 *
		 * @return ActivityState object
		 */

		ActivityState saveActivityState () const
		{
		       	return ActivityState (_activities.top ());
		}

		/** @internal
		 * Restore activity state.
		 *
		 * Restores the activity state to the point it was when the
		 * given ActivityState object was passed. Note that this
		 * function assumes that the commentator is currently in a more
		 * deeply-nested configuration than when the object was created;
		 * if it is not, the method will give up.
		 */

		void restoreActivityState (ActivityState state);

		/** @internal
		 * @name Configuration
		*/

		//@{

		/** @internal
		 * Output format.
		 * OUTPUT_CONSOLE - output human-readable and tailor-made for the console
		 * OUTPUT_PIPE - output made for piping into other processes
		 *
		 * This affects only the brief report
		 */
		enum OutputFormat {
		       	OUTPUT_CONSOLE,
		       	OUTPUT_PIPE
		};

		enum EstimationMethod {
			BEST_ESTIMATE,
			POLY_ESTIMATE,
			EXPO_ESTIMATE
		};

		enum TimingDelay {
			SHORT_TIMING,
			LONG_TIMING
		};

		/** @internal
		 * Set maximum message depth.
		 * Sets the maximum activity depth, as defined by
		 * \c Commentator::start and \c Commentator::stop, at which messages
		 * of all classes will be printed.
		 * @param depth Maximum depth at which to print messages; set to
		 *              -1 to print messages at any depth and 0 to
		 *              disable entirely
		 */
		void setMaxDepth (long depth);

		/** @internal
		 * Set maximum detail level.
		 * Sets the maximum detail level at which to print messages
		 * @param level Maximum detail level at which to print messages;
		 *              set to -1 to print messages at all levels and 0
		 *              to disable entirely
		 */
		void setMaxDetailLevel (long level);

		/** @internal
		 * Register a new message class.
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
		MessageClass &registerMessageClass (const char *msg_class,
						    std::ostream &stream,
						    unsigned long max_depth = 1,
						    unsigned long max_level = 2);

		/** @internal
		 * Clone an existing message class.
		 * Clone the message class to construct a new message class with
		 * the same parameters.
		 * @param new_msg_class Name of the new class
		 * @param msg_class Name of class to clone
		 * @return A reference to the new message class object
		 */
		MessageClass &cloneMessageClass (const char *new_msg_class, const char *msg_class);

		/** @internal
		 * Clone an existing message class, specifying a new output stream.
		 * Clone the message class to construct a new message class with
		 * the same parameters.
		 * @param new_msg_class Name of the new class
		 * @param msg_class Name of class to clone
		 * @param stream New stream to which to direct output
		 * @return A reference to the new message class object
		 */
		MessageClass &cloneMessageClass (const char *new_msg_class,
						 const char *msg_class,
						 std::ostream &stream);

		/** @internal
		 * Retrieve a message class by name.
		 * @param msg_class Name of message class
		 * @return Reference to message class object
		 */
		MessageClass &getMessageClass (const char *msg_class);

		/** @internal
		 * Precise control over printing.
		 * Specifies that all messages up to the given depth and the
		 * given detail level should be printed
		 *
		 * This is similar to MessageClass::setPrintParameters but
		 * affects all message classes simultaneously.
		 * @param depth Depth up to which directive is valid
		 * @param level Detail level up to which directive is valid
		 * @param fn Fully qualified name of the function for which this
		 *           is valid; may be 0, in which case this is valid for
		 *           everything
		 */
		void setPrintParameters (unsigned long depth,
					 unsigned long level,
					 const char *fn = (const char *) 0);

		/** @internal
		 * Set parameters for the brief report.
		 * @param format Output format
		 * @param show_timing Show the CPU time for each toplevel activity
		 * @param show_progress Show a counter of the progress as each
		 *                      activity progresses
		 * @param show_est_time Show estimated time remaining
		 */
		void setBriefReportParameters (OutputFormat format,
					       bool show_timing,
					       bool show_progress,
					       bool show_est_time);

		/** @internal
		 * Determine whether a message will be printed.
		 * @param depth Activity depth
		 * @param level Message level
		 * @param msg_class Type of message
		 * @param fn Fully qualified function name (0 if not applicable)
		 * @return true if the message with the given characteristics
		 *         will be printed as things are currently configured
		 */
		bool isPrinted (unsigned long depth,
				unsigned long level,
				const char *msg_class,
				const char *fn = (const char *) 0);

		/** @internal
		 * Determine whether a message will be printed.
		 * This variant uses the current activity depth rather than
		 * specifying it explicitly.
		 * @param level Message level
		 * @param msg_class Type of message
		 * @param fn Fully qualified function name (0 if not applicable)
		 * @return true if the message with the given characteristics
		 *         will be printed as things are currently configured
		 */
		bool isPrinted (unsigned long level,
				const char *msg_class,
				const char *fn = (const char *) 0)
		{
		       	return isPrinted (_activities.size (), level, msg_class, fn);
		}

		/** @internal
		 * Determine whether the stream given is the null stream.
		 * @param str Stream to check
		 * @return true if stream is the null stream; false otherwise
		 */
		bool isNullStream (const std::ostream &str)
		{
			return &str == &cnull;
		}

		/** @internal
		 * Set output stream for brief report.
		 * @param stream Output stream
		 */
		void setBriefReportStream (std::ostream &stream);

		/** @internal
		 * Set output stream for all reports other than the brief.
		 * report
		 * @param stream Output stream
		 */
		void setReportStream (std::ostream &stream);

		/** @internal
		 * Set the output stream for a given message class.
		 * @param msg_class Message class
		 * @param stream Output stream
		 */
		void setMessageClassStream (const char *msg_class, std::ostream &stream);

		/** @internal
		 * Set default report file.
		 * @param filename of file
		 */
		void setDefaultReportFile (const char *filename);

		//@} Configuration

		/** @internal
		 * @name Legacy commentator interface
		 * These routines provide compatibility with the old commentator
		 * interface. They are deprecated.
		 * @deprecated
		 */

		//@{

		/** @internal
		 * Start an activity.
		 * @param id String identifier of activity
		 * @param msg Message to print
		 * @param msglevel Level of message
		 * @param msgclass Class of message
		 */
		void start (const char *id, const char *msg, long msglevel, const char *msgclass)
		{
			start (id);
			report ((MessageLevel) msglevel, msgclass) << msg << std::endl;
		}

		void start (std::string id, const char *msg, long msglevel, const char *msgclass)
		{
			start(id.c_str(),msg,msglevel,msgclass);
		}

		/** @internal
		 * Stop an activity.
		 * @param msg Message to print
		 */
		 /* @param msglevel Level of message
		 * @param msgclass Class of message
		 * @param time_type Type of timing to use
		 */
		void stop (const char *msg
			   , long //msglevel
			   , const char * //msgclass
			   , long //time_type
			  )
		{
			stop (msg);
		}

		/** @internal
		 * Report progress.
		 * @param msg Message to print
		 * @param msglevel Level of message
		 * @param k Number of steps completed
		 * @param n Total number of steps in operation
		 */
		void progress (const char *msg, long msglevel, long k, long n)
		{
			progress (k, n);
			report ((MessageLevel) msglevel, INTERNAL_DESCRIPTION) << msg << std::endl;
		}

		/** @internal
		 * General reporting.
		 * @param msg Message to print
		 * @param msglevel Level of message
		 * @param msgclass Class of message
		 */
		void report (const char *msg, long msglevel, const char *msgclass)
		{
			report ((MessageLevel) msglevel, msgclass) << msg << std::endl;
		}

		/** @internal
		 * Test whether message is printed.
		 * @param msglevel Level of message
		 * @param msgclass Class of message
		 */
		bool printed (long msglevel, const char *msgclass)
		{
			return isPrinted (_activities.size (), (MessageLevel) msglevel, msgclass);
		}

		//@} Legacy commentator interface

		/** @internal
		 * Use this stream to disable a message class entirely.
		*/
		std::ofstream cnull;

	private:
#if 0
		// Null std::ostream prints nothing
		struct nullstreambuf : public std::streambuf {
			nullstreambuf() {};
			// GV modidied seek_dir twice
			std::streampos seekoff(std::streambuf::off_type, std::ios::seekdir, std::ios::openmode)
			{
				return 0;
			}
			std::streampos seekpos(std::streambuf::pos_type, std::ios::openmode)
			{
				return 0;
			}
			std::streampos sys_seek(std::streambuf::off_type, std::ios::seekdir)
			{
				return 0;
			}
			std::streamsize showmanyc ()
			{
				return 0;
			}
			void imbue(const std::locale &loc)
			{}
		};
#endif

	protected:
		struct StepsAndTime {
			StepsAndTime (long k = 0, double t = 0.0) :
				_time(t), _steps(k)
			{}

			double                   _time;
			long                     _steps;
		};

		typedef std::deque<StepsAndTime> Estimator;

		struct Activity {
			Activity (const char *desc, const char *fn, unsigned long len) :
				_desc (desc), _fn (fn), _len (len), _progress (0)
			{}

			const char              *_desc;
			const char              *_fn;
			unsigned long            _len;
			unsigned long            _progress;
		 Givaro::Timer                    _timer;
			Estimator                _estimate;
		};

		std::stack<Activity *>           _activities;      // Stack of activity structures

		struct C_str_Less {
			bool operator() (const char* x, const char * y) const {
				return strcmp(x,y)<0;
			}
		};

		std::map<const char *, MessageClass *, C_str_Less>
		_messageClasses;
		EstimationMethod                 _estimationMethod;     // Activity estimator

		// Brief report parameters
		OutputFormat                     _format;
		bool                             _show_timing;
		bool                             _show_progress;
		bool                             _show_est_time;

		unsigned int                     _last_line_len;      // Length of last printed line in the brief report

		std::ofstream                    _report;

		std::string                      _iteration_str;     // String referring to current iteration -- HACK

		// Functions for the brief report
		virtual void printActivityReport  (Activity &activity);
		virtual void updateActivityReport (Activity &activity);
		virtual void finishActivityReport (Activity &activity, const char *msg);
	};

	/** @internal
	 * Message class object.
	 * This object encapsulates the configuration of a given message class
	 */
	class MessageClass {
	public:
		friend class Commentator;

		/**  @internal
		 * Constructor.
		 *
		 * Constructs a new MessageClass object with the given type,
		 * outputing to the given stream. All other parameters are set
		 * to factory defaults.
		 * @param comm Commentator for this message class
		 * @param msg_class Name of message class
		 * @param stream Output stream to which to send messages of this
		 *               class
		 * @param max_depth Default maximal recursion depth at which to
		 *                  print messages of this class (default 2)
		 * @param max_level Default maximal detail level at which to
		 *                  print messages of this class (default 1)
		 */
		MessageClass (const Commentator &comm,
			      const char *msg_class,
			      std::ostream &stream,
			      unsigned long max_depth = 1,
			      unsigned long max_level = 2);

		/** @internal
		 * Copy constructor.
		*/
		MessageClass (const MessageClass &message_class) :
			_msg_class (message_class._msg_class),
			_smart_streambuf (message_class._smart_streambuf.comm (),
					  message_class._smart_streambuf.stream ()),
			_stream (&_smart_streambuf),
			_configuration (message_class._configuration),
			_max_level (message_class._max_level),
			_max_depth (message_class._max_depth)
		{}

		/** @internal
		 * Set maximum message depth.
		 * Sets the maximum activity depth, as defined by
		 * Commentator::start and Commentator::stop, at which messages
		 * of this class will be printed.
		 * @param depth Maximum depth at which to print messages; set to
		 *              -1 to print messages at any depth and 0 to
		 *              disable entirely
		 */
		void setMaxDepth (long depth);

		/** @internal
		 * Set maximum detail level.
		 * Sets the maximum detail level at which to print messages
		 * @param level Maximum detail level at which to print messages;
		 *              set to -1 to print messages at all levels and 0
		 *              to disable entirely
		 */
		void setMaxDetailLevel (long level);

		/** @internal
		 * Precise control over printing.
		 * Specifies that all messages up to the given depth and the
		 * given activity level should be printed
		 * @param depth Maximum depth at which to print
		 * @param level Maximum detail level at which to print
		 * @param fn Fully qualified name of the function for which this
		 *           is valid; may be 0, in which case this is valid for
		 *           everything
		 */
		void setPrintParameters (unsigned long depth, unsigned long level, const char *fn);

		/** @internal
		 * Determine whether a given message will be printed when of
		 * this message class.
		 * @param depth Activity depty of the message
		 * @param level Detail level of the message
		 * @param fn Fully qualified function name (0 if not applicable;
		 *           default 0)
		 * @return true if a message of the given characteristics is
		 *         printed
		 */
		bool isPrinted (unsigned long depth, unsigned long level, const char *fn = (const char *) 0);

	private:
		typedef std::map <const char *, std::list<std::pair <unsigned long, unsigned long> > > Configuration;

		class smartStreambuf : public std::streambuf {
			const Commentator &_comm;
			std::ostream &_stream;
			bool _indent_next;

		public:
			smartStreambuf (const Commentator &Comm, std::ostream & Stream) :
				_comm (Comm), _stream (Stream), _indent_next (true)
			{}

			int sync ();
			int overflow (int ch);
			std::streamsize xsputn (const char *text, std::streamsize n);
			int writeData (const char *text, std::streamsize n);

			const Commentator &comm () const
			{
			       	return _comm;
			}
			std::ostream &stream () const
			{
			       	return _stream;
			}
		};

		const char     *_msg_class; // Name of this message class
		smartStreambuf  _smart_streambuf;
		std::ostream    _stream;    // Stream to which to output data

		Configuration   _configuration;

		unsigned long   _max_level;
		unsigned long   _max_depth;

		// Constructor that gives existing configuration to copy
		MessageClass (const Commentator &comm, const char *msg_class, std::ostream &stream, Configuration configuration);

		void fixDefaultConfig ();
		bool checkConfig (std::list <std::pair <unsigned long, unsigned long> > &config, unsigned long depth, unsigned long level);
		void dumpConfig () const;   // Dump the contents of configuration to stderr
	};

	// Default static commentator is now common to enabled or disabled
// 	extern Commentator commentator;
}

// #ifdef LinBoxSrcOnly
// #include "linbox/util/commentator.C"
// #endif
#include "linbox/util/commentator.inl"

#define aside commentator().report(Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTON)
// Usage: "aside << stuff" or "ostream& report = aside; report << stuff"

#else //DISABLE_COMMENTATOR

#define aside NoStream()
struct NoStream {};
template<typename T>
NoStream & operator<< (NoStream& o, const T & x) { return o; }
#if 0
#  define Commentator CommentatorDisabled
#  define MessageClass MessageClassDisabled
#  define commentator commentatorDisabled
#endif

// This provides a "null" commentator that should compile totally out of the
// program when DISABLE_COMMENTATOR is defined. All code making use of the
// commentator should disappear.


namespace LinBox
{
	// Forward declaration
	class Commentator;

	class MessageClass {
	public:
		friend class Commentator;

		inline MessageClass (const char *, std::ostream &, unsigned long = 1, unsigned long = 2)
	       	{}
		inline MessageClass ()
	       	{}
		inline void setMaxDepth (long )
	       	{}
		inline void setMaxDetailLevel (long )
	       	{}
		inline void setPrintParameters (unsigned long, unsigned long , const char *)
	       	{}
		inline bool isPrinted (unsigned long , unsigned long , const char * = (const char *) 0)
		{ return false; }
	};

	class ActivityState {};

	class Commentator {
	public:
#if 0
		inline Commentator () :
			cnull (new nullstreambuf)
		{}
		inline Commentator () :
			cnull ("/dev/null")
		{}
#endif
		inline Commentator (std::ostream& out = std::cerr) :
			cnull (0)
		{}
		inline  ~Commentator ()
	       	{}
		inline void start (const char *, const char * = (const char *) 0, unsigned long = 0)
	       	{}
		inline void startIteration (unsigned int , unsigned long = 0)
		{}
		inline void stop (const char *, const char * = (const char *) 0, const char * = (const char *) 0)
		{}
		inline void progress (long = -1, long = -1)
	       	{}

		enum MessageLevel {
			LEVEL_ALWAYS       =  0,
			LEVEL_IMPORTANT    =  1,
			LEVEL_NORMAL       =  2,
			LEVEL_UNIMPORTANT  =  3
		};

		inline std::ostream &report (long level = LEVEL_IMPORTANT, const char *msg_class = INTERNAL_DESCRIPTION)
		{
		       	return cnull;
		}

		inline void indent (std::ostream &)
		{}

		enum OutputFormat {
		       	OUTPUT_CONSOLE,
			OUTPUT_PIPE
		};

		enum EstimationMethod {
			BEST_ESTIMATE,
			POLY_ESTIMATE,
			EXPO_ESTIMATE
		};

		enum TimingDelay {
			SHORT_TIMING,
			LONG_TIMING
		};

		inline void setMaxDepth (long )
	       	{}
		inline void setMaxDetailLevel (long )
		{}
		inline MessageClass &registerMessageClass (const char *, std::ostream &, unsigned long = 1, unsigned long = 2)
		{ return _msgcls; }
		inline MessageClass &cloneMessageClass (const char *, const char *)
		{ return _msgcls; }
		inline MessageClass &cloneMessageClass (const char *, const char *, std::ostream &)
		{ return _msgcls; }
		inline MessageClass &getMessageClass (const char *)
		{ return _msgcls; }
		inline void setPrintParameters (unsigned long , unsigned long , const char * = (const char *) 0)
	       	{}
		inline void setBriefReportParameters (OutputFormat , bool , bool , bool )
	       	{}
		inline bool isPrinted (unsigned long , unsigned long , const char *, const char * = (const char *) 0)
		{ return false; }
		inline bool isPrinted (unsigned long , const char *, const char * = (const char *) 0)
		{ return false; }
		inline bool isNullStream (const std::ostream &)
		{ return true; }
		inline void setBriefReportStream (std::ostream &)
		{}
		inline void setReportStream (std::ostream &)
	       	{}
		inline void setMessageClassStream (const char *, std::ostream &)
	       	{}
		inline void setDefaultReportFile (const char *)
	       	{}
		inline void start (const char *, const char *, long , const char *)
	       	{}
		inline void stop (const char *, long , const char *, long)
	       	{}
		inline void progress (const char *, long , long , long )
	       	{}
		inline void report (const char *, long , const char *)
	       	{}
		inline bool printed (long , const char *)
		{ return false; }

		ActivityState saveActivityState () const
		{ return ActivityState(); }
		void restoreActivityState (ActivityState state)
		{}

		std::ofstream cnull;

	private:
#if 0
		// Null std::ostream prints nothing
		struct nullstreambuf : public std::streambuf {
			nullstreambuf () {};
			~nullstreambuf () {};
			inline std::streampos seekoff (std::streambuf::off_type, std::ios::seekdir, std::ios::openmode)
			{
			       	return 0;
			}
			inline std::streampos seekpos (std::streambuf::pos_type, std::ios::openmode)
			{
			       	return 0;
			}
			inline std::streampos sys_seek (std::streambuf::off_type, std::ios::seekdir)
			{
			       	return 0;
			}
			inline std::streamsize showmanyc ()
			{
			       	return 0;
			}
			inline void imbue (const std::locale &loc) {}
		};
#endif

		MessageClass _msgcls;
	}; // class (disabled) Commentator

// 	// Default global commentator
// 	extern Commentator commentator;
// 	//static Commentator commentator;
}

#endif // DISABLE_COMMENTATOR

namespace LinBox
{
	// Default static commentator
    Commentator& commentator() {
        static Commentator internal_static_commentator;
        return internal_static_commentator;
    }
    Commentator& commentator(std::ostream& stream) {
        static Commentator internal_static_commentator(stream);
        return internal_static_commentator;
    }
}




#include "linbox/util/args-parser.h"
namespace LinBox
{
    void parseArguments (int argc, char **argv, Argument *args, bool printDefaults = true) {
        for (int i = 1; i < argc; ++i) {
            if (argv[i][0] == '-') {
                if (argv[i][1] == 0) {
                    LinBox::commentator().setReportStream (std::cout);
                    LinBox::commentator().setBriefReportStream (std::cout);
                } else {
                        // Skip the argument next to "-xxx"
                        // except if next argument is a switch
                    if ( ((i+1) < argc) &&
                         (argv[i+1][0] != '-') ) {
                        ++i;
                    }
                }
            } else {
                LinBox::commentator().setDefaultReportFile (argv[i]);
                LinBox::commentator().setBriefReportStream(std::cout);
            }
        }
        FFLAS::parseArguments(argc,argv,args,printDefaults);
    }
}

#endif // __LINBOX_commentator_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
