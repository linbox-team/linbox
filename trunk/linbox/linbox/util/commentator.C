/* -*- mode: C++; style: linux -*- */

/* linbox/util/commentator.C
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

#include "linbox-config.h"

#include <string>
#include <strstream>
#include <cmath>
#include <cstring>
#include <cstdlib>

#include "linbox/util/commentator.h"
#include "linbox/util/debug.h"

namespace LinBox 
{
	using namespace std;

	// -----------------------------------------------------
	// Mathematical routines
	// -----------------------------------------------------
	double nroot (double a, long r, double precision)
	{
		long rm = r - 1 ;
		double c1 = double (rm) / double (r), c2 = a / double (r);
		double g = 1, pgr = 1, err = a - 1;
    
		while (err > precision) {
			g = g * c1 + c2 / pgr;
			pgr = pow (g, (double) rm);
			err = a - pgr * g;
			if (err < 0)
				err = -err;
		}

		return g;
	}

	long isnpower (long& l, long a)
	{
		long r = 2;
		double g;

		while ((g = nroot (a, r, 0.1)) >= 2) {
			l = (long) floor (g);
			if (g-double (l) > 0.1)
				++l;
			if (pow ((double) l, (double) r) == a)
				return r;
			++r;
		}

		return 0;
	}

	Commentator::Commentator () 
		: _estimationMethod (BEST_ESTIMATE), _format (OUTPUT_CONSOLE),
		  _show_timing (true), _show_progress (true), _show_est_time (true),
		  cnull (new nullstreambuf)
	{
		registerMessageClass (BRIEF_REPORT,         cout);
		registerMessageClass (PROGRESS_REPORT,      _report);
		registerMessageClass (TIMING_MEASURE,       _report);
		registerMessageClass (TIMING_ESTIMATE,      _report);
		registerMessageClass (PARTIAL_RESULT,       _report);
		registerMessageClass (INTERNAL_WARNING,     _report, 10, 3);
		registerMessageClass (INTERNAL_ERROR,       _report, 10, 3);
		registerMessageClass (INTERNAL_DESCRIPTION, _report);
	}

	void Commentator::start (const char *description, const char *fn = (const char *) 0, unsigned long len = 0) 
	{
		if (fn == (const char *) 0 && _activities.size () > 0)
			fn = _activities.top ()->_fn;

		report (LEVEL_IMPORTANT, INTERNAL_DESCRIPTION) << "Starting activity: " << description << endl;

		Activity *new_act = new Activity (description, fn, len);
		_activities.push (new_act);

		if (isPrinted (_activities.size (), LEVEL_IMPORTANT, BRIEF_REPORT, fn))
			printActivityReport (*new_act);

		new_act->_timer.start ();
	}

	void Commentator::startIteration (unsigned int iter, unsigned long len = 0) 
	{
		char buf[80];

		ostrstream str (buf, 80);

		str << "Iteration " << iter << ends;

		start (buf, (const char *) 0, len);
	}

	void Commentator::stop (const char *msg, const char *long_msg = (const char *) 0, const char *fn = (const char *) 0) 
	{
		float usertime;

		linbox_check (_activities.top () != (Activity *) 0);
		linbox_check (msg != (const char *) 0);

		if (long_msg == (const char *) 0)
			long_msg = msg;

		_activities.top ()->_timer.stop ();

		if (fn != (const char *) 0 &&
		    _activities.size () > 0 &&
		    _activities.top ()->_fn != (const char *) 0 &&
		    strcmp (fn, _activities.top ()->_fn) != 0)
		{
			report (LEVEL_IMPORTANT, INTERNAL_WARNING)
				<< "Activity report mismatch. Check that start () and stop () calls are paired correctly." << endl;
		}

		if (isPrinted (_activities.size (), LEVEL_IMPORTANT, BRIEF_REPORT, fn))
			finishActivityReport (*(_activities.top ()), msg);

		usertime = _activities.top ()->_timer.usertime ();

		delete _activities.top ();
		_activities.pop ();
		ostream &output = report (LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		output.precision (3);
		output << "Finished activity (" << long_msg << "): " << usertime << endl;
	}

	void Commentator::progress (long k = -1, long len = -1) 
	{
		linbox_check (_activities.top () != (Activity *) 0);

		if (k == -1)
			_activities.top ()->_progress++;
		else
			_activities.top ()->_progress = k;

		if (len != -1)
			_activities.top ()->_len = len;

		if (_activities.top ()->_progress > _activities.top ()->_len)
			_activities.top ()->_len = _activities.top ()->_progress;

		report (LEVEL_IMPORTANT, PROGRESS_REPORT)
			<< "Progress: " << _activities.top ()->_progress << " out of " << _activities.top ()->_len << endl;

		if (_show_progress && isPrinted (_activities.size (), LEVEL_IMPORTANT, BRIEF_REPORT, _activities.top ()->_fn))
			updateActivityReport (*(_activities.top ()));
	}

	ostream &Commentator::report (long level, const char *msg_class) 
	{
		linbox_check (msg_class != (const char *) 0);

		if (!isPrinted (_activities.size (), level, msg_class,
				(_activities.size () > 0) ? _activities.top ()->_fn : (const char *) 0))
			return cnull;

		MessageClass &messageClass = getMessageClass (msg_class);

		indent (messageClass._stream);

		return messageClass._stream;
	}

	void Commentator::indent (ostream &stream) 
	{
		unsigned int i;

		for (i = 0; i < _activities.size (); i++)
			stream << "  ";
	}

	void Commentator::setMaxDepth (long depth) 
	{
		map <const char *, MessageClass *, LessThanString>::iterator i;

		for (i = _messageClasses.begin (); i != _messageClasses.end (); i++)
			(*i).second->setMaxDepth (depth);
	}

	void Commentator::setMaxDetailLevel (long level) 
	{
		map <const char *, MessageClass *, LessThanString>::iterator i;

		for (i = _messageClasses.begin (); i != _messageClasses.end (); i++)
			(*i).second->setMaxDetailLevel (level);
	}

	MessageClass &Commentator::registerMessageClass (const char *msg_class, ostream &stream, unsigned long max_depth = 1, unsigned long max_level = 2) 
	{
		linbox_check (msg_class != (const char *) 0);

		MessageClass *new_obj = new MessageClass (msg_class, stream, max_depth, max_level);
		_messageClasses[msg_class] = new_obj;
		return *new_obj;
	}

	MessageClass &Commentator::cloneMessageClass (const char *new_msg_class, const char *msg_class) 
	{
		linbox_check (new_msg_class != (const char *) 0);
		linbox_check (msg_class != (const char *) 0);

		MessageClass *new_obj = new MessageClass (getMessageClass (msg_class));
		new_obj->_msg_class = new_msg_class;
		_messageClasses[new_msg_class] = new_obj;
		return *new_obj;
	}

	MessageClass &Commentator::cloneMessageClass (const char *new_msg_class, const char *msg_class, ostream &stream) 
	{
		linbox_check (new_msg_class != (const char *) 0);
		linbox_check (msg_class != (const char *) 0);

		MessageClass &old_obj = getMessageClass (msg_class);
		MessageClass *new_obj = new MessageClass (new_msg_class, stream, old_obj._configuration);
		_messageClasses[new_msg_class] = new_obj;
		return *new_obj;
	}

	MessageClass &Commentator::getMessageClass (const char *msg_class)
		{ return *_messageClasses[msg_class]; }

	void Commentator::setPrintParameters (long low_depth, long high_depth, long max_level, const char *fn = (const char *) 0) 
	{
		map <const char *, MessageClass *, LessThanString>::iterator i;

		for (i = _messageClasses.begin (); i != _messageClasses.end (); i++)
			(*i).second->setPrintParameters (low_depth, high_depth, max_level, fn);
	}

	void Commentator::setBriefReportParameters (OutputFormat format, bool show_timing, bool show_progress, bool show_est_time) 
	{
		_format        = format;
		_show_timing   = show_timing;
		_show_progress = show_progress;
		_show_est_time = show_est_time;
	}

	bool Commentator::isPrinted (long depth, long level, const char *msg_class, const char *fn = (const char *) 0)
	{
		if (_messageClasses.find (msg_class) == _messageClasses.end ())
			return false;

		MessageClass &messageClass = getMessageClass (msg_class);

		return messageClass.isPrinted (depth, level, fn);
	}

	void Commentator::setBriefReportStream (ostream &stream) 
		{ setMessageClassStream (BRIEF_REPORT, stream); }

	void Commentator::setReportStream (ostream &stream) 
	{
		setMessageClassStream (PROGRESS_REPORT,      stream);
		setMessageClassStream (TIMING_MEASURE,       stream);
		setMessageClassStream (TIMING_ESTIMATE,      stream);
		setMessageClassStream (PARTIAL_RESULT,       stream);
		setMessageClassStream (INTERNAL_ERROR,       stream);
		setMessageClassStream (INTERNAL_WARNING,     stream);
		setMessageClassStream (INTERNAL_DESCRIPTION, stream);
	}

	void Commentator::setMessageClassStream (const char *msg_class, ostream &stream) 
	{
		MessageClass *old_msg_class = _messageClasses[msg_class];
		cloneMessageClass (msg_class, msg_class, stream);
		delete old_msg_class;

	void Commentator::setDefaultReportFile (const char *filename) 
	{
		_report.open (filename);
	}

	void Commentator::printActivityReport (Activity &activity) 
	{
		MessageClass &messageClass = getMessageClass (BRIEF_REPORT);

		unsigned int i;
		if (_format == OUTPUT_CONSOLE) {
			messageClass._stream << activity._desc << "...";
			for (i = 0; i < _activities.size () - 1; i++)
				messageClass._stream << "  ";


			if (messageClass.isPrinted (_activities.size () + 1, LEVEL_IMPORTANT, activity._fn))
				messageClass._stream << endl;
			else if (_show_progress && activity._len > 0) {
			else if (_show_progress)
				messageClass._stream << "    0%";
			messageClass._smart_streambuf.stream ().flush ();
			messageClass._stream.flush ();
		else if (_format == OUTPUT_PIPE &&
		else if (_format == OUTPUT_PIPE) {
			for (i = 0; i < _activities.size () - 1; i++)
				messageClass._stream << "  ";


			if (_show_progress)
				messageClass._stream << endl;
		}
	}

	void Commentator::updateActivityReport (Activity &activity) 
	{
		MessageClass &messageClass = getMessageClass (BRIEF_REPORT);
		unsigned int i, old_len;
		unsigned int i;

		if (_format == OUTPUT_CONSOLE) {
			if (!messageClass.isPrinted (_activities.size (), LEVEL_IMPORTANT, activity._fn)) {
			if (!messageClass.isPrinted (_activities.size () + 1, LEVEL_IMPORTANT, activity._fn)) {
					for (i = 0; i < _last_line_len; i++)
					for (i = 0; i < 4; i++)
					str.width (3);
					messageClass._stream.width (3);
					messageClass._stream.precision (0);
					messageClass._stream << percent << '%';
			}
			else if (messageClass.isPrinted (_activities.size () - 1, LEVEL_UNIMPORTANT, activity._fn)) {
			else if (messageClass.isPrinted (_activities.size (), LEVEL_UNIMPORTANT, activity._fn)) {
				if (_show_progress) {
					for (i = 0; i < _activities.size () - 1; i++)
						messageClass._stream << "  ";

					messageClass._stream.precision (0);
					messageClass._stream << percent << "% done";
				if (_show_est_time)
					if (_show_est_time)
						messageClass._stream << " (" << activity._estimate.front ()._time
								     << " remaining)";
#endif
					messageClass._stream << endl;
				}
#if 0
				else if (_show_est_time) {
					for (i = 0; i < _activities.size () - 1; i++)
						messageClass._stream << "  ";

							     << " remaining" << endl;
#endif
				}
			}

			messageClass._smart_streambuf.stream ().flush ();
			messageClass._stream.flush ();
		else if (_format == OUTPUT_PIPE) {
			if (_show_progress) {
				messageClass._stream << floor (percent + 0.5) << "% done";
				for (i = 0; i < _activities.size () - 1; i++)
					messageClass._stream << "  ";

				messageClass._stream.precision (2);
				messageClass._stream << percent << "% done";
				if (_show_est_time)
					messageClass._stream << " (" << activity._estimate.front ()._time
							     << " remaining)";
#endif
				messageClass._stream << endl;
			}
#if 0
			else if (_show_est_time)
			else if (_show_est_time) {
				for (i = 0; i < _activities.size () - 1; i++)
					messageClass._stream << "  ";

						     << " remaining" << endl;
#endif
			}
		}
	}

	void Commentator::finishActivityReport (Activity &activity, const char *msg) 
	{
		MessageClass &messageClass = getMessageClass (BRIEF_REPORT);
		unsigned int i;

		if (_format == OUTPUT_CONSOLE) {
			if (!messageClass.isPrinted (_activities.size () + 1, LEVEL_IMPORTANT, activity._fn)) {
				if (_show_progress)
					for (i = 0; i < _last_line_len; i++)
					for (i = 0; i < 6; i++)

				messageClass._stream << msg;
				messageClass._stream << msg << endl;
			else if (messageClass.isPrinted (_activities.size (), LEVEL_UNIMPORTANT, activity._fn)) {
				for (i = 0; i < _activities.size (); i++)
				for (i = 0; i < _activities.size () - 1; i++)

				messageClass._stream << "Done: " << msg;
				messageClass._stream << "Done: " << msg << endl;

			messageClass._smart_streambuf.stream ().flush ();
			messageClass._stream.flush ();
		else if (_format == OUTPUT_PIPE) {
			for (i = 0; i < _activities.size (); i++)
			for (i = 0; i < _activities.size () - 1; i++)

			if (((_show_progress || _show_est_time) && activity._len > 0) ||
			messageClass._stream << "Done: " << msg << endl;
	}

	MessageClass::MessageClass (const Commentator &comm,
	MessageClass::MessageClass (const char *msg_class, ostream &stream, unsigned long max_depth = 1, unsigned long max_level = 2) 
		: _msg_class (msg_class), _stream (stream), _max_level (max_level), _max_depth (max_depth)
		fixDefaultConfig ();
	}

	void MessageClass::setMaxDepth (long depth) 
	{
		_max_depth = (unsigned long) depth;
		fixDefaultConfig ();
	}

	void MessageClass::setMaxDetailLevel (long level)
	{
		_max_level = (unsigned long) level;
		fixDefaultConfig ();
	}

	void MessageClass::setPrintParameters (unsigned long depth, unsigned long level, const char *fn) 
	void MessageClass::setPrintParameters (long low_depth, long high_depth, long max_level, const char *fn) 
		if (fn == (const char *) 0)
			fn = "";

		list <pair <unsigned long, unsigned long> > &config = _configuration[fn];
		list <pair <unsigned long, unsigned long> >::iterator i, j;
		list <pair <unsigned long, unsigned long> >::iterator i, j, k;

		long save_level;
		i = config.begin ();

		// Iterate through preceeding elements in the list and remove
		while (i != config.end () && (*i).first < (unsigned long) low_depth) i++;
		j = i;
		// Insert our new directive into the list
		if (i != config.end ())
			save_level = (*i).second;
		else
			save_level = -1;
		// Iterate through following elements in the list and remove any
		while (j != config.end () && (*j).first < (unsigned long) high_depth) {
			k = j;
			j++;
			config.erase (k);

		// End result: The list should be monotonically increasing in
		config.insert (i, pair <unsigned long, unsigned long> (low_depth, max_level));

		if (!(j != config.end () && (*j).first == (unsigned long) high_depth))
			config.insert (i, pair <unsigned long, unsigned long> (high_depth, save_level));

	bool MessageClass::isPrinted (unsigned long depth, unsigned long level, const char *fn)
	bool MessageClass::isPrinted (long depth, long level, const char *fn = (const char *) 0)
		if (checkConfig (_configuration[""], depth, level))
			return true;
		else if (fn != (const char *) 0)
			return checkConfig (_configuration[fn], depth, level);
		else
			return false;
	}

	MessageClass::MessageClass (const Commentator &comm,
	MessageClass::MessageClass (const char *msg_class, ostream &stream, Configuration configuration) 
		: _msg_class (msg_class), _stream (stream), _configuration (configuration)

	void MessageClass::fixDefaultConfig () 
	{
		list <pair <unsigned long, unsigned long> > &config = _configuration[""];


		while (config.size () > 0 && config.front ().first <= _max_depth)
			config.pop_front ();

		config.push_front (pair <unsigned long, unsigned long> (_max_depth, _max_level));

	bool MessageClass::checkConfig (list <pair <unsigned long, unsigned long> > &config, unsigned long depth, unsigned long level) 
	bool MessageClass::checkConfig (list <pair <unsigned long, unsigned long> > &config, long depth, long level) 
		list <pair <unsigned long, unsigned long> >::iterator i;

		i = config.begin ();
		while (i != config.end ()) {
			if (depth < i->first) {
			if ((unsigned long) depth <= (*i).first && (unsigned long) level <= (*i).second) {
				return true;
			} else if ((unsigned long) depth <= (*i).first) {
				return false;

		}

		return false;
	}

	void MessageClass::dumpConfig () const
	{
		Configuration::const_iterator i;
		list <pair <unsigned long, unsigned long> >::const_iterator j;

		for (i = _configuration.begin (); i != _configuration.end (); i++) {
			cerr << "Configuration (" << (*i).first << "):" << endl;

			for (j = (*i).second.begin (); j != (*i).second.end (); j++)
				cerr << "  Depth: " << (*j).first << ", Level: " << (*j).second << endl;

			cerr << endl;
		}
	}

	// Default global commentator
	Commentator commentator;
}
