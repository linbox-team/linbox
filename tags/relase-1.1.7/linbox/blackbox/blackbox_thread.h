/* Copyright (C) 2010 LinBox
 * Written by 
 * zhendong wan 
 *
 *
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
 */

#ifndef __LINBOX_blackbox_thread_H
#define __LINBOX_blackbox_thread_H

/* create a thread, which is bound to  a lwp to run matrix apply
 */

#include <linbox/vector/subvector.h>
#include <pthread.h>
#include <signal.h>
#include <string.h>

namespace LinBox 
{

 	/** built on posix threads
	\ingroup blackbox

 	 This is a thread interface, built on posix threads.
	  */
	class Thread {
	
	public:

#ifdef DEBUG
		static int count;
#endif
	
		// attributes for posix thread
        	static pthread_attr_t attr;

		// a signal set
        	static sigset_t sigset; 

        	static const int SIGAPPLY;

		static void terminate_thread (const Thread* t) {
		
			pthread_cancel(t -> getpid());

			pthread_join (t -> getpid(), 0);

		}

		// signal a thread
		static void signal_thread (const Thread* t) {

			pthread_kill ( t -> getpid(), SIGAPPLY);

		}

		// wait for a signal
		static void wait_thread (const Thread* t) {

			int err = 0, sig = 0;

			err = sigwait (&sigset, &sig);

			if (err == -1) {

				std::cerr << strerror (errno);


				exit (1);
			}
			
		}

		// notify the caller
		static void notify_parent (const Thread* t) {

			pthread_kill ( t -> getppid(), SIGAPPLY);

		}

		// terminate a thread
		static void terminate_thread (const Thread& t) {
		
			pthread_cancel(t. getpid());

			pthread_join (t. getpid(), 0);

		}

		static void signal_thread (const Thread& t) {

			pthread_kill ( t. getpid(), SIGAPPLY);

		}

		static void wait_thread (const Thread& t) {

			int err = 0, sig = 0;

			err = sigwait (&sigset, &sig);

			if (err == -1) {

				std::cerr << strerror (errno);


				exit (1);
			}
			
		}

		static void notify_parent (const Thread& t) {

			pthread_kill ( t. getppid(), SIGAPPLY);

		}

	protected:

		pthread_t pid;

		pthread_t ppid;
		
	public:	

		/** return the unique id associate with the thread.*/
		pthread_t getpid ( ) const {
		
			return pid;
		}
		
		/** set the unique id associate with the thread.*/
		void setpid (pthread_t _pid) {
			
			pid = _pid;
		}
		
		/** return caller's id */
		pthread_t getppid ( ) const {

                        return ppid;
                }

		/** set the caller's id */
                void setppid (pthread_t _ppid) {

                        ppid = _ppid;
                }

		/** run the thread */
		virtual void run (void) = 0;
		
		virtual ~Thread () {
#ifdef DEBUG
			--count;
#endif

		}

		Thread () {
#ifdef DEBUG

			++ count;
#endif
		}

	};

#ifdef DEBUG

	int Thread::count  = 0;
#endif

	/* initialization of static member of Thread */
	const int Thread::SIGAPPLY = SIGRTMIN;

	pthread_attr_t Thread::attr = (

			pthread_attr_init (&(Thread::attr)),

			pthread_attr_setscope (&(Thread::attr), PTHREAD_SCOPE_SYSTEM),

			Thread::attr
	);

	sigset_t Thread::sigset = ( 

			sigemptyset (&(Thread::sigset)),

			sigaddset (&(Thread::sigset), Thread::SIGAPPLY),

			(Thread::sigset)
	);

	
	/* a posix thread bound to an object of Thread */
	void* runThread (void* arg) {
		
		Thread* t = (Thread*) arg;

		pthread_setcancelstate (PTHREAD_CANCEL_ENABLE, NULL);

		pthread_setcanceltype (PTHREAD_CANCEL_ASYNCHRONOUS, NULL);

		t -> run ( );

		return 0;
	}

	// blackbox base
	// a pointer to the address of BB
	// a pointer to the address of output
	// a pointer to the address of input
	struct BBBase : public virtual Thread {

		typedef enum {Apply, ApplyTranspose} BBType;

		BBBase (const void* addr = NULL, void* outaddr = NULL, 
			const void* inaddr = NULL, BBType type = Apply) :

			Thread (),
		        bb_addr (addr), 
			bb_outaddr (outaddr), 
			bb_inaddr (inaddr),
			bb_type(type){
			

#ifdef DEBUG
				std::cout << "A blackbox was created with info:\n";

				std::cout << "\tAddress of thread " << bb_addr << '\n';

				std::cout << "\tAddrees of output " << bb_outaddr << '\n';
			
				std::cout << "\tAddress of input " << bb_inaddr << '\n';

#endif
		}

		virtual ~BBBase () {
		}

		const void* bb_addr;

		void*  bb_outaddr;

		const void*  bb_inaddr;

		BBType bb_type;

        };

        class LessTypeInfo {

		public:

                bool operator() (const std::type_info* info1, const std::type_info* info2) const {

                        return info1 -> before (*info2);
                }
        };

        typedef std::vector<BBBase*> BB_list;

        typedef std::map <const std::type_info*, BB_list, LessTypeInfo> BB_list_list;

	template <class Matrix, class Out, class In>

	class BBThread : public virtual BBBase {
		
		public:

		typedef Subvector<typename Out::iterator> SubOut;

		typedef Subvector<typename In::const_iterator> SubIn;
		
		BBThread(const Matrix* _m = NULL, Out* _outaddr = NULL,
			 const In* _inaddr = NULL, BBBase::BBType type= BBBase::Apply) :
			 Thread(), BBBase ((const void*)_m, (void*)_outaddr, (const void*)_inaddr, type) {

			
#ifdef DEBUG
				std::cout << "A BB thread was created with properties:\n";

				std::cout << "Matrix: ";

				_m -> write (std::cout);

				std::cout << '\n';

				std::cout << "Out address " << bb_outaddr << '\n';

				std::cout << "In address " << bb_inaddr << '\n';

#endif
		}
		
		inline void run (void) {

			const Matrix& matrix = *((const Matrix*)bb_addr);

			while (true) {

				//wait for signal
				Thread::wait_thread (this);

				switch (bb_type) {

					// Apply request
					case BBBase::Apply :

						matrix. apply ( *(SubOut*)bb_outaddr, 
						   	        *(const In*)bb_inaddr );

						break;

					// ApplyTranspose request
					case BBBase::ApplyTranspose :

						matrix. applyTranspose (*(std::vector<typename Matrix::Field::Element>*)bb_outaddr,
									*(const SubIn*)bb_inaddr);
					
						break;

					// othre choices
					default :

						std::cerr << "Unknown command type\n";

						break;
				}

				Thread::notify_parent (this);
			}
		}
		
		virtual ~BBThread () {

		delete (const Matrix*)bb_addr;

		}
	};

	/*- create posix thread corresponding a BBThread, and run it */
	template<class Matrix, class Out, class In>

		BBThread<Matrix, Out, In>* createBBThread ( const Matrix* m, Out* out, const In* in) {
		
		pthread_t id;
		
		BBThread<Matrix, Out, In>* t = new BBThread<Matrix, Out, In>(m, out, in);
		
		int ret = pthread_create (&id, &(Thread::attr), runThread, t);
		
		if (ret != 0) {
#ifdef DEBUG
			std::cout << "Number of thread created: " << Thread::count << '\n';
#endif

			std::cout << "Not enough resource\n";

			std::cerr << strerror(ret);

			exit (1);
		}
		
		t -> setpid (id);

		t -> setppid (pthread_self ());
		
		return t;
	}
	
}
#endif //__LINBOX_blackbox_thread_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
