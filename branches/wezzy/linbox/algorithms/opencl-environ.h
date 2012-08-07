/* linbox/algorithms/opencl-environ.h
 * Copyright (C) 2012 Matthew Wezowicz
 *
 * Written by Matthew Wezowicz <mwezz@udel.edu>
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

#include <cstdlib>
#include <pthread.h>

#include <algorithm>
#include <istream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

#include <CL/cl.h>

#include "named-mutex.h"

#ifndef __LINBOX_opencl_environ_H
#define __LINBOX_opencl_environ_H

namespace LinBox{

	/**
	 * Container for all pertenant information needed to use an OpenCL device,
	 * compile kernels for the device, track resource usage, and gain exclusive
	 * access to the device
	 *
	 * Complies to the OpenCL 1.0 spec
	 */
	class OpenCLEnviron{

	private:
		//Platform info
		std::string platformName;
		double platformVersion;
		std::vector<std::string> platformExtensions;

		//Device info
		std::vector<std::string> deviceExtensions;

		//OpenCL variables
		cl_context context;
		cl_device_id device;
		cl_command_queue commandQueue;
		cl_int errcode;

		//Memory stats
		unsigned long memCapacity;
		unsigned long maxBufferSize;
		unsigned long memInUse;

		//Device type flags
		bool CPU;
		bool GPU;
		bool ACCEL;
		bool Other;

		//Lock to control access to environ
#ifndef __MPI_SHARED
		pthread_mutex_t* deviceLock;
#else
		NamedMutex* deviceLock;
#endif

		//ID number
		unsigned int IDnum;

		OpenCLEnviron();
		OpenCLEnviron(const OpenCLEnviron&);
		OpenCLEnviron& operator=(const OpenCLEnviron&);

	public:
		OpenCLEnviron(
			std::string platformName,
			double platformVersion,
			std::vector<std::string> platformExtensions,
			cl_context context,
			cl_device_id device,
#ifndef __MPI_SHARED
			pthread_mutex_t* deviceLock,
#else
			NamedMutex* deviceLock,
#endif
			unsigned int IDnum) :
			//Init list starts here
			platformName(platformName),
			platformVersion(platformVersion),
			platformExtensions(platformExtensions),
			context(context),
			device(device),
			errcode(CL_SUCCESS),
			memCapacity(0L),
			maxBufferSize(0L),
			memInUse(0L),
			CPU(false),
			GPU(false),
			ACCEL(false),
			Other(true),
			deviceLock(deviceLock),
			IDnum(IDnum){

			//Initialize command queue
			clCreateCommandQueue(
				context,
				device,
				CL_QUEUE_PROFILING_ENABLE,
				&errcode);

			//Get memory stats
			errcode = clGetDeviceInfo(
				device,
				CL_DEVICE_GLOBAL_MEM_SIZE,
				sizeof(cl_ulong),
				&memCapacity,
				NULL);
			errcode = clGetDeviceInfo(
				device,
				CL_DEVICE_MAX_MEM_ALLOC_SIZE,
				sizeof(cl_ulong),
				&maxBufferSize,
				NULL);

			//Get device type
			cl_device_type type;

			errcode = clGetDeviceInfo(
				device,
				CL_DEVICE_TYPE,
				sizeof(cl_device_type),
				&type,
				NULL);

			CPU = (type == CL_DEVICE_TYPE_CPU);
			GPU = (type == CL_DEVICE_TYPE_GPU);
			ACCEL = (type == CL_DEVICE_TYPE_ACCELERATOR);
			Other = ((CPU || GPU || ACCEL) == false);

			//Enumerate device extensions
			size_t sizeRet;

			errcode = clGetDeviceInfo(
				device,
				CL_DEVICE_EXTENSIONS,
				0,
				NULL,
				&sizeRet);

			char* devExt = (char*)malloc(sizeRet);

			errcode = clGetDeviceInfo(
				device,
				CL_DEVICE_EXTENSIONS,
				sizeRet,
				devExt,
				NULL);

			//Parse out individual extensions into vector
			std::string devExt2(devExt);
			std::stringstream strstream(devExt2);
			std::istream_iterator<std::string> it(strstream);
			std::istream_iterator<std::string> end;

			for(; it != end; it++){
				deviceExtensions.push_back(*it);
			}

			//Deallocate memory
			free(devExt);
		}

		~OpenCLEnviron(){}

		/**
		 * getPlatformName() returns string with platform name
		 */
		std::string getPlatformName() const{
			return platformName;
		}

		/**
		 * getPlatformVersion() returns double with platform version
		 */
		double getPlatformVersion() const{
			return platformVersion;
		}

		/**
		 * getPlatformExtensions() returns a constant reference to the vector with
		 * the platform extensions supported by this environ
		 */
		const std::vector<std::string>& getPlatformExtensions() const{
			return platformExtensions;
		}

		/**
		 * getDeviceExtensions() returns a constant reference to the vector with
		 * the device extensions supported by this environ
		 */
		const std::vector<std::string>& getDeviceExtensions() const{
			return deviceExtensions;
		}

		/**
		 * getIDNum() returns the unsigned int with the OpenCLEnviron ID number
		 */
		unsigned int getIDNum() const{
			return IDnum;
		}

		/**
		 * getContext() accessor for the context associated with environ
		 */
		cl_context getContext() const{
			return context;
		}

		/**
		 * getDevice() accessor for the device associated with environ
		 */
		cl_device_id getDevice() const{
			return device;
		}

		/**
		 * getCommandQueue() accessor for thecommand queue associated with environ
		 */
		cl_command_queue getCommandQueue() const{
			return commandQueue;
		}

		/**
		 * getErrorCode() accessor for the error code associated with environ
		 */
		cl_int getErrorCode() const{
			return errcode;
		}

		/**
		 * setErrorCode() accessor for setting error code associated with environ
		 */
		void setErrorCode(cl_int err){
			errcode = err;
		}

		/**
		 * getMemCapacity() accessor for the memory capacity associated with environ
		 */
		unsigned long getMemCapacity() const{
			return memCapacity;
		}

		/**
		 * getMaxBufferSize() accessor for the maximum buffer size associated
		 * with environ
		 */
		unsigned long getMaxBufferSize() const{
			return maxBufferSize;
		}

		/**
		 * getMemInUse() accessor for the memory currently in use by the environ
		 */
		unsigned long getMemInUse() const{
			return memInUse;
		}

		/**
		 * getMemAvailable() accessor for the memory available in the environ
		 */
		unsigned long getMemAvailable() const{
			return memCapacity - memInUse;
		}

		/**
		 * memAllocated() setter for updating the memory levels
		 */
		void memAllocated(unsigned long alloc){
			memInUse += alloc;
		}

		/**
		 * memDeallocated() setter for updating the memory levels
		 */
		void memDeallocated(unsigned long dealloc){
			memInUse -= dealloc;
		}

		/**
		 * isCPU() accessor for the device type in the environ
		 */
		bool isCPU() const{
			return CPU;
		}

		/**
		 * isGPU() accessor for the device type in the environ
		 */
		bool isGPU() const{
			return GPU;
		}

		/**
		 * isAccelerator() accessor for the device type in the environ
		 */
		bool isAccelerator() const{
			return ACCEL;
		}

		/**
		 * isOtherDevice() accessor for the device type in the environ
		 */
		bool isOtherDevice() const{
			return Other;
		}

		/**
		 * acquireEnviron() accessor for gaining exclusive access to the environ
		 */
		void acquireEnviron(){
#ifndef __MPI_SHARED
			pthread_mutex_lock(deviceLock);
#else
			deviceLock->lock();
			memInUse = deviceLock->updateLocalValue<unsigned long>("memInUse");
#endif
		}

		/**
		 * releaseEnviron() accessor for releasing exclusive access to the environ
		 */
		void releaseEnviron(){
#ifndef __MPI_SHARED
			pthread_mutex_unlock(deviceLock);
#else
			deviceLock->updateGlobalValue<unsigned long>("memInUse", memInUse);
			deviceLock->unlock();
#endif
		}

		/**
		 * getDeviceLock() accessor for underlying lock in the environ
		 */
#ifndef __MPI_SHARED
		pthread_mutex_t* getDeviceLock(){
#else
		NamedMutex* getDeviceLock(){
#endif
			return deviceLock;
		}
	};

} /* end of namespace LinBox */

#endif /* __LINBOX_opencl_environ_H */
