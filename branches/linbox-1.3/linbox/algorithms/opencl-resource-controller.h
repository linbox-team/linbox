/* linbox/algorithms/opencl-resource-controller.h
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

#include <pthread.h>

#include <algorithm>
#include <istream>
#include <iterator>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <CL/cl.h>

//#include "named-mutex.h"
#include "opencl-environ.h"

#ifndef __LINBOX_opencl_resource_controller_H
#define __LINBOX_opencl_resource_controller_H

namespace LinBox{

	/**
	 * Enumerate all of the platforms currently available on the system
	 */
	std::vector<cl_platform_id> enumPlatforms(){
		//Variables
		cl_int errcode;
		cl_platform_id* platforms;
		cl_uint numPlatforms;

		//Get number of platforms on system
		errcode = clGetPlatformIDs(0, NULL, &numPlatforms);

		//Return empty vector if error
		if(errcode != CL_SUCCESS || numPlatforms <= 0){
			std::vector<cl_platform_id> ret;
			return ret;
		}

		//Allocate memory for platforms IDs
		platforms = new cl_platform_id[numPlatforms];

		//Read in platform IDs
		errcode = clGetPlatformIDs(numPlatforms, platforms, NULL);

		//Allocate vector for platform IDs
		std::vector<cl_platform_id> ret;

		//Copy platform IDs into vector
		for(int i = 0; i < (int)numPlatforms; i++){
			ret.push_back(platforms[i]);
		}

		//Deallocate memory
		delete[] platforms;

		return ret;
	}

	/**
	 * Get the platform name associated with the platform
	 */
	std::string getPlatformName(cl_platform_id platform){
		//Variables
		char* tempBuffer;
		cl_int errcode;
		size_t sizeRet;

		//Get size of platform name
		errcode = clGetPlatformInfo(platform, CL_PLATFORM_NAME, 0, NULL, &sizeRet);

		//Allocate memory for platform name
		tempBuffer = new char[sizeRet];

		//Read in platform name
		errcode = clGetPlatformInfo(
			platform,
			CL_PLATFORM_NAME,
			sizeRet,
			tempBuffer,
			NULL);

		//Convert platform name to std::string
		std::string ret(tempBuffer);

		//Deallocate memory
		delete[] tempBuffer;

		return ret;
	}

	/**
	 * Get the platform version associated with the platform
	 */
	double getPlatformVersion(cl_platform_id platform){
		//Variables
		char* tempBuffer;
		cl_int errcode;
		double ret;
		size_t first;
		size_t second;
		size_t sizeRet;

		//Get size of platform version string
		errcode = clGetPlatformInfo(
			platform,
			CL_PLATFORM_VERSION,
			0,
			NULL,
			&sizeRet);

		//Allocate memory for platform version string
		tempBuffer = new char[sizeRet];

		//Read in platform name string
		errcode = clGetPlatformInfo(
			platform,
			CL_PLATFORM_VERSION,
			sizeRet,
			tempBuffer,
			NULL);

		//Convert platform version string to std::string
		std::string versionString(tempBuffer);

		//Isolate the version number in string and convert to double
		first = versionString.find_first_of(" ", 0);
		second = versionString.find_first_of(" ", (first + 1));
		std::string tempString = versionString.substr(first, (second - first));
		std::stringstream strstream(tempString);
		strstream >> ret;

		//Deallocate memory
		delete[] tempBuffer;

		return ret;
	}

	/**
	 * Get the platform extensions associated with the platform
	 */
	std::vector<std::string> getPlatformExtensions(cl_platform_id platform){
		//Variables
		char* tempBuffer;
		cl_int errcode;
		size_t sizeRet;

		//Get size of platform extensions string
		errcode = clGetPlatformInfo(
			platform,
			CL_PLATFORM_EXTENSIONS,
			0,
			NULL,
			&sizeRet);

		//Allocate memory for platform extensions string
		tempBuffer = new char[sizeRet];

		//Read in platform extensions string
		errcode = clGetPlatformInfo(
			platform,
			CL_PLATFORM_EXTENSIONS,
			sizeRet,
			tempBuffer,
			NULL);

		//Convert platfrom extensions string to std::string
		std::string tempString(tempBuffer);

		//Parse out the individual extensions and copy into vector
		std::stringstream strstream(tempString);
		std::istream_iterator<std::string> it(strstream);
		std::istream_iterator<std::string> end;
		std::vector<std::string> ret;
		for(; it != end; it++){
			ret.push_back(*it);
		}

		//Deallocate memory
		delete[] tempBuffer;

		return ret;
	}

	/**
	 * Enumerate all of the devices currently available on the platform
	 */
	std::vector<cl_device_id> enumDevices(cl_platform_id platform){
		//Variables
		cl_device_id* devices;
		cl_int errcode;
		cl_uint numDevices;

		//Get number of devices in platform
		errcode = clGetDeviceIDs(
			platform,
			CL_DEVICE_TYPE_ALL,
			0,
			NULL,
			&numDevices);

		//Return empty vector if error
		if(errcode != CL_SUCCESS || numDevices <= 0){
			std::vector<cl_device_id> ret;
			return ret;
		}

		//Allocate memory for device IDs
		devices = new cl_device_id[numDevices];

		//Read in device IDs
		errcode = clGetDeviceIDs(
			platform,
			CL_DEVICE_TYPE_ALL,
			numDevices,
			devices,
			NULL);

		//Allocate vector for device IDs
		std::vector<cl_device_id> ret;

		//Copy device IDs into vector
		for(int i = 0; i < (int)numDevices; i++){
			ret.push_back(devices[i]);
		}

		//Deallocate memory
		delete[] devices;

		return ret;
	}

	/**
	 * Create an OpenCL context from a platfrom and device
	 */
	cl_context createContext(cl_platform_id platform, cl_device_id device){
		//Variables
		cl_context ret;
		cl_int errcode;
		cl_context_properties props[3] = {
			CL_CONTEXT_PLATFORM,
			(cl_context_properties)platform,
			0};

		//Create context
		ret = clCreateContext(props, 1, &device, NULL, NULL, &errcode);

		return ret;
	}

	class OpenCLResourceController{

	private:
		static std::vector<OpenCLEnviron*>* environs;

#ifndef __MPI_SHARED
		static std::map<unsigned int, unsigned int>* IDsToInstances;
		static pthread_mutex_t controllerLock;
#else
		static std::map<unsigned int, std::string>* IDsToNames;
		static std::map<std::string, unsigned int>* namesToInstances;
		static std::vector<NamedMutex*>* mutexs;
		static NamedMutex controllerLock;
#endif

		static bool initialized;

		void init(){
			unsigned int ID = 0;

			environs = new std::vector<OpenCLEnviron*>;
#ifndef __MPI_SHARED
			IDsToInstances = new std::map<unsigned int, unsigned int>;
#else
			mutexs = new std::vector<NamedMUtex*>;
			IDsToNames = new std::map<unsigned int, std::string>;
			namesToInstances = new std::map<std::string, unsigned int>;
#endif

			std::vector<cl_platform_id> platforms = enumPlatforms();

			for(int i = 0; i < (int)platforms.size(); i++){
				std::string platformName = getPlatformName(platforms[i]);

				double platformVersion = getPlatformVersion(platforms[i]);

				std::vector<std::string> platformExtensions =
					getPlatformExtensions(platforms[i]);

				std::vector<cl_device_id> devices = enumDevices(platforms[i]);

				for(int j = 0; j < (int)devices.size(); j++){
					cl_context context = createContext(platforms[i], devices[j]);

#ifndef __MPI_SHARED
					//Allocate and initialize device lock
					pthread_mutex_t* tempMutex = new pthread_mutex_t;

					pthread_mutex_init(tempMutex, NULL);

					//Build OpenCLEnviron
					OpenCLEnviron* tempEnviron = new OpenCLEnviron(
						platformName,
						platformVersion,
						platformExtensions,
						context,
						devices[j],
						tempMutex,
						ID);

					//Update map of ID numbers to allocated instances
					(*IDsToInstances)[ID] = 0;

					ID++;
#else
					//Build NamedMutex name
					std::stringstream strstream;
					strstream << platformName << "_" << j;
					std::string mutexName(strstream.str());

					//Construct named mutex
					NamedMutex* tempMutex = new NamedMutex(mutexName);
					mutexs->push_back(tempMutex);

					//Build OpenCLEnviron with local ID number
					OpenCLEnviron* tempEnviron = new OpenCLEnviron(
						platformName,
						platformVersion,
						platformExtensions,
						context,
						devices[j],
						tempMutex,
						ID);

					//Update local copy of map with
					//NamedMutex names to allocated instances
					unsigned int allocated =
						controllerLock.updateLocalValue<unsigned int>(mutexName);
					(*namesToInstances)[mutexName] = allocated;

					//Update local map of ID numbers to NamedMutex names
					(*IDsToNames)[ID] = mutexName;

					//Update global copy of map with
					//NamedMutex names to allocated instances
					controllerLock.updateGlobalValue<unsigned int>(mutexName, allocated);

					ID++;
#endif

					environs->push_back(tempEnviron);
				}
			}

			initialized = true;
		}

	public:

		OpenCLResourceController(){
#ifndef __MPI_SHARED
			pthread_mutex_lock(&controllerLock);

			init();
			//removeDuplicates();

			pthread_mutex_unlock(&controllerLock);
#else
			controllerLock.lock();

			init();
			//removeDuplicates();

			controllerLock.unlock();
#endif
		}

		~OpenCLResourceController(){
#ifndef __MPI_SHARED
			for(int i = 0; i < (int)environs->size(); i++){
				pthread_mutex_destroy(environs->at(i)->getDeviceLock());
				delete environs->at(i)->getDeviceLock();
				delete environs->at(i);
			}
#else
			for(int i = 0; i < mutexs->size(); i++){
				delete mutexs->at(i);
			}
			for(int i = 0; i < environs->size(); i++){
				delete environs->at(i);
			}
#endif
		}

		const std::vector<OpenCLEnviron*>* accessEnvirons(){
			return environs;
		}

		void acquireController(){
#ifndef __MPI_SHARED
			pthread_mutex_lock(&controllerLock);
#else
			controllerLock.lock();

			//Update local copy of map with
			//NamedMutex names to allocated instances
			///!TODO
#endif
		}

		void releaseController(){
#ifndef __MPI_SHARED
			pthread_mutex_unlock(&controllerLock);
#else
			//Update global copy of map with
			//NamedMutex names to allocated instances
			///!TODO

			controllerLock.unlock();
#endif
		}

		unsigned int getInstanceCount(unsigned int ID){
#ifndef __MPI_SHARED
			return (*IDsToInstances)[ID];
#else
			return (*namesToInstances)[(*IDsToNames)[ID]];
#endif
		}

		unsigned int allocateInstance(unsigned int ID){
#ifndef __MPI_SHARED
			(*IDsToInstances)[ID] += 1;
			return (*IDsToInstances)[ID];
#else
			(*namesToInstances)[(*IDsToNames)[ID]] += 1;
			return (*namesToInstances)[(*IDsToNames)[ID]];
#endif
		}

		unsigned int deallocateInstance(unsigned int ID){
#ifndef __MPI_SHARED
			(*IDsToInstances)[ID] -= 1;
			if((*IDsToInstances)[ID] < 0){ //BB Toujours faux
				(*IDsToInstances)[ID] = 0;
			}
			return (*IDsToInstances)[ID];
#else
			(*namesToInstances)[(*IDsToNames)[ID]] -= 1;
			if((*namesToInstances)[(*IDsToNames)[ID]] < 0){
				(*namesToInstances)[(*IDsToNames)[ID]] = 0;
			}
			return (*namesToInstances)[(*IDsToNames)[ID]];
#endif
		}

	}; /* end of class OpenCLResourceController */

	std::vector<OpenCLEnviron*>* OpenCLResourceController::environs = NULL;
#ifndef __MPI_SHARED
	std::map<unsigned int, unsigned int>*
		OpenCLResourceController::IDsToInstances = NULL;
	pthread_mutex_t OpenCLResourceController::controllerLock =
		PTHREAD_MUTEX_INITIALIZER;
#else
	std::map<unsigned int, std::string>*
		OpenCLResourceController::IDsToNames = NULL;
	std::map<std::string, unsigned int>*
		OpenCLResourceController::namesToInstances = NULL;
	std::vector<NamedMutex*>*
		OpenCLResourceController::mutexs = NULL;
	NamedMutex OpenCLResourceController::controllerLock =
		NamedMutex("OpenCLResourceContollerLock");
#endif
	bool OpenCLResourceController::initialized(false);

	OpenCLResourceController& accessOpenCLResourceController(){
		static OpenCLResourceController controller;
		return controller;
	}

} /* end of namespace LinBox */

#endif // __LINBOX_opencl_resource_controller_H

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

