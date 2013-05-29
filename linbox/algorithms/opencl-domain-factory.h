/* linbox/algorithms/opencl-domain-factory.h
 * Copyright (C) 2011-2012 Matthew Wezowicz
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

#ifndef __LINBOX_opencl_matrix_domain_factory_H
#define __LINBOX_opencl_matrix_domain_factory_H

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <pthread.h>

#include <vector>

#include "linbox/algorithms/opencl-domain.h"
#include "linbox/algorithms/opencl-environ.h"
#include "linbox/algorithms/opencl-resource-controller.h"
#include "linbox/algorithms/opencl-kernels/opencl-domain-kernels.inl"

#include "CL/cl.h"

#define NUM_KERNELS 20

namespace LinBox{

	const char* dpKernelSources[] = {
			matrixMulKernelModular1DP, matrixMulKernelModular8DP,
			matrixMulKernelModular32DP, matrixMulKernelModular1024DP,
			matrixMuladdKernelModular1DP, matrixMuladdKernelModular8DP,
			matrixMuladdKernelModular32DP, matrixMuladdKernelModular1024DP,
			matrixAxpyKernelModular1DP, matrixAxpyKernelModular8DP,
			matrixAxpyKernelModular32DP, matrixAxpyKernelModular1024DP,
			matrixMaxpyKernelModular1DP, matrixMaxpyKernelModular8DP,
			matrixMaxpyKernelModular32DP, matrixMaxpyKernelModular1024DP,
			matrixAxmyKernelModular1DP, matrixAxmyKernelModular8DP,
			matrixAxmyKernelModular32DP, matrixAxmyKernelModular1024DP};

		const char* dpKernelNames[] = {
			"matrixMulKernelModular1DP", "matrixMulKernelModular8DP",
			"matrixMulKernelModular32DP", "matrixMulKernelModular1024DP",
			"matrixMuladdKernelModular1DP", "matrixMuladdKernelModular8DP",
			"matrixMuladdKernelModular32DP", "matrixMuladdKernelModular1024DP",
			"matrixAxpyKernelModular1DP", "matrixAxpyKernelModular8DP",
			"matrixAxpyKernelModular32DP", "matrixAxpyKernelModular1024DP",
			"matrixMaxpyKernelModular1DP", "matrixMaxpyKernelModular8DP",
			"matrixMaxpyKernelModular32DP", "matrixMaxpyKernelModular1024DP",
			"matrixAxmyKernelModular1DP", "matrixAxmyKernelModular8DP",
			"matrixAxmyKernelModular32DP", "matrixAxmyKernelModular1024DP"};

		const char* spKernelSources[] = {
			matrixMulKernelModular1SP, matrixMulKernelModular16SP,
			matrixMulKernelModular32SP, matrixMulKernelModular1024SP,
			matrixMuladdKernelModular1SP, matrixMuladdKernelModular16SP,
			matrixMuladdKernelModular32SP, matrixMuladdKernelModular1024SP,
			matrixAxpyKernelModular1SP, matrixAxpyKernelModular16SP,
			matrixAxpyKernelModular32SP, matrixAxpyKernelModular1024SP,
			matrixMaxpyKernelModular1SP, matrixMaxpyKernelModular16SP,
			matrixMaxpyKernelModular32SP, matrixMaxpyKernelModular1024SP,
			matrixAxmyKernelModular1SP, matrixAxmyKernelModular16SP,
			matrixAxmyKernelModular32SP, matrixAxmyKernelModular1024SP};

		const char* spKernelNames[] = {
			"matrixMulKernelModular1SP", "matrixMulKernelModular16SP",
			"matrixMulKernelModular32SP", "matrixMulKernelModular1024SP",
			"matrixMuladdKernelModular1SP", "matrixMuladdKernelModular16SP",
			"matrixMuladdKernelModular32SP", "matrixMuladdKernelModular1024SP",
			"matrixAxpyKernelModular1SP", "matrixAxpyKernelModular16SP",
			"matrixAxpyKernelModular32SP", "matrixAxpyKernelModular1024SP",
			"matrixMaxpyKernelModular1SP", "matrixMaxpyKernelModular16SP",
			"matrixMaxpyKernelModular32SP", "matrixMaxpyKernelModular1024SP",
			"matrixAxmyKernelModular1SP", "matrixAxmyKernelModular16SP",
			"matrixAxmyKernelModular32SP", "matrixAxmyKernelModular1024SP"};

	class OpenCLMatrixDomainFactory{

	protected:
		/*--- This class is a Singleton ---*/
		//No one can construct the class
		OpenCLMatrixDomainFactory();
		//No one can copy the class
		OpenCLMatrixDomainFactory(const OpenCLMatrixDomainFactory&);
		//No one can assign the class
		OpenCLMatrixDomainFactory& operator=(const OpenCLMatrixDomainFactory&);
		//No one can destroy the class
		~OpenCLMatrixDomainFactory();

		struct oclEnviron{
			//OpenCL specific variables
			cl_context context;
			cl_device_id device;
			cl_command_queue commandQue;
			cl_int errcode;

			//Storage for memory levels
			unsigned long memCapacity;
			unsigned long maxBufferSize;

			//Type and status flags
			bool GPUcontainer;
			bool CPUcontainer;
			bool setupCorrect;
			bool doubleSupported;

			//Storage for kernels
			cl_kernel* dpKernels;
			bool* dpKernelsAvailable;
			cl_kernel* spKernels;
			bool* spKernelsAvailable;

			//Mutex
			pthread_mutex_t* deviceLock;
		};

		public:
		typedef struct oclEnviron oclEnviron;

		protected:

		//Storage for all OpenCL evironments created by the factory
		static std::vector<oclEnviron>* environs;
		static std::vector<int>* instances;

		//Initialization flag
		static bool initialized;

		//Storage for the error code
		static cl_int errcode;

		//Factory mutex
		static pthread_mutex_t factoryLock;

		/**
		 * @internal
		 * Picks the platform used
		 */
		static cl_int oclGetPlatformID(cl_platform_id& selectedPlatform){
			//Allocate temporary char array for platform name
			char chBuffer[256];
			cl_uint numPlatforms;
			cl_platform_id* platforms;
			selectedPlatform = NULL;

			//Get OpenCL platform count
			//return custom error codes if there are no platforms
			//or could not get number of platforms
			errcode = clGetPlatformIDs(0, NULL, &numPlatforms);
			if(errcode != CL_SUCCESS){
				return -1000;
			}
			if(numPlatforms == 0){
				return -2000;
			}

			//Allocate space to store cl_platform_id's
			platforms = (cl_platform_id*)malloc(numPlatforms * sizeof(cl_platform_id));
			errcode = clGetPlatformIDs(numPlatforms, platforms, NULL);

			//Search through the platforms looking for a prefered one specified by a
			//string
			for(int i = 0; i < (int)numPlatforms; i++){
				errcode = clGetPlatformInfo(
					platforms[i],
					CL_PLATFORM_NAME,
					256,
					(void*)chBuffer,
					NULL);

				if(errcode == CL_SUCCESS && !(strcmp(chBuffer,"NVIDIA CUDA"))){
					selectedPlatform = platforms[i];
					break;
				}
			}

			//If prefered platform could not be found use first platform
			if(selectedPlatform == NULL){
				selectedPlatform = platforms[0];
			}

			//Clean up memory
			free(platforms);

			return CL_SUCCESS;
		}

		/**
		 * @internal
		 * Computes a score for all devices
		 */
		static std::vector<long> oclComputeDeviceScores(
			cl_device_id* devices,
			cl_uint numDevices){

			std::vector<long> ret;

			//For all devices, query information about them and compute a score
			for(int i = 0; i < (int)numDevices; i++){

				cl_device_type type;
				errcode = clGetDeviceInfo(
					devices[i],
					CL_DEVICE_TYPE,
					sizeof(cl_device_type),
					&type,
					NULL);

				bool GPU = (type == CL_DEVICE_TYPE_GPU);

				cl_uint computeUnits;
				errcode = clGetDeviceInfo(
					devices[i],
					CL_DEVICE_MAX_COMPUTE_UNITS,
					sizeof(cl_uint),
					&computeUnits,
					NULL);

				cl_uint clockFrequency;
				errcode = clGetDeviceInfo(
					devices[i],
					CL_DEVICE_MAX_CLOCK_FREQUENCY,
					sizeof(cl_uint),
					&clockFrequency,
					NULL);

				cl_ulong maxGlobalMemoryAllocSize;
				errcode = clGetDeviceInfo(
					devices[i],
					CL_DEVICE_MAX_MEM_ALLOC_SIZE,
					sizeof(cl_ulong),
					&maxGlobalMemoryAllocSize,
					NULL);

				cl_ulong globalMemory;
				errcode = clGetDeviceInfo(
					devices[i],
					CL_DEVICE_GLOBAL_MEM_SIZE,
					sizeof(cl_ulong),
					&globalMemory,
					NULL);

				long score = (computeUnits * clockFrequency);
				score *= (GPU ? 8 : 1);
				score += (maxGlobalMemoryAllocSize / (1024 * 1024));
				score += (globalMemory / (1024 * 1024));

				ret.push_back(score);
			}

			return ret;
		}

		/**
		 * @internal
		 * Creates a kernel given a file name and kernel name
		 * Returns the kernel
		 */
		static cl_kernel oclCreateKernel(
			const char* kernel,
			const char* kernelName,
			cl_context& context){

			//find length of the kernel
			size_t kernelLength = strlen(kernel);

			//Create program from file
			cl_program program = clCreateProgramWithSource(
				context,
				1,
				(const char**)&kernel,
				&kernelLength,
				&errcode);

			//Build the program into executable
			if(errcode == CL_SUCCESS){
				const char* options = {"-cl-mad-enable -cl-no-signed-zeros "
				                       "-cl-finite-math-only\0"};
				errcode = clBuildProgram(program, 0, NULL, options, NULL, NULL);
			}

			//Create kernel from executable
			cl_kernel tempKernel = NULL;
			if(errcode == CL_SUCCESS){
				tempKernel = clCreateKernel(program, kernelName, &errcode);
			}

			//Releasing program
			clReleaseProgram(program);

			//Return kernel
			return tempKernel;
		}

		/**
		 * @internal
		 * Builds a single oclEnviron
		 */
		static oclEnviron& oclBuildEnviron(
			oclEnviron& environ,
			cl_platform_id& platform,
			cl_device_id& device){

			//Create OpenCL context
			cl_context_properties properties[3] = {CL_CONTEXT_PLATFORM,
			                                       (cl_context_properties)platform,
			                                       0};
			environ.context = clCreateContext(
				properties,
				1,
				&device,
				NULL,
				NULL,
				&errcode);

			//Assign device to environ
			environ.device = device;

			//Create OpenCL command queue
			environ.commandQue = clCreateCommandQueue(
				environ.context,
				device,
				0,
				&errcode);

			//Get amount of memory that the device has
			errcode = clGetDeviceInfo(
				device,
				CL_DEVICE_GLOBAL_MEM_SIZE,
				sizeof(cl_ulong),
				&(environ.memCapacity),
				NULL);
			errcode = clGetDeviceInfo(
				device,
				CL_DEVICE_MAX_MEM_ALLOC_SIZE,
				sizeof(cl_ulong),
				&(environ.maxBufferSize),
				NULL);

			//Get the device type for the selected device
			cl_device_type type;
			errcode = clGetDeviceInfo(
				device,
				CL_DEVICE_TYPE,
				sizeof(cl_device_type),
				&type,
				NULL);

			//Set container type
			environ.GPUcontainer = (type == CL_DEVICE_TYPE_GPU);
			environ.CPUcontainer = (type == CL_DEVICE_TYPE_CPU);


			//Determine double precision support for device by getting a
			//listing of all OpenCL device extensions supported by the device
			//and searching through them for the string "cl_khr_fp64"
			environ.doubleSupported = false;

			size_t sizeReturn;

			errcode = clGetDeviceInfo(
				device,
				CL_DEVICE_EXTENSIONS,
				0,
				NULL,
				&sizeReturn);
			char* deviceExtensions = (char*)malloc(sizeReturn);
			errcode = clGetDeviceInfo(
				device,
				CL_DEVICE_EXTENSIONS,
				sizeReturn,
				deviceExtensions,
				NULL);

			char* exten[200];

			exten[0] = strtok(deviceExtensions, " ");

			if(exten[0] != NULL){
				for(int i = 1; i < 200; i++){
					exten[i] = strtok(NULL, " ");
					if(exten[i] == NULL){
						break;
					}
				}

				for(int i = 0; i < 200; i++){
					if(exten[i] == NULL){
						break;
					}
					if(!strcmp(exten[i], "cl_khr_fp64")){
						environ.doubleSupported = true;
					}
				}
			}

			free(deviceExtensions);

			//Copy current errcode to environ as starting errcode
			environ.errcode = errcode;
			if(environ.errcode != CL_SUCCESS){
				environ.setupCorrect = false;
			}
			else{
				environ.setupCorrect = true;
			}

			//Allocate memory for kernels and kernel flags
			environ.dpKernels = (cl_kernel*)malloc(NUM_KERNELS * sizeof(cl_kernel));
			environ.spKernels = (cl_kernel*)malloc(NUM_KERNELS * sizeof(cl_kernel));
			environ.dpKernelsAvailable = (bool*)malloc(NUM_KERNELS * sizeof(bool));
			environ.spKernelsAvailable = (bool*)malloc(NUM_KERNELS * sizeof(bool));

			//Set kernel flags to default false
			for(int i = 0; i < NUM_KERNELS; i++){
				environ.dpKernelsAvailable[i] = false;
				environ.spKernelsAvailable[i] = false;
			}

			//Compile all of the kernels
			if(errcode == CL_SUCCESS){
				for(int i = 0; i < NUM_KERNELS; i++){
					environ.dpKernels[i] = oclCreateKernel(
						dpKernelSources[i],
						dpKernelNames[i],
						environ.context);
					if(errcode == CL_SUCCESS){environ.dpKernelsAvailable[i] = true;}
				}

				for(int i = 0; i < NUM_KERNELS; i++){
					environ.spKernels[i] = oclCreateKernel(
						spKernelSources[i],
						spKernelNames[i],
						environ.context);
					if(errcode == CL_SUCCESS){environ.spKernelsAvailable[i] = true;}
				}
			}

			//Allocate and initialize the device mutex
			environ.deviceLock = (pthread_mutex_t*)malloc(sizeof(pthread_mutex_t));
			pthread_mutex_init(environ.deviceLock, NULL);

			return environ;
		}

		/**
		 * @internal
		 * Builds the vector of oclEnvirons
		 */
		static std::vector<oclEnviron>* oclGetEnvirons(
			std::vector<oclEnviron>* environs,
			cl_platform_id& platform,
			cl_device_id* devices,
			cl_uint numDevices){

			//Compute a score for all devices
			std::vector<long> rankings = oclComputeDeviceScores(devices, numDevices);

			//Find the max score and calculate a lower bound for device selection
			long maxScore = *(std::max_element(rankings.begin(), rankings.end()));
			long lowerBound = (long)((double)maxScore * 0.90);

			//Build OpenCL compute environments only for devices within 10% of
			//the top device. This should keep run times for identical computations
			//relatively similar across different intstances of OpenCLMatrixDomain
			//and eliminate the primary graphics adapter if there is one.
			for(int i = 0; i < (int)numDevices; i++){
				if(rankings.at(i) < lowerBound){
					continue;
				}

				oclEnviron temp;

				temp = oclBuildEnviron(temp, platform, devices[i]);

				environs->push_back(temp);
			}

			return environs;
		}

		/**
		 * @internal
		 * Initializes the OpenCL compute environment
		 */
		static void oclEnvironInit(){
			//Declare OpenCL specific variables
			cl_platform_id platform;
			cl_uint numDevices = 0;
			cl_device_id* devices;

			//Get the platform to be used
			errcode = oclGetPlatformID(platform);

			//Proceed only if successful with previous phase
			if(errcode == CL_SUCCESS){
				//Get number of GPU's in the platform
				errcode = clGetDeviceIDs(
					platform,
					CL_DEVICE_TYPE_GPU,
					0,
					NULL,
					&numDevices);
			}

			//Allocate the oclEnviron and instance count arrays
			environs = new std::vector<oclEnviron>();
			instances = new std::vector<int>();

			//Proceed only if successful with previous phase
			if(errcode == CL_SUCCESS && numDevices != 0){
				//Allocate memory for device array and read the devices into it
				devices = (cl_device_id*)malloc(numDevices * sizeof(cl_device_id));
				errcode = clGetDeviceIDs(
					platform,
					CL_DEVICE_TYPE_GPU,
					numDevices,
					devices,
					NULL);

				//Build environs from devices
				environs = oclGetEnvirons(environs, platform, devices, numDevices);

				//Set intial instance count to zero
				for(int i = 0; i < (int)environs->size(); i++){
					instances->push_back(0);
				}

				free(devices);
			}
			//If there are no devices or if an error occured make a default
			//oclEnviron that disables the OpenCLMatrixDomain and defaults it
			//to BlasMatrixDomain methods.
			else{
				oclEnviron temp;

				temp.context = NULL;
				temp.device = NULL;
				temp.commandQue = NULL;
				temp.errcode = CL_INVALID_PLATFORM;
				temp.memCapacity = 0UL;
				temp.maxBufferSize = 0UL;
				temp.GPUcontainer = false;
				temp.CPUcontainer = false;
				temp.setupCorrect = false;
				temp.doubleSupported = false;

				temp.dpKernels = (cl_kernel*)malloc(NUM_KERNELS * sizeof(cl_kernel));
				temp.spKernels = (cl_kernel*)malloc(NUM_KERNELS * sizeof(cl_kernel));
				temp.dpKernelsAvailable = (bool*)malloc(NUM_KERNELS * sizeof(bool));
				temp.spKernelsAvailable = (bool*)malloc(NUM_KERNELS * sizeof(bool));

				for(int i = 0; i < NUM_KERNELS; i++){
					temp.dpKernels[i] = NULL;
					temp.spKernels[i] = NULL;
					temp.dpKernelsAvailable[i] = false;
					temp.spKernelsAvailable[i] = false;
				}

				temp.deviceLock = NULL;

				environs->push_back(temp);
				instances->push_back(0);
			}

		}

		/**
		 * @internal
		 * Handles the releasing of resources at program exit
		 */
		static void oclResourceCleanUp(){
			//Release all reources held in each oclEnviron
			for(int i = 0; i < (int)environs->size(); i++){
				oclEnviron current = environs->at(i);

				//Release the kernels
				for(int j = 0; j < NUM_KERNELS; j++){
					errcode = clReleaseKernel(current.dpKernels[j]);
					errcode = clReleaseKernel(current.spKernels[j]);
				}

				//Deallocate the array memory
				free(current.dpKernels);
				free(current.spKernels);
				free(current.dpKernelsAvailable);
				free(current.spKernelsAvailable);

				//Destroy mutex and deallocate memory
				pthread_mutex_destroy(current.deviceLock);
				free(current.deviceLock);

				//Release the command queue and context
				errcode = clReleaseCommandQueue(current.commandQue);
				errcode = clReleaseContext(current.context);
			}


			//Delete vectors
			delete instances;
			delete environs;

			//Destroy factory mutex
			pthread_mutex_destroy(&factoryLock);

			//Set state to unitialized incase exit is aborted
			initialized = false;
		}

	public:
		template <class Field>
		static void oclMatrixDomainInstance(OpenCLMatrixDomain<Field>* target){

			//Check if the OpenCL environment has been initialized
			//It will only need to be initialized when the first OpenCLMatrixDomain
			//instance is created. Double check locking ensures Singleton
			if(!initialized){
				pthread_mutex_lock(&factoryLock);

				if(!initialized){
					oclEnvironInit();
					initialized = true;
					atexit(oclResourceCleanUp);
				}

				pthread_mutex_unlock(&factoryLock);
			}

			pthread_mutex_lock(&factoryLock);

			//Selected least used oclEnviron
			int leastUsedIndex = 0;
			for(int i = 0; i < (int)instances->size(); i++){
				if(instances->at(i) < instances->at(leastUsedIndex)){
					leastUsedIndex = i;
				}
			}

			//Increment use count
			(instances->at(leastUsedIndex))++;

			//Copy all of the data required for the OpenCLMatrixDomain instance to
			//function
			target->context = environs->at(leastUsedIndex).context;
			target->device = environs->at(leastUsedIndex).device;
			target->commandQue = environs->at(leastUsedIndex).commandQue;
			target->errcode = environs->at(leastUsedIndex).errcode;

			target->memCapacity = environs->at(leastUsedIndex).memCapacity;
			target->maxBufferSize = environs->at(leastUsedIndex).maxBufferSize;

			target->GPUcontainer = environs->at(leastUsedIndex).GPUcontainer;
			target->CPUcontainer = environs->at(leastUsedIndex).CPUcontainer;
			target->setupCorrect = environs->at(leastUsedIndex).setupCorrect;
			target->doubleSupported = environs->at(leastUsedIndex).doubleSupported;

			for(int i = 0; i < 20; i++){
				target->dpKernels[i] = environs->at(leastUsedIndex).dpKernels[i];
				target->spKernels[i] = environs->at(leastUsedIndex).spKernels[i];

				target->dpKernelsAvailable[i] =
					environs->at(leastUsedIndex).dpKernelsAvailable[i];
				target->spKernelsAvailable[i] =
					environs->at(leastUsedIndex).spKernelsAvailable[i];
			}

			//Assign an ID number the OpenCLMatrixDomain instance to be used for
			//locking and releasing the OpenCL resources
			target->IDnum = leastUsedIndex;

			//Point OpenCLMatrixDomain to the mutex
			target->deviceLock = environs->at(leastUsedIndex).deviceLock;

			pthread_mutex_unlock(&factoryLock);
		}

		static void oclMatrixDomainDeallocate(unsigned int IDnum){
			pthread_mutex_lock(&factoryLock);

			(instances->at(IDnum))--;

			pthread_mutex_unlock(&factoryLock);
		}

		static int oclGetNumberOfDevices(){

			//Check if the OpenCL environment has been initialized
			//It will only need to be initialized when the first OpenCLMatrixDomain
			//instance is created. Double check locking ensures Singleton
			if(!initialized){
				pthread_mutex_lock(&factoryLock);

				if(!initialized){
					oclEnvironInit();
					initialized = true;
					atexit(oclResourceCleanUp);
				}

				pthread_mutex_unlock(&factoryLock);
			}

			return (int) environs->size();

		}
	};

	//Static data fields used by OpenCLMatrixDomainFactory
	std::vector<OpenCLMatrixDomainFactory::oclEnviron>* OpenCLMatrixDomainFactory::environs = NULL;
	std::vector<int>* OpenCLMatrixDomainFactory::instances = NULL;
	bool OpenCLMatrixDomainFactory::initialized = false;
	cl_int OpenCLMatrixDomainFactory::errcode = CL_SUCCESS;
	pthread_mutex_t OpenCLMatrixDomainFactory::factoryLock = PTHREAD_MUTEX_INITIALIZER;
}

#endif // __LINBOX_opencl_domain_fatory

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
