/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/algorithms/opencl-domain-factory.h
 * Copyright (C) 2011 Matthew Wezowicz
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
 * Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

#ifndef __LINBOX_opencl_matrix_domain_factory_H
#define __LINBOX_opencl_matrix_domain_factory_H

#include <new>
#include <cstring>
#include "linbox/algorithms/opencl-domain.h"
#include "linbox/algorithms/opencl-kernels/opencl-domain-kernels.inl"

#include "CL/cl.hpp"

namespace LinBox{

	class OpenCLMatrixDomainFactory {

	protected:

		//OpenCL specific variables
		static cl_context context;
		static cl_device_id device;
		static cl_command_queue commandQue;
		static cl_int errcode;

		//Storage for memory levels
		static unsigned long memCapacity;
		static unsigned long maxBufferSize;

		//Type and Status flags
		static bool GPUcontainer;
		static bool CPUcontainer;
		static bool setupCorrect;
		static bool doubleSupported;
		static bool initialized;

		//Storage for kernels
		static cl_kernel* dpKernels;
		static bool* dpKernelsAvailable;
		static cl_kernel* spKernels;
		static bool* spKernelsAvailable;

		//Count of instances
		static int countOpenCLMatrixDomain;

		/**
		 * @internal
		 * Picks the platform used for the container
		 */
		static cl_int oclGetPlatformID(cl_platform_id& selectedPlatform){
			//Allocate temporary char array for platform name
			char chBuffer[256];
			cl_uint numPlatforms;
			cl_platform_id* platforms;
			selectedPlatform = NULL;

			//OpenCL Platform count return custom error codes if there are no platforms
			//or could not get number of platforms
			errcode = clGetPlatformIDs(0, NULL, &numPlatforms);
			if(errcode != CL_SUCCESS){
				return -1000;
			}
			if(numPlatforms == 0){
				return -2000;
			}

			//Allocate space to store cl_platform_id's
			platforms = (cl_platform_id*)operator new(numPlatforms * sizeof(cl_platform_id));
			errcode = clGetPlatformIDs(numPlatforms, platforms, NULL);

			//Search through the platforms looking of a prefered one specified by a string
			for(unsigned int i = 0; i < numPlatforms; i++){
				errcode= clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, 256, &chBuffer, NULL);

				if(errcode == CL_SUCCESS && !(strcmp(chBuffer, "NVIDIA CUDA"))){
					selectedPlatform = platforms[i];
					break;
				}
			}

			//If prefered platform could not be found use first platforms
			if(selectedPlatform == NULL){
				selectedPlatform = platforms[0];
			}

			//Clean up memory
			delete platforms;

			return CL_SUCCESS;
		}

		/**
		 * @internal
		 * Picks the device used for the container
		 */
		static cl_device_id oclDeviceSelector(cl_int numDevices, cl_device_id* devices){
			//If there is only one device on the selected platform
			if(numDevices == 1){
				cl_device_type type;
				errcode = clGetDeviceInfo(devices[0], CL_DEVICE_TYPE, sizeof(cl_device_type), &type, NULL);

				//Set container type
				GPUcontainer = (type == CL_DEVICE_TYPE_GPU);
				CPUcontainer = (type == CL_DEVICE_TYPE_CPU);

				return devices[0];
			}

			//Allocate space for the calculation of device scores
			int rankings[numDevices];
			int selected = 0;

			//Query device info and compute weighted score
			for(int i = 0; i < numDevices; i++){
				cl_device_type type;
				errcode = clGetDeviceInfo(devices[selected], CL_DEVICE_TYPE, sizeof(cl_device_type),
				&type, NULL);

				//Temp set container type
				GPUcontainer = (type == CL_DEVICE_TYPE_GPU);
				CPUcontainer = (type == CL_DEVICE_TYPE_CPU);

				//Get the number of processors, clock rate, device type, maximum buffere size, and total memory
				cl_uint computeUnits;
				errcode = clGetDeviceInfo(devices[i], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint),
					&computeUnits, NULL);

				cl_uint clockFrequency;
				errcode = clGetDeviceInfo(devices[i], CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(cl_uint),
					&clockFrequency, NULL);

				errcode = clGetDeviceInfo(devices[i], CL_DEVICE_TYPE, sizeof(cl_device_type),
					&type, NULL);

				cl_ulong maxGlobalMemoryAllocSize;
				errcode = clGetDeviceInfo(devices[i], CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(cl_ulong),
					&maxGlobalMemoryAllocSize, NULL);

				cl_ulong globalMemory;
				errcode = clGetDeviceInfo(devices[i], CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong),
					&globalMemory, NULL);

				//Calculate score for device -- Subject to modification
				rankings[i] = ((computeUnits * clockFrequency) * (type == CL_DEVICE_TYPE_GPU ? 4 : 1) +
					(maxGlobalMemoryAllocSize / (1024 * 1024)) + (globalMemory / (1024 * 1024)));
			}

			//Select the device with the highest score
			for(int i = 0; i < numDevices; i++){
				if(rankings[selected] < rankings[i]){
					selected = i;
				}
			}

			//Get the device type for the selected device
			cl_device_type type;
			errcode = clGetDeviceInfo(devices[selected], CL_DEVICE_TYPE, sizeof(cl_device_type),
				&type, NULL);

			//Set container type
			GPUcontainer = (type == CL_DEVICE_TYPE_GPU);
			CPUcontainer = (type == CL_DEVICE_TYPE_CPU);

			//Determine double precision support for device by getting a listing of all
			//OpenCL device extensions supported by the device and searching through them for
			//the string "cl_khr_fp64"
			doubleSupported = false;

			size_t sizeReturn;

			errcode = clGetDeviceInfo(devices[selected], CL_DEVICE_EXTENSIONS, 0, NULL, &sizeReturn);
			char* deviceExtensions = (char*) operator new(sizeReturn);
			errcode = clGetDeviceInfo(devices[selected], CL_DEVICE_EXTENSIONS, sizeReturn, deviceExtensions, NULL);

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
						doubleSupported = true;
					}
				}
			}

			delete deviceExtensions;

			return devices[selected];
		}

		/**
		 * @internal
		 * Loads the contents of the specified file into memory
		 * Returns a pointer to a char array and the length of the file
		 */
		// static char* readFileContents(const char* fileName, int& length){

			// //Open file
			// std::ifstream input;
			// input.open(fileName);

			// //Get file size
			// input.seekg(0, std::ios::end);
			// length = input.tellg();
			// input.seekg(0, std::ios::beg);

			// //Allocate memory for file contents
			// char* buffer = new char[(length + 1)];

			// //Read the file
			// input.read(buffer, length);

			// //Close the file
			// input.close();

			// //Append null character to end of file contents
			// buffer[length] = '\0';

			// //Return the buffer
			// return buffer;
		// }

		/**
		 * @internal
		 * Creates a kernel given a file name and kernel name
		 * Returns the kernel
		 */
		static cl_kernel oclCreateKernel(const char* kernel, const char* kernelName){

			//Load the file and get length
			//int fileLength;
			//char* fileContents = readFileContents(fileName, fileLength);

			//find length of the kernel
			size_t kernelLength = strlen(kernel);
			
			//Create program from file
			//size_t kernelLength = fileLength;
			cl_program program = clCreateProgramWithSource(context, 1,
				(const char**)&kernel, &kernelLength, &errcode);

			//Build the program into executable
			if(errcode == CL_SUCCESS){
				errcode = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
			}
			/*
			printf(kernelName);
			printf("\n");
			size_t ret;
			errcode = clGetProgramBuildInfo(program,device,CL_PROGRAM_BUILD_LOG,0,NULL,&ret);
			char* log = (char*)operator new(ret);
			errcode = clGetProgramBuildInfo(program,device,CL_PROGRAM_BUILD_LOG,ret,log,NULL);
			printf(log);
			delete log;
			*/

			//Create kernel from executable
			cl_kernel tempKernel = NULL;
			if(errcode == CL_SUCCESS){
				tempKernel = clCreateKernel(program,
					kernelName, &errcode);
			}

			//Releasing program
			clReleaseProgram(program);

			//Return kernel
			return tempKernel;
		}


		/**
		 * @internal
		 * Initializes the OpenCL compute environment
		 */
		static void oclEnvironInit(){

			//Declare OpenCL specific variables
			cl_platform_id platform;
			cl_uint numDevices;
			cl_device_id* devices;

			//Get platform id
			errcode = oclGetPlatformID(platform);

			//Proceed only if successful with previous phase
			if(errcode == CL_SUCCESS){
				//Get number of devices
				errcode = clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, 0, 0, &numDevices);
			}

			//Allocate memory for device array
			devices = (cl_device_id*)operator new(numDevices * sizeof(cl_device_id));

			//Proceed only if successful with previous phase
			if(errcode == CL_SUCCESS){
				//Get device information
				errcode = clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, numDevices, devices, NULL);
			}

			//Proceed only if successful with previous phase
			if(errcode == CL_SUCCESS){
				//Select Device
				device = oclDeviceSelector(numDevices, devices);
			}

			//Proceed only if successful with previous phase
			if(errcode == CL_SUCCESS){
				//Create compute context
				context = clCreateContext(0, 1, &device, NULL, NULL, &errcode);
			}

			//Proceed only if successful with previous phase
			if(errcode == CL_SUCCESS){
				//Create command queue for device
				commandQue = clCreateCommandQueue(context, device, 0, &errcode);
			}

			delete devices;

			cl_ulong memSize;
			cl_ulong maxGlobalMemoryAllocSize;

			//Proceed only if successful with previous phase
			if(errcode == CL_SUCCESS){
				//Get amount of memory that the device has
				errcode = clGetDeviceInfo(device, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(memSize),
					&memSize, NULL);
				errcode = clGetDeviceInfo(device, CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(cl_ulong),
					&maxGlobalMemoryAllocSize, NULL);
			}

			memCapacity = (unsigned long)memSize;
			maxBufferSize = (unsigned long)maxGlobalMemoryAllocSize;

			//Initialize al kernel compilation flags to false
			for(int i = 0; i < 22; i++){
				dpKernelsAvailable[i] = false;
				spKernelsAvailable[i] = false;
			}

			//Compile all of the kernels
			if(errcode == CL_SUCCESS){
				dpKernels[0] = oclCreateKernel(addKernelModularDP, "addKernelModularDP");
				if(errcode == CL_SUCCESS){dpKernelsAvailable[0] = true;}

				dpKernels[1] = oclCreateKernel(subKernelModularDP, "subKernelModularDP");
				if(errcode == CL_SUCCESS){dpKernelsAvailable[1] = true;}

				dpKernels[2] = oclCreateKernel(matrixMulKernelModular1DP, "matrixMulKernelModular1DP");
				if(errcode == CL_SUCCESS){dpKernelsAvailable[2] = true;}

				dpKernels[3] = oclCreateKernel(matrixMulKernelModular8DP, "matrixMulKernelModular8DP");
				if(errcode == CL_SUCCESS){dpKernelsAvailable[3] = true;}

				dpKernels[4] = oclCreateKernel(matrixMulKernelModular32DP, "matrixMulKernelModular32DP");
				if(errcode == CL_SUCCESS){dpKernelsAvailable[4] = true;}

				dpKernels[5] = oclCreateKernel(matrixMulKernelModular1024DP, "matrixMulKernelModular1024DP");
				if(errcode == CL_SUCCESS){dpKernelsAvailable[5] = true;}

				dpKernels[6] = oclCreateKernel(matrixMuladdKernelModular1DP, "matrixMuladdKernelModular1DP");
				if(errcode == CL_SUCCESS){dpKernelsAvailable[6] = true;}

				dpKernels[7] = oclCreateKernel(matrixMuladdKernelModular8DP, "matrixMuladdKernelModular8DP");
				if(errcode == CL_SUCCESS){dpKernelsAvailable[7] = true;}

				dpKernels[8] = oclCreateKernel(matrixMuladdKernelModular32DP, "matrixMuladdKernelModular32DP");
				if(errcode == CL_SUCCESS){dpKernelsAvailable[8] = true;}

				dpKernels[9] = oclCreateKernel(matrixMuladdKernelModular1024DP, "matrixMuladdKernelModular1024DP");
				if(errcode == CL_SUCCESS){dpKernelsAvailable[9] = true;}

				dpKernels[10] = oclCreateKernel(matrixAxpyKernelModular1DP, "matrixAxpyKernelModular1DP");
				if(errcode == CL_SUCCESS){dpKernelsAvailable[10] = true;}

				dpKernels[11] = oclCreateKernel(matrixAxpyKernelModular8DP, "matrixAxpyKernelModular8DP");
				if(errcode == CL_SUCCESS){dpKernelsAvailable[11] = true;}

				dpKernels[12] = oclCreateKernel(matrixAxpyKernelModular32DP, "matrixAxpyKernelModular32DP");
				if(errcode == CL_SUCCESS){dpKernelsAvailable[12] = true;}

				dpKernels[13] = oclCreateKernel(matrixAxpyKernelModular1024DP, "matrixAxpyKernelModular1024DP");
				if(errcode == CL_SUCCESS){dpKernelsAvailable[13] = true;}

				dpKernels[14] = oclCreateKernel(matrixMaxpyKernelModular1DP, "matrixMaxpyKernelModular1DP");
				if(errcode == CL_SUCCESS){dpKernelsAvailable[14] = true;}

				dpKernels[15] = oclCreateKernel(matrixMaxpyKernelModular8DP, "matrixMaxpyKernelModular8DP");
				if(errcode == CL_SUCCESS){dpKernelsAvailable[15] = true;}

				dpKernels[16] = oclCreateKernel(matrixMaxpyKernelModular32DP, "matrixMaxpyKernelModular32DP");
				if(errcode == CL_SUCCESS){dpKernelsAvailable[16] = true;}

				dpKernels[17] = oclCreateKernel(matrixMaxpyKernelModular1024DP, "matrixMaxpyKernelModular1024DP");
				if(errcode == CL_SUCCESS){dpKernelsAvailable[17] = true;}

				dpKernels[18] = oclCreateKernel(matrixAxmyKernelModular1DP, "matrixAxmyKernelModular1DP");
				if(errcode == CL_SUCCESS){dpKernelsAvailable[18] = true;}

				dpKernels[19] = oclCreateKernel(matrixAxmyKernelModular8DP, "matrixAxmyKernelModular8DP");
				if(errcode == CL_SUCCESS){dpKernelsAvailable[19] = true;}

				dpKernels[20] = oclCreateKernel(matrixAxmyKernelModular32DP, "matrixAxmyKernelModular32DP");
				if(errcode == CL_SUCCESS){dpKernelsAvailable[20] = true;}

				dpKernels[21] = oclCreateKernel(matrixAxmyKernelModular1024DP, "matrixAxmyKernelModular1024DP");
				if(errcode == CL_SUCCESS){dpKernelsAvailable[21] = true;}


				spKernels[0] = oclCreateKernel(addKernelModularSP, "addKernelModularSP");
				if(errcode == CL_SUCCESS){spKernelsAvailable[0] = true;}

				spKernels[1] = oclCreateKernel(subKernelModularSP, "subKernelModularSP");
				if(errcode == CL_SUCCESS){spKernelsAvailable[1] = true;}

				spKernels[2] = oclCreateKernel(matrixMulKernelModular1SP, "matrixMulKernelModular1SP");
				if(errcode == CL_SUCCESS){spKernelsAvailable[2] = true;}

				spKernels[3] = oclCreateKernel(matrixMulKernelModular16SP, "matrixMulKernelModular16SP");
				if(errcode == CL_SUCCESS){spKernelsAvailable[3] = true;}

				spKernels[4] = oclCreateKernel(matrixMulKernelModular32SP, "matrixMulKernelModular32SP");
				if(errcode == CL_SUCCESS){spKernelsAvailable[4] = true;}

				spKernels[5] = oclCreateKernel(matrixMulKernelModular1024SP, "matrixMulKernelModular1024SP");
				if(errcode == CL_SUCCESS){spKernelsAvailable[5] = true;}

				spKernels[6] = oclCreateKernel(matrixMuladdKernelModular1SP, "matrixMuladdKernelModular1SP");
				if(errcode == CL_SUCCESS){spKernelsAvailable[6] = true;}

				spKernels[7] = oclCreateKernel(matrixMuladdKernelModular16SP, "matrixMuladdKernelModular16SP");
				if(errcode == CL_SUCCESS){spKernelsAvailable[7] = true;}

				spKernels[8] = oclCreateKernel(matrixMuladdKernelModular32SP, "matrixMuladdKernelModular32SP");
				if(errcode == CL_SUCCESS){spKernelsAvailable[8] = true;}

				spKernels[9] = oclCreateKernel(matrixMuladdKernelModular1024SP, "matrixMuladdKernelModular1024SP");
				if(errcode == CL_SUCCESS){spKernelsAvailable[9] = true;}

				spKernels[10] = oclCreateKernel(matrixAxpyKernelModular1SP, "matrixAxpyKernelModular1SP");
				if(errcode == CL_SUCCESS){spKernelsAvailable[10] = true;}

				spKernels[11] = oclCreateKernel(matrixAxpyKernelModular16SP, "matrixAxpyKernelModular16SP");
				if(errcode == CL_SUCCESS){spKernelsAvailable[11] = true;}

				spKernels[12] = oclCreateKernel(matrixAxpyKernelModular32SP, "matrixAxpyKernelModular32SP");
				if(errcode == CL_SUCCESS){spKernelsAvailable[12] = true;}

				spKernels[13] = oclCreateKernel(matrixAxpyKernelModular1024SP, "matrixAxpyKernelModular1024SP");
				if(errcode == CL_SUCCESS){spKernelsAvailable[13] = true;}

				spKernels[14] = oclCreateKernel(matrixMaxpyKernelModular1SP, "matrixMaxpyKernelModular1SP");
				if(errcode == CL_SUCCESS){spKernelsAvailable[14] = true;}

				spKernels[15] = oclCreateKernel(matrixMaxpyKernelModular16SP, "matrixMaxpyKernelModular16SP");
				if(errcode == CL_SUCCESS){spKernelsAvailable[15] = true;}

				spKernels[16] = oclCreateKernel(matrixMaxpyKernelModular32SP, "matrixMaxpyKernelModular32SP");
				if(errcode == CL_SUCCESS){spKernelsAvailable[16] = true;}

				spKernels[17] = oclCreateKernel(matrixMaxpyKernelModular1024SP, "matrixMaxpyKernelModular1024SP");
				if(errcode == CL_SUCCESS){spKernelsAvailable[17] = true;}

				spKernels[18] = oclCreateKernel(matrixAxmyKernelModular1SP, "matrixAxmyKernelModular1SP");
				if(errcode == CL_SUCCESS){spKernelsAvailable[18] = true;}

				spKernels[19] = oclCreateKernel(matrixAxmyKernelModular16SP, "matrixAxmyKernelModular16SP");
				if(errcode == CL_SUCCESS){spKernelsAvailable[19] = true;}

				spKernels[20] = oclCreateKernel(matrixAxmyKernelModular32SP, "matrixAxmyKernelModular32SP");
				if(errcode == CL_SUCCESS){spKernelsAvailable[20] = true;}

				spKernels[21] = oclCreateKernel(matrixAxmyKernelModular1024SP, "matrixAxmyKernelModular1024SP");
				if(errcode == CL_SUCCESS){spKernelsAvailable[21] = true;}
			}

			//Check if everthing is setup correctly
			if(errcode != CL_SUCCESS){
				setupCorrect = false;
			}
			else{
				setupCorrect = true;
			}
		}

	public:
		template <class Field>
		static void oclMatrixDomainCreate(OpenCLMatrixDomain<Field>* target){

			//Check if the OpenCL environment has been initialized
			//It will only need to be initialized when the first OpenCLMatrixDomain instance
			//is created otherwise it will just reuse the environment
			if(!initialized){
				oclEnvironInit();
				initialized = true;
			}

			//Copy all of the data required for the OpenCLMatrixDomain instance to function
			target->context = context;
			target->device = device;
			target->commandQue = commandQue;
			target->errcode = errcode;

			target->memCapacity = memCapacity;
			target->maxBufferSize = maxBufferSize;

			target->GPUcontainer = GPUcontainer;
			target->CPUcontainer = CPUcontainer;
			target->setupCorrect = setupCorrect;
			target->doubleSupported = doubleSupported;

			for(int i = 0; i < 22; i++){
				target->dpKernels[i] = dpKernels[i];
				target->dpKernelsAvailable[i] = dpKernelsAvailable[i];
				target->spKernels[i] = spKernels[i];
				target->spKernelsAvailable[i] = spKernelsAvailable[i];
			}

			//Assign an ID number the OpenCLMatrixDomain instance to be used for locking and
			//releasing the OpenCL resources
			target->IDnum = countOpenCLMatrixDomain;

			//Increase count of OpenCLMatrixDomain instances
			countOpenCLMatrixDomain++;
		}

		static void oclMatrixDomainDestroy(unsigned int IDnum){
			//Decrease count of OpenCLMatrixDomain instances
			countOpenCLMatrixDomain--;
		}
	};

	//Initialization of static members to default values
	cl_context OpenCLMatrixDomainFactory::context = NULL;
	cl_device_id OpenCLMatrixDomainFactory::device = NULL;
	cl_command_queue OpenCLMatrixDomainFactory::commandQue = NULL;
	cl_int OpenCLMatrixDomainFactory::errcode = CL_SUCCESS;

	unsigned long OpenCLMatrixDomainFactory::memCapacity = 0L;
	unsigned long OpenCLMatrixDomainFactory::maxBufferSize = 0L;

	bool OpenCLMatrixDomainFactory::GPUcontainer = false;
	bool OpenCLMatrixDomainFactory::CPUcontainer = false;
	bool OpenCLMatrixDomainFactory::setupCorrect = false;
	bool OpenCLMatrixDomainFactory::doubleSupported = false;
	bool OpenCLMatrixDomainFactory::initialized = false;

	cl_kernel* OpenCLMatrixDomainFactory::dpKernels = (cl_kernel*)operator new(22 * sizeof(cl_kernel));
	bool* OpenCLMatrixDomainFactory::dpKernelsAvailable = (bool*)operator new(22 * sizeof(bool));
	cl_kernel* OpenCLMatrixDomainFactory::spKernels = (cl_kernel*)operator new(22 * sizeof(cl_kernel));
	bool* OpenCLMatrixDomainFactory::spKernelsAvailable = (bool*)operator new(22 * sizeof(bool));

	int OpenCLMatrixDomainFactory::countOpenCLMatrixDomain = 0;

}

#endif