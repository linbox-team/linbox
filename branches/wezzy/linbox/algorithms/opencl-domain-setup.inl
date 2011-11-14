/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/algorithms/opencl-domain-setup.inl
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

#ifndef __LINBOX_opencl_matrix_domain_setup_INL
#define __LINBOX_opencl_matrix_domain_setup_INL

#include "CL/cl.hpp"
//#include "helper_functions.hpp" -- For debugging only

#include <iostream>
#include <new>
#include <cstring>

namespace LinBox{

	/**
	 * @internal
	 * Picks the platform used for the container
	 */
	template<class Field>
	cl_int OpenCLMatrixDomain<Field>::oclGetPlatformID(cl_platform_id* selectedPlatform){

		char chBuffer[1024];
		cl_uint numPlatforms;
		cl_platform_id* platforms;
		*selectedPlatform = NULL;

		//OpenCL Platform count
		errcode = clGetPlatformIDs(0, NULL, &numPlatforms);
		if(errcode != CL_SUCCESS){

			return -1000;
		}
		if(numPlatforms == 0){

			return -2000;
		}

		platforms = (cl_platform_id*)operator new(numPlatforms * sizeof(cl_platform_id));
		errcode = clGetPlatformIDs(numPlatforms, platforms, NULL);

		for(unsigned int i = 0; i < numPlatforms; i++){
			errcode= clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, 1024, &chBuffer, NULL);

			if(errcode == CL_SUCCESS && !(strcmp(chBuffer, "NVIDIA CUDA"))){
				*selectedPlatform = platforms[i];

				GPUcontainer = false;
				CPUcontainer = true;
				break;
			}
		}
		if(*selectedPlatform == NULL){
			*selectedPlatform = platforms[0];

			GPUcontainer = false;
			CPUcontainer = false;
		}
		delete platforms;
		return CL_SUCCESS;
	}

	/**
	 * @internal
	 * Picks the device used for the container
	 */
	template<class Field>
	cl_device_id OpenCLMatrixDomain<Field>::oclDeviceSelector(cl_int numDevices, cl_device_id* devices){

		if(numDevices == 1){
			cl_device_type type;
			errcode = clGetDeviceInfo(devices[0], CL_DEVICE_TYPE, sizeof(cl_device_type), &type, NULL);

			//Set container type
			GPUcontainer = (type == CL_DEVICE_TYPE_GPU);
			CPUcontainer = (type == CL_DEVICE_TYPE_CPU);

			return devices[0];
		}

		int rankings[numDevices];
		int selected = 0;

		//Query device info and compute weighted score
		for(int i = 0; i < numDevices; i++){
			cl_uint computeUnits;
			errcode = clGetDeviceInfo(devices[i], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint),
				&computeUnits, NULL);

			cl_uint clockFrequency;
			errcode = clGetDeviceInfo(devices[i], CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(cl_uint),
				&clockFrequency, NULL);

			cl_device_type type;
			errcode = clGetDeviceInfo(devices[i], CL_DEVICE_TYPE, sizeof(cl_device_type),
				&type, NULL);

			cl_ulong maxGlobalMemoryAllocSize;
			errcode = clGetDeviceInfo(devices[i], CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(cl_ulong),
				&maxGlobalMemoryAllocSize, NULL);

			cl_ulong globaMemory;
			errcode = clGetDeviceInfo(devices[i], CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong),
				&globaMemory, NULL);

			rankings[i] = ((computeUnits * clockFrequency) * (type == CL_DEVICE_TYPE_GPU ? 4 : 1) +
				(maxGlobalMemoryAllocSize / (1024 * 1024)) + (globaMemory / (1024 * 1024)));
		}

		for(int i = 0; i < numDevices; i++){
			if(rankings[selected] < rankings[i]){
				selected = i;
			}
		}


		cl_device_type type;
		errcode = clGetDeviceInfo(devices[selected], CL_DEVICE_TYPE, sizeof(cl_device_type),
			&type, NULL);

		//Set container type
		GPUcontainer = (type == CL_DEVICE_TYPE_GPU);
		CPUcontainer = (type == CL_DEVICE_TYPE_CPU);

		//Determine double precision support for device
		doubleSupported = false;

		size_t sizeReturn;

		errcode = clGetDeviceInfo(devices[selected], CL_DEVICE_EXTENSIONS, 0 , NULL, &sizeReturn);
		char* device_extensions = (char*) operator new(sizeReturn);
		errcode = clGetDeviceInfo(devices[selected], CL_DEVICE_EXTENSIONS, sizeReturn, device_extensions, NULL);

		char* exten[200];

		exten[0] = strtok(device_extensions, " ");

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

		delete device_extensions;

		return devices[selected];
	}

	/**
	 * @internal
	 * Loads the contents of the specified file into memory
	 * Returns a pointer to a char array and the length of the file
	 */
	template<class Field>
	char* OpenCLMatrixDomain<Field>::readFileContents(const char* fileName, int &length){

		//Open file
		std::ifstream input;
		input.open(fileName);

		//Get file size
		input.seekg(0, std::ios::end);
		length = input.tellg();
		input.seekg(0, std::ios::beg);

		//Allocate memory for file contents
		char* buffer = new char[(length + 1)];

		//Read the file
		input.read(buffer, length);

		//Close the file
		input.close();

		//Append null character to end of file contents
		buffer[length] = '\0';

		//Return the buffer
		return buffer;
	}

	/**
	 * @internal
	 * Creates a kernel given a file name and kernel name
	 * Returns the kernel
	 */
	template<class Field>
	cl_kernel OpenCLMatrixDomain<Field>::oclCreateKernel(const char* fileName, const char* kernelName){

		//Load the file and get length
		int fileLength;
		char* fileContents = readFileContents(fileName, fileLength);

		//Create program from file
		size_t kernelLength = fileLength;
		cl_program program = clCreateProgramWithSource(context,
				1, (const char**)&fileContents,
				&kernelLength, &errcode);

		//Build the program into executable
		errcode = clBuildProgram(program, 0,
				NULL, NULL, NULL, NULL);

		//Create kernel from executable
		cl_kernel tempKernel = clCreateKernel(program,
						kernelName, &errcode);

		//Releasing program
		errcode = clReleaseProgram(program);

		//Return kernel
		return tempKernel;
	}

	/**
	 * @internal
	 * Initializes the OpenCL compute environment
	 */
	template<class Field>
	void OpenCLMatrixDomain<Field>::oclDomainInit(){

		//Declare OpenCL specific variables
		cl_platform_id platform;
		cl_uint numDevices;
		cl_device_id* devices;

		//Get platform id
		errcode = oclGetPlatformID(&platform);

		//Proceed only if successful with previous phase
		if(errcode == CL_SUCCESS){
			//Get number of devices
			if(GPUcontainer){
				errcode = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 0, 0, &numDevices);
			}
			else{
				errcode = clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, 0, 0, &numDevices);
			}
		}

		//Allocate memory for device array
		devices = (cl_device_id*)operator new(numDevices * sizeof(cl_device_id));

		//Proceed only if successful with previous phase
		if(errcode == CL_SUCCESS){
			//Get device information
			if(GPUcontainer){
				errcode = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, numDevices, devices, NULL);
			}
			else{
				errcode = clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, numDevices, devices, NULL);
			}
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

		if(errcode == CL_SUCCESS){
			//Provides the code to load and create the kernels on the compiling platform
			//Does not allow for compile once and reuse the host code binary unless
			//the file paths are identical
			#include "linbox/algorithms/opencl-kernels/opencl-domain-file-paths.inl"
		}

		if(errcode != CL_SUCCESS){
			setupCorrect = false;
			//std::cout << "False\n";
		}
		else{
			setupCorrect = true;
			//std::cout << "True\n";
		}
	}

	/**
	 * @internal
	 * Releases OpenCL cumpute resources
	 */
	template<class Field>
	void OpenCLMatrixDomain<Field>::oclDomainTearDown(){
		//Release all of the kernels
		errcode = clReleaseKernel(dpKernels[0]);
		errcode = clReleaseKernel(dpKernels[1]);
		errcode = clReleaseKernel(dpKernels[2]);
		errcode = clReleaseKernel(dpKernels[3]);
		errcode = clReleaseKernel(dpKernels[4]);
		errcode = clReleaseKernel(dpKernels[5]);
		errcode = clReleaseKernel(dpKernels[6]);
		errcode = clReleaseKernel(dpKernels[7]);
		errcode = clReleaseKernel(dpKernels[8]);
		errcode = clReleaseKernel(dpKernels[9]);

		errcode = clReleaseKernel(spKernels[0]);
		errcode = clReleaseKernel(spKernels[1]);
		errcode = clReleaseKernel(spKernels[2]);
		errcode = clReleaseKernel(spKernels[3]);
		errcode = clReleaseKernel(spKernels[4]);
		errcode = clReleaseKernel(spKernels[5]);
		errcode = clReleaseKernel(spKernels[6]);
		errcode = clReleaseKernel(spKernels[7]);
		errcode = clReleaseKernel(spKernels[8]);
		errcode = clReleaseKernel(spKernels[9]);

		//Release the command queue
		errcode = clReleaseCommandQueue(commandQue);

		//Release the compute context
		errcode = clReleaseContext(context);
	}

}; //end of namespace LinBox

#endif // __LINBOX_opencl_matrix_domain_setup_INL