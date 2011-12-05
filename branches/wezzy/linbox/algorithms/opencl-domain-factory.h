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

#include "CL/cl.hpp"
//#include "helper_functions.hpp" -- For debugging only

#include <iostream>
#include <new>
#include <cstring>

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
		static cl_kernel* spKernels;

		//Count of instances
		static int countOpenCLMatrixDomain;

		/**
		 * @internal
		 * Picks the platform used for the container
		 */
		static cl_int oclGetPlatformID(cl_platform_id& selectedPlatform){

			char chBuffer[256];
			cl_uint numPlatforms;
			cl_platform_id* platforms;
			selectedPlatform = NULL;

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
				errcode= clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, 256, &chBuffer, NULL);

				if(errcode == CL_SUCCESS && !(strcmp(chBuffer, "NVIDIA CUDA"))){
					selectedPlatform = platforms[i];
					break;
				}
			}
			if(selectedPlatform == NULL){
				selectedPlatform = platforms[0];
			}
			delete platforms;
			return CL_SUCCESS;
		}

		/**
		 * @internal
		 * Picks the device used for the container
		 */
		static cl_device_id oclDeviceSelector(cl_int numDevices, cl_device_id* devices){

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
				cl_device_type type;
				errcode = clGetDeviceInfo(devices[selected], CL_DEVICE_TYPE, sizeof(cl_device_type),
				&type, NULL);

				//Temp set container type
				GPUcontainer = (type == CL_DEVICE_TYPE_GPU);
				CPUcontainer = (type == CL_DEVICE_TYPE_CPU);

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

				rankings[i] = ((computeUnits * clockFrequency) * (type == CL_DEVICE_TYPE_GPU ? 4 : 1) +
					(maxGlobalMemoryAllocSize / (1024 * 1024)) + (globalMemory / (1024 * 1024)));
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
		static char* readFileContents(const char* fileName, int& length){

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
		static cl_kernel oclCreateKernel(const char* fileName, const char* kernelName){

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

			if(errcode == CL_SUCCESS){
				dpKernels[0] = oclCreateKernel("/home/mwezz/Downloads/linbox/linbox/algorithms/opencl-kernels/kernel_add_matrix_dp.cl", "vector_sum_kernel");
				dpKernels[1] = oclCreateKernel("/home/mwezz/Downloads/linbox/linbox/algorithms/opencl-kernels/kernel_sub_matrix_dp.cl", "vector_sum_kernel");
				dpKernels[2] = oclCreateKernel("/home/mwezz/Downloads/linbox/linbox/algorithms/opencl-kernels/kernel_modulus_dp.cl", "matrix_mul_kernel");
				dpKernels[3] = oclCreateKernel("/home/mwezz/Downloads/linbox/linbox/algorithms/opencl-kernels/kernel_partial_8_dp.cl", "matrix_mul_kernel");
				dpKernels[4] = oclCreateKernel("/home/mwezz/Downloads/linbox/linbox/algorithms/opencl-kernels/kernel_partial_32_dp.cl", "matrix_mul_kernel");
				dpKernels[5] = oclCreateKernel("/home/mwezz/Downloads/linbox/linbox/algorithms/opencl-kernels/kernel_partial_1024_dp.cl", "matrix_mul_kernel");
				dpKernels[6] = oclCreateKernel("/home/mwezz/Downloads/linbox/linbox/algorithms/opencl-kernels/kernel_muladd_modulus_dp.cl", "matrix_mul_kernel");
				dpKernels[7] = oclCreateKernel("/home/mwezz/Downloads/linbox/linbox/algorithms/opencl-kernels/kernel_muladd_partial_8_dp.cl", "matrix_mul_kernel");
				dpKernels[8] = oclCreateKernel("/home/mwezz/Downloads/linbox/linbox/algorithms/opencl-kernels/kernel_muladd_partial_32_dp.cl", "matrix_mul_kernel");
				dpKernels[9] = oclCreateKernel("/home/mwezz/Downloads/linbox/linbox/algorithms/opencl-kernels/kernel_muladd_partial_1024_dp.cl", "matrix_mul_kernel");

				spKernels[0] = oclCreateKernel("/home/mwezz/Downloads/linbox/linbox/algorithms/opencl-kernels/kernel_add_matrix_sp.cl", "vector_sum_kernel");
				spKernels[1] = oclCreateKernel("/home/mwezz/Downloads/linbox/linbox/algorithms/opencl-kernels/kernel_sub_matrix_sp.cl", "vector_sum_kernel");
				spKernels[2] = oclCreateKernel("/home/mwezz/Downloads/linbox/linbox/algorithms/opencl-kernels/kernel_modulus_sp.cl", "matrix_mul_kernel");
				spKernels[3] = oclCreateKernel("/home/mwezz/Downloads/linbox/linbox/algorithms/opencl-kernels/kernel_partial_16_sp.cl", "matrix_mul_kernel");
				spKernels[4] = oclCreateKernel("/home/mwezz/Downloads/linbox/linbox/algorithms/opencl-kernels/kernel_partial_32_sp.cl", "matrix_mul_kernel");
				spKernels[5] = oclCreateKernel("/home/mwezz/Downloads/linbox/linbox/algorithms/opencl-kernels/kernel_partial_1024_sp.cl", "matrix_mul_kernel");
				spKernels[6] = oclCreateKernel("/home/mwezz/Downloads/linbox/linbox/algorithms/opencl-kernels/kernel_muladd_modulus_sp.cl", "matrix_mul_kernel");
				spKernels[7] = oclCreateKernel("/home/mwezz/Downloads/linbox/linbox/algorithms/opencl-kernels/kernel_muladd_partial_16_sp.cl", "matrix_mul_kernel");
				spKernels[8] = oclCreateKernel("/home/mwezz/Downloads/linbox/linbox/algorithms/opencl-kernels/kernel_muladd_partial_32_sp.cl", "matrix_mul_kernel");
				spKernels[9] = oclCreateKernel("/home/mwezz/Downloads/linbox/linbox/algorithms/opencl-kernels/kernel_muladd_partial_1024_sp.cl", "matrix_mul_kernel");
			}

			if(errcode != CL_SUCCESS){
				setupCorrect = false;
				//std::cout << "Factory False\n";
			}
			else{
				setupCorrect = true;
				//std::cout << "Factory True\n";
			}
		}

	public:
		static void oclDomainCreate(cl_context& cont, cl_device_id& dev,
			cl_command_queue& queue, cl_int& err, unsigned long& memCap,
			unsigned long& maxBuf, bool& GPUflag, bool& CPUflag, bool& correct,
			bool& doubleSupport, cl_kernel* dpKern, cl_kernel* spKern){

			if(!initialized){
				oclEnvironInit();
				initialized = true;
			}

			cont = context;
			dev = device;
			queue = commandQue;
			err = errcode;

			memCap = memCapacity;
			maxBuf = maxBufferSize;

			GPUflag = GPUcontainer;
			CPUflag = CPUcontainer;
			correct = setupCorrect;
			doubleSupport = doubleSupported;

			for(int i = 0; i < 10; i++){
				dpKern[i] = dpKernels[i];
				spKern[i] = spKernels[i];
			}

			countOpenCLMatrixDomain++;
		}

		static void oclDomainRelease(){
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
	cl_kernel* OpenCLMatrixDomainFactory::dpKernels = (cl_kernel*)operator new(10 * sizeof(cl_kernel));
	cl_kernel* OpenCLMatrixDomainFactory::spKernels = (cl_kernel*)operator new(10 * sizeof(cl_kernel));;
	int OpenCLMatrixDomainFactory::countOpenCLMatrixDomain = 0;

}

#endif