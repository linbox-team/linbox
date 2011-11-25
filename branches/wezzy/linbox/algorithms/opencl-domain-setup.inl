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
#include "linbox/algorithms/opencl-domain-factory.h"

#include <iostream>
#include <new>
#include <cstring>

namespace LinBox{

	/**
	 * @internal
	 * Initializes the OpenCL compute environment
	 */
	template<class Field>
	void OpenCLMatrixDomain<Field>::oclDomainInit(){

		OpenCLMatrixDomainFactory::oclDomainCreate(context, device,
			commandQue, errcode, memCapacity, maxBufferSize, GPUcontainer,
			CPUcontainer, setupCorrect, doubleSupported, dpKernels, spKernels);

		if(!setupCorrect){
			setupCorrect = false;
			//::cout << "False\n";
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

		OpenCLMatrixDomainFactory::oclDomainRelease();
	}

}; //end of namespace LinBox

#endif // __LINBOX_opencl_matrix_domain_setup_INL