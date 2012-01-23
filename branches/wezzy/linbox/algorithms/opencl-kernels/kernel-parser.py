"""
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/algorithms/opencl-domain-factory.h
 * Copyright (C) 2011-2012 Matthew Wezowicz
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
"""

import string
import time

def parse_kernel(file_name, output):
    # Open the kernel file
    fp = open(file_name)

    # Begin parsing of the kernel file
    parsed = []

    # Search for the kernel name in the file
    # Once found generate the const char* kernel name line
    for line in fp.readlines():
        line = string.split(line)
        if(len(line) > 0 and line[0] == '__kernel'):
            parsed.append("const char*")
            temp = string.split(line[2],'(')
            parsed.append(temp[0])
            parsed.append("= {\n")
            break
    
    # Join the line with spaces
    parsed = string.join(parsed)

    # Write the const char* kernel name line to the header file
    output.write("\t")
    output.write(parsed)

    # Seek to beginning of the kernel file
    fp.seek(0,0)

    # Skip over the comment block at the beginning of the kernel file
    # Write rest of file to the const char* block
    for line in fp.readlines():
        temp = string.split(line)
        if(len(temp) > 0 and (temp[0] == '/*' or temp[0] == '*' or temp[0] == '*/')):
            continue
        output.write("\t\t\"")
        output.write(string.rstrip(line))
        output.write("\\n\"\n")

    # End const char* block
    output.write("\t};\n\n")

    # Close kernel file
    fp.close()

    return

def gen_timestamp():
    # Get the time
    timestruct = time.localtime()

    # Begin generating the timestamp string
    timestamp = ""
    if(timestruct.tm_mon < 10):
        timestamp = "0" + str(timestruct.tm_mon) + "/"
    else:
        timestamp = str(timestruct.tm_mon) + "/"

    if(timestruct.tm_mday < 10):
        timestamp += "0" + str(timestruct.tm_mday) + "/"
    else:
        timestamp += str(timestruct.tm_mday) + "/"

    timestamp += str(timestruct.tm_year) + " "

    if(timestruct.tm_hour < 10):
        timestamp += "0" + str(timestruct.tm_hour) + ":"
    else:
        timestamp += str(timestruct.tm_hour) + ":"

    if(timestruct.tm_min < 10):
        timestamp += "0" + str(timestruct.tm_min) + ":"
    else:
        timestamp += str(timestruct.tm_min) + ":"

    if(timestruct.tm_sec < 10):
        timestamp += "0" + str(timestruct.tm_sec)
    else:
        timestamp += str(timestruct.tm_sec)

    return timestamp

def main():
    # Open file containing list of all selected kernels
    fp = open("file-list.txt","r")

    # Create file to write the generated header into
    output = open("opencl-domain-kernels.inl","w")

    # Add standard comment block to
    output.write("/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */\n"
                 "// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s\n"
                 "/* linbox/algorithms/opencl-domain.h\n"
                 " * Copyright (C) 2011 Matthew Wezowicz\n"
                 " *\n"
                 " * This library is free software; you can redistribute it and/or\n"
                 " * modify it under the terms of the GNU Lesser General Public\n"
                 " * License as published by the Free Software Foundation; either\n"
                 " * version 2 of the License, or (at your option) any later version.\n"
                 " *\n"
                 " * This library is distributed in the hope that it will be useful,\n"
                 " * but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
                 " * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU\n"
                 " * Lesser General Public License for more details.\n"
                 " *\n"
                 " * You should have received a copy of the GNU Lesser General Public\n"
                 " * License along with this library; if not, write to the\n"
                 " * Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,\n"
                 " * Boston, MA 02110-1301, USA.\n"
                 " */\n\n")

    # Add timestamp for generation
    output.write(" /*###---    Generated On     ---###*/\n")
    output.write(" /*###--- " + gen_timestamp() + " ---###*/\n\n")
    
    # Add header file lines
    output.write("#ifndef __LINBOX_opencl_matrix_domain_kernels_INL\n")
    output.write("#define __LINBOX_opencl_matrix_domain_kernels_INL\n\n")
    output.write("namespace LinBox{\n\n")

    # For every file in the file list generate a const char* and add to the header
    for file_name in fp.readlines():
        parse_kernel(string.rstrip(file_name),output)

    # Add end of header lines
    output.write("} /* end of namespace LinBox */\n\n")
    output.write("#endif /* __LINBOX_opencl_matrix_domain_kernels_INL */\n")

    # Close the files
    fp.close()
    output.close()

    return

main()
