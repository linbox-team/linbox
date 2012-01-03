import string

def parse_kernel(file_name, output):
    fp = open(file_name)

    parsed = []

    for line in fp.readlines():
        line = string.split(line)
        if(len(line) > 0 and line[0] == '__kernel'):
            parsed.append("const char*")
            temp = string.split(line[2],'(')
            parsed.append(temp[0])
            parsed.append("= {\n")
            break

    parsed = string.join(parsed)

    output.write("\t")
    output.write(parsed)

    fp.seek(0,0)

    for line in fp.readlines():
        temp = string.split(line)
        if(len(temp) > 0 and (temp[0] == '/*' or temp[0] == '*' or temp[0] == '*/')):
            continue
        output.write("\t\t\"")
        output.write(string.rstrip(line))
        output.write("\\n\"\n")

    output.write("\t};\n\n")

    fp.close()

    return

def main():
    fp = open("file-list.txt","r")
    output = open("opencl-domain-kernels.inl","w")

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
    output.write("#ifndef __LINBOX_opencl_matrix_domain_kernels_INL\n")
    output.write("#define __LINBOX_opencl_matrix_domain_kernels_INL\n\n")
    output.write("namespace LinBox{\n\n")

    for file_name in fp.readlines():
        parse_kernel(string.rstrip(file_name),output)

    output.write("} /* end of namespace LinBox */\n\n")
    output.write("#endif /* __LINBOX_opencl_matrix_domain_kernels_INL */")

    fp.close()
    output.close()

    return

main()
