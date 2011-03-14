# Copyright (c) the LinBox group
# This file is part of LinBox
# see COPYING for licence


fd := FileTools[Text][Open]("maple_version.txt", create=true, overwrite=true);
v:=substring(kernelopts(version),7..SearchText(".",kernelopts(version))-1);
FileTools[Text][WriteString](fd,convert(v,string));
