#!/usr/bin/perl

$class_c_name = $ARGV[0];
@class_comp = split (/_/, $class_c_name);

$class_c_uc_name = uc $class_c_name;
$class_cpp_name = join ('', map (ucfirst, @class_comp));
$class_filename_base = join ('-', @class_comp);

$parent_c_name = $ARGV[1];
@parent_comp = split (/_/, $parent_c_name);

$parent_c_uc_name = uc $parent_c_name;
$parent_cpp_name = join ('', map (ucfirst, @parent_comp));
$parent_filename_base = join ('-', @parent_comp);

write_file ("template.c", "$class_filename_base.c");
write_file ("template.h", "$class_filename_base.h");

sub write_file
  {
    my ($input, $output) = @_;

    open FILE_INPUT, "<$input";
    open FILE_OUTPUT, ">$output";

    while (<FILE_INPUT>)
      {
	next if /\#include \"parent-class-name.h\"/ &&
	  ($parent_c_name =~ /gtk|gnome/);

	s/WidgetClassName/$class_cpp_name/g;
	s/WIDGET_CLASS_NAME/$class_c_uc_name/g;
	s/widget_class_name/$class_c_name/g;
	s/widget-class-name/$class_filename_base/g;

	s/ParentClassName/$parent_cpp_name/g;
	s/PARENT_CLASS_NAME/$parent_c_uc_name/g;
	s/parent_class_name/$parent_c_name/g;
	s/parent-class-name/$parent_filename_base/g;

	print FILE_OUTPUT;
      }
  }
