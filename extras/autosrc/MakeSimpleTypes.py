File=open("../../src/mapgd_0.4/datatypes/key_defs.txt")
header=open("../../src/mapgd_0.4/datatypes/simple_types.h", 'w')
cc=open("../../src/mapgd_0.4/datatypes/simple_types.cc", 'w')

basic_type={}
File.readline()

header.write("#ifndef _SIMPLE_TYPES_H_\n#define _SIMPLE_TYPES_H_\n\n")
header.write("#include \"key.h\"\n")
header.write("#include \"basic_types.h\"\n\n")
cc.write("#include \"simple_types.h\"\n\n")

for line in File:
	line=line.strip('\n').split('\t')
	class_name=line[0]
	class_id=line[1]
	class_type=line[2]
	class_desc=''.join(line[3:])
	if not(class_type.isupper() ):
		basic_type[class_name]=class_type
	else:
		if class_type=="SYNOM":
			continue
		print "static constexpr uint8_t "+class_name.lower()+"="+class_id+";"
		header.write("\n")
		header.write("class "+class_name+" :  public "+class_type+" {\n")
		header.write("public:\n")
		header.write("\t"+class_name+"();\n")
       		header.write("};\n")
		cc.write(class_name+"::"+class_name+"(void)\n")
		cc.write("{\n")
		cc.write("\tkeynum_="+class_name+";\n")
		cc.write("\tkeyname_="+class_id+";\n")
		cc.write("\tkeydesc_="+class_desc+";\n")
		cc.write("}\n")
header.write("#endif\n")
