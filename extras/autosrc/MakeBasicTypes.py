File=open("../src/mapgd_0.4/datatypes/key_defs.txt")
header=open("../src/mapgd_0.4/datatypes/basic_types.h", 'w')
cc=open("../src/mapgd_0.4/datatypes/basic_types.cc", 'w')

basic_type={}
File.readline()

header.write("#ifndef _BASIC_TYPES_H_\n#define _BASIC_TYPES_H_\n\n")
cc.write("#include \"basic_types.h\"\n\n")
for line in File:
	line=line.strip('\n').split('\t')
	class_name=line[0]
	class_id=line[1]
	class_type=line[2]
	class_desc=''.join(line[3:])
	if not(class_type.isupper() ):
		if (float(class_id)<4):
			continue
		header.write("\n")
		header.write("class "+class_name+" :  public data {\n")
		header.write("public:\n")
		header.write("\tsize_t sizeb(void){return sizeof("+class_type+");};\n")
		header.write("\tsize_t sizet(void){return 7;};\n")
       	 	header.write("\t"+class_name+"();\n")
		header.write("}\n")

		cc.write(class_name+"::"+class_name+"(void)\n")
		cc.write("\t"+class_name+"();\n")
		cc.write("}\n\n")
header.write("#endif\n")
