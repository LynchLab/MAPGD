#ifndef REF_H_
#define REF_H_

class ref: public data {
private:	
	gt_t value;									//!< an array to allow sorted access to quartets.
public:
	size_t size(void){return sizeof(gt_t);};						//!< all class of type data need to decleare a size function that returns the total size of ...
	const std::type_info& type(void ){return typeid(*this);};
	ref (){
		key_=new key("REF", 21, "the ancestraol allele at this locus", this);
	};
};
#endif
