/* a class to provide fast retreval of ? */ 

class index 
{
private:
	fastmap <std::string> id0_;
public:
	uint64_t get_offset();	
}

template <class T> class fastmap
{
private:
        std::map <T, count_t> id_t_;		//
        std::list <T> id_;			//

        uint16_t lastid_;			//initilize to 0-1;
        T * lastid_t_;				//initilize to "";
	unint16_t maxid_;
public:
	uint16_t encodeid(const T &);		//
	const T & decodeid(?);			//
}

template <class T> class lookup
{
private:
        uint16_t * id_;				//
        T * id_t_;
public:
	uint16_t encodeid(const T &);		//
	const T & decodeid(uint16_t &);		//
}
