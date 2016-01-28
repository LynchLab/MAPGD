class Data {
public:
	virtual Data(std::vector <std::string> &);
	virtual static const file_name;	
        virtual std::string header(void) const;
        virtual std::string table_name(void) const;	//!< The destination table in the Db.
        virtual size_t size(void) const;        	//!< The size of the ? in bytes.
};
