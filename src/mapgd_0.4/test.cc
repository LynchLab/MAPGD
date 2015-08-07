
class BaseClass{
public:
	virtual void * get(void)=0;
};

class Class {
public:
	const char* get(void) {return "HI!";};
};

int main (int argc, char *argv []){
	Class A;
};
