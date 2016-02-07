#include <iostream>
#include <string>
#include <map>


class Base
{
public:
	Base () {};
        static Base * Create (const std::string &Name)
        {
 		return m_factory_items[Name]();
        }

	virtual ~Base(){};
	virtual void shout(void)=0;
	static std::map <std::string, Base *(*)(void)> m_factory_items;
private:
};

std::map <std::string, Base *(*)(void)> Base::m_factory_items;

class Registration   
{
public:
	Registration (const std::string &str, Base*(*fn)(void) ){
		if(Base::m_factory_items.find(str)!=Base::m_factory_items.end() ){
			std::cerr << "A class of name " << str << " already exist. Access to the previous class is now lost\n"; 
		}
		Base::m_factory_items[str]=fn;
		name_=str;
	}
	~Registration (void){
		Base::m_factory_items.erase(name_);
	}
private:
	std::string name_;
};


class A : public Base 
{
private:
	static const std::string name;
	static Base* Create(){
		return new A;
	}
	static Registration registered;
public:
	A()
	{
		std::cout << "Made Class A\n";		
	}
	~A()
	{
		std::cout << "Killed Class A\n";
	}
	void shout()
	{
		std::cout << "A!!\n";
	}
};

Registration A::registered("A", A::Create);

class B : public Base 
{
private:
	static Base* Create(){
		return new B;
	}
	static Registration registered1;
	static Registration registered2;
public:
	B()
	{
		std::cout << "Made Class B\n";		
	}
	~B()
	{
		std::cout << "Killed Class B\n";
	}
	void shout()
	{
		std::cout << "B!!\n";
	}
};
Registration B::registered1("B", B::Create);
Registration B::registered2("A", B::Create);


int main(int argc, char *argv[]){
	Base *a=Base::Create("A");
	Base *b=Base::Create("B");
	a->shout();
	b->shout();
	delete a;
	delete b;
};
