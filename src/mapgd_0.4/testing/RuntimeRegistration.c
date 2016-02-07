    #include <iostream>
    #include <string>
    #include <map>

    class D {
	public:
	virtual D(){};
	virtual ~D(){};
	protected:
    }

    class FactoryBase
    {
    protected:
        FactoryBase (const std::string &name, D* item)
        {
            m_factory_items [name] = item;
        }

    public:
        virtual ~FactoryBase ()
        {
        }

        static D * Create ()
        {
            return (m_factory_items [T::Name]->Create () );
        }

    private:
        virtual D * Create () = 0;

    private:
        static std::map <const std::string, D *> m_factory_items;
    };

    std::map <const std::string, D *>
        FactoryBase::m_factory_items;

    template <class T>
        class FactoryItem : public FactoryBase
    {
    public:
        FactoryItem () :
            FactoryBase (T::Name)
        {
            std::cout << "Registering class: " << T::Name << std::endl;
        }

        virtual ~FactoryItem ()
        {
        }

    private:
        virtual void *Create ()
        {
            return new T;
        }
    };

    class A
    {
    public:
        A ()
        {
            std::cout << "Creating A" << std::endl;
        }

        virtual ~A ()
        {
            std::cout << "Deleting A" << std::endl;
        }

        static const std::string
            Name;

    private:
        static FactoryItem <A>
            m_registration;
    };

    const std::string
        A::Name ("A");

    FactoryItem <A>
        A::m_registration;

    class B
    {
    public:
        B ()
        {
            std::cout << "Creating B" << std::endl;
        }

        virtual ~B ()
        {
            std::cout << "Deleting B" << std::endl;
        }

        static const std::string
            Name;

    private:
        static FactoryItem <B>
            m_registration;
    };

    const std::string
        B::Name ("B");

    FactoryItem <B>
        B::m_registration;

    int main (int argc, char *argv [])
    {
        A
            *item_a = FactoryBase::Create <A> ();

        B
            *item_b = FactoryBase::Create <B> ();

        delete item_a;
        delete item_b;
    }
