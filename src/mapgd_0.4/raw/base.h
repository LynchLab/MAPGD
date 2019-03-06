#ifndef BASE_H_
#define BASE_H_

#include <iostream>
#include "typedef.h"
#include <cfloat>
#include <iomanip>

// PLEASE LIMIT LINE LENGTH TO 79 CHARACTERS----------------------------------/

/// A class converts human readable bases to bit flags
/*!
 * Uses a bit flags to represent the IUPCA encoded bases. Bits 1,2,3 and 4 
 * represent A, C, G and T. An ambiguous base will set all applicable flags. 
 */
class Base { 
private:
    /// indicates the site has been masked.
    /*! Some calculations check if a site is masked before use.
    */
    bool masked;
public:
    Base();

    /// constructor from IUPCA base code. 
    Base(const char &);

    /// constructor from Base. 
    Base(const gt_t &);

    /// Represents a single base.
    gt_t base;

    /// use the << operator to write Allele.
    friend std::ostream& operator << (std::ostream&, const Base&);    
    /// use the >> operator to read Allele.
    friend std::istream& operator >> (std::istream&, Base&);    

    static char btoc(const gt_t &);
    static gt_t ctob(const char &);

    //Logical comparisons are 'could be' 'must be' 'can't be'
    bool operator==(const Base& rhs) const
    {
        return base == rhs.base;
    }

};

#endif
