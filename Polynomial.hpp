#ifndef _POLYNOMIAL_H
#define _POLYNOMIAL_H

#include <cstddef>
#include <iosfwd>
#include <functional>
#include <map>
#include <cmath>
#include <vector>
#include <algorithm>
#include "absvalue_wrapper.hpp"
#include "add_mult_identity.hpp"
#include <iostream>  // for Hentrich interface


void pause();

template<typename T>
class Polynomial
{
    typedef std::map<size_t, T, std::greater<size_t> > coefficientsMap;

    public:
        /* Constructors */
        // Default constructor
        Polynomial<T>();

        // Constructors for constant polynomials
        Polynomial<T>(T coefficient);

        // Default constructor with initial coefficient length provided
       // Polynomial<T>(const size_t coefficients);

        // Copy constructor from another polynom
        Polynomial<T>(const Polynomial<T>& poly);

        // MWH
        // Copy Constructor from a std::vector<T>
        Polynomial<T>(const std::vector<T>& coeffs);

        /* Destructor */
        ~Polynomial<T>();

        /* Operators */
        // Assignment operator
        const Polynomial<T>& operator = (const Polynomial<T>& poly);

        // Subscript getter operator to get an nth coefficient
        T operator [] (const size_t index) const;
        /*T operator [] (const size_t index) const;
        T& operator [] (const size_t index);*/

        // Arithmetic operators
        template<typename U>
        friend Polynomial<U> operator + (const Polynomial<U>& a, const Polynomial<U>& b);
        template<typename U>
        friend Polynomial<U> operator - (const Polynomial<U>& a, const Polynomial<U>& b);

        template<typename U>
        friend Polynomial<U> operator * (const Polynomial<U>& a, const Polynomial<U>& b);

        template<typename U>
        friend Polynomial<U> operator / (const Polynomial<U>& a, const Polynomial<U>& b);
        template<typename U>
        friend Polynomial<U> operator % (const Polynomial<U>& a, const Polynomial<U>& b);

         // equality, inequality
        template<typename U>
        friend bool operator == (const Polynomial<U>& a, const Polynomial<U>& b);

        template<typename U>
        friend bool operator != (const Polynomial<U>& a, const Polynomial<U>& b);

        // Ordering
        template<typename U>
        friend bool operator < (const Polynomial<U>& a, const Polynomial<U>& b);
        template<typename U>
        friend bool operator <= (const Polynomial<U>& a, const Polynomial<U>& b);
        template<typename U>
        friend bool operator > (const Polynomial<U>& a, const Polynomial<U>& b);
        template<typename U>
        friend bool operator >= (const Polynomial<U>& a, const Polynomial<U>& b);

        // Printing operator
        template<typename U>
        friend std::ostream& operator << (std::ostream& o, const Polynomial<U>& poly);

        /* Member functions */
        // Get the degree
        size_t degree() const;

        // Get the leading coefficient (shortcut to the degree-th member's)
        T leadingCoefficient() const
        {
            return this->getMember(this->degree());
        }

        // Get the nth coefficient
        T getMember(const size_t index) const;

        // Set the nth coefficient
        void setMember(const size_t index, const T coefficient);

        // Calculate the polynomial function's value for variable 't'
        T at(const T t) const;

        // Algebraic derivative (prime)
        Polynomial<T> derive() const;

        // Arithmetical methods
        void add(const Polynomial<T>& poly);
        void subtract(const Polynomial<T>& poly);
        void multiply(const Polynomial<T>& poly);
        bool divide(const Polynomial<T>& divisor, Polynomial<T>& quotient, Polynomial<T>& remainder) const;

        // Equality
        bool equals(const Polynomial<T>& poly) const;

        // Get if the polynomial is a nullpolynomial (that is: every coefficient is zero)
        bool isNull() const;

        // Get if the polynomials is a constant polynomials
        // (that is: a number... only one coefficient on the zeroth power, and it is not zero)
        bool isConstant() const;

    private:
        // Array for coefficients
        // (First part of the pair is the power of the indeterminate,
        // the second part is the coefficient.)
        typename Polynomial<T>::coefficientsMap m_coefficients;

        // Degree of the polynomial
        size_t m_degree;

        // Determine if there is an indeterminate of the given power index
        bool hasMember(const size_t index) const;

        // Internal cleanup function.
        void _performCleanup();
};

template<typename T>
Polynomial<T>::Polynomial()
{
    //cout << "Polynomial initialized." << endl;
    this->_performCleanup(); // Make the object into default state
    //cout << "Polynomial initialized." << endl;
    // Noop.
}

template<typename T>
Polynomial<T>::Polynomial(T coefficient)
{
    //cout << "Polynomial initialized as constant " << coefficient << endl;
    this->setMember(0, coefficient);
    this->_performCleanup();
}

/* template<typename T>
Polynomial<T>::Polynomial(const size_t coefficients)
{
        this->m_coefficients.reserve(coefficients);
}  */

template<typename T>
Polynomial<T>::Polynomial(const Polynomial<T>& poly)
{
    //cout << "Polynomial copied from " << poly << endl;
    // Copy the map.
    this->m_coefficients = poly.m_coefficients;
    this->_performCleanup();
}

// MWH
// Copy Constructor from a std::vector<T>
template<typename T>
Polynomial<T>::Polynomial(const std::vector<T>& coeffs)  {
    //cout << "Polynomial copied from " << coeffs << endl;
    // Copy the elements.
   // this->m_coefficients = <=== 
   Polynomial<T> extracted_coefficients((size_t)(coeffs.size()));
   for (size_t power = 0; power < coeffs.size(); ++power) {
        extracted_coefficients.setMember( power, (T)(coeffs.at(power)) );
   }
   this->_performCleanup();
}


#include <sstream>

template<typename T>
Polynomial<T>::~Polynomial()
{
    /*stringstream ss;
    ss << *this;

    cout << "The polynomial " << ss.str() << " was destructed." << endl;*/
}

template<typename T>
const Polynomial<T>& Polynomial<T>::operator=(const Polynomial<T>& poly)
{
    this->m_coefficients = poly.m_coefficients;
    this->_performCleanup();
    return *this;
}

template<typename T>
T Polynomial<T>::operator[] (const size_t index) const
{
    return getMember(index);
}

template<typename T>
bool operator==(const Polynomial<T>& a, const Polynomial<T>& b)
{
    return a.equals(b);
}

template<typename T>
bool operator!=(const Polynomial<T>& a, const Polynomial<T>& b)
{
    return !a.equals(b);
}

template<typename T>
std::ostream& operator << (std::ostream& o, const Polynomial<T>& poly)
{
    if (poly.isNull())
        return o << "nullPolynomial";
    else if (poly.isConstant())
        return o << poly.getMember(0); // << " (constant)";

    bool firstCoeff = true;
    for (typename Polynomial<T>::coefficientsMap::const_iterator cit = poly.m_coefficients.cbegin();
        cit != poly.m_coefficients.cend(); ++cit)
    {
        // cit->first: power
        // cit->second: coefficient

        // Omit the members with 0 coefficient
        if (cit->second == id_additive<T>::value) continue;

        if (!firstCoeff)
            o << " ";

        // Print the signum of the coefficient
        // (so it would neatly look as -5x^2 + 3x - 2)
        if (cit->second < id_additive<T>::value)
            o << "- ";
        else if (cit->second > id_additive<T>::value && !firstCoeff)
            o << "+ ";

        // Only print the coefficient if it is not 1. The 1 multiplier as it's multiplicative identity can be omitted.
        // Because of templating, we use a custom implementation here.
        //
        // But this should only happen for the non-constant member...
        if (cit->first != 0)
        {
            if (id_multiplicative_exists<T>::value)
            {
                if (abs_value<T>::known && abs_value<T>::abs(cit->second) != id_multiplicative<T>::value)
                    o << abs_value<T>::abs(cit->second);
                else if (!abs_value<T>::known && cit->second != id_multiplicative<T>::value)
                    o << cit->second;
            }
            else
                o << cit->second;
        }
        else
        {
            // If the number is negative, print only its absolute value (the - was printed earlier)
            if (abs_value<T>::known && cit->second < abs_value<T>::abs(cit->second))
                o << abs_value<T>::abs(cit->second);
            else
                o << cit->second;
        }

        // Don't print the power for the first-order member and don't print x for the constant member
        if (cit->first > 1)
            o << "x^" << cit->first;
        if (cit->first == 1)
            o << "x";

        if (firstCoeff) firstCoeff = false;
    }

    return o;
}

template<typename T>
size_t Polynomial<T>::degree() const
{
    return this->m_degree;
}

template<typename T>
T Polynomial<T>::getMember(const size_t index) const
{
    if (this->hasMember(index))
        return this->m_coefficients.at(index);
    else
        return id_additive<T>::value; // Non-existant members are just not stored: they are mathematically there with 0 coefficient.
}

template<typename T>
void Polynomial<T>::setMember(const size_t index, const T coefficient)
{
    this->m_coefficients[index] = coefficient; // Add or reassign a member
    this->_performCleanup();
}

template<typename T>
T Polynomial<T>::at(const T t) const
{
    // Using Horner's rule:
    // Consider the polynomial x^4 - 3x^3 + 0x^2 - x + 6 at -1.
    //
    //     1  -3   0   1   6
    //   |------------------
    // -1| 1  -4   4  -3   9
    //       (^)
    // The first element of the table (the LC) is copied down.
    // After that, every element is calculated as follows: the previous element * the number to the left (the -1)
    // + the element above. (So at (^) we have 1 * (-1) + (-3) = -1 - 3 = -4, etc.

    T result = this->leadingCoefficient();

    //                                                 V have to make sure we don't underflow.
    for (size_t power = this->degree() - 1; power >= 0 && power <= this->degree(); --power)
        result = result * t + this->getMember(power);

    return result;
}

template<typename T>
Polynomial<T> Polynomial<T>::derive() const
{
    if (this->isNull() || this->isConstant())
        return Polynomial<T>();

    // If we consider the polynomial as f0 + f1x + f2x^2 + ..., the
    //  algebraic derivative (prime) is f1 + 2*f2x + 3*f3x^2 + ...
    // (Where 2*f2 actually means f2 + f2, as f2 is just a variable of type T)
    Polynomial<T> derivative;
    for (size_t power = 0; power < this->degree(); ++power)
    {
        // The i-th coefficient of the derivate is (i + 1)* the (i+1)th coefficient
        T curr_coeff = this->getMember(power + 1);

        // Don't run the loop for coefficients which are 0
        if (curr_coeff == id_additive<T>::value) continue;

        for (size_t i = 1; i < power + 1; ++i)
            curr_coeff += this->getMember(power + 1);

        derivative.setMember(power, curr_coeff);

        // The coefficient for the degree-th member is zero, the prime "eats" that part.
        // Also, the constant part is eaten too.
    }

    return derivative;
}

template<typename T>
void Polynomial<T>::add(const Polynomial<T>& poly)
{
    //          2x^2 + 3x + 1
    // (+) x^3 + x^2 - 2x - 2
    // ----------------------
    //     x^3 +3x^2 +  x - 1
    //
    // So we have to member-by-member add the coefficients to get the added polynomial

    for (size_t p = 0; p <= std::max(this->degree(), poly.degree()); ++p)
        this->setMember(p, this->getMember(p) + poly.getMember(p));
}

template<typename T>
void Polynomial<T>::subtract(const Polynomial<T>& poly)
{
    // Subtraction works just as so
    for (size_t p = 0; p <= std::max(this->degree(), poly.degree()); ++p)
        this->setMember(p, this->getMember(p) - poly.getMember(p));
}

template<typename T>
void Polynomial<T>::multiply(const Polynomial<T>& poly)
{
    // If either is a nullpolynomial, the multiple is trivially a nullpolynomial
    if (this->isNull() || poly.isNull())
        this->m_coefficients.clear();
    // If both is a constant, the multiple is trivially the constants' multiple
    else if (this->isConstant() && poly.isConstant())
    {
        T multiple = this->leadingCoefficient() * poly.leadingCoefficient();

        this->m_coefficients.clear();
        this->setMember(0, multiple);
    }
    // If either polynomials is a constant one (but not both), we can still ease the multiplication
    else if (this->isConstant() || poly.isConstant())
    {
        // Select which polynomial is the constant and the more complex one;
        const Polynomial* constantOne = (this->isConstant() ? this : &poly);
        const Polynomial* complexOne = (!this->isConstant() ? this : &poly);
        Polynomial<T> multiple;

        // Only multiply the coefficients of the complex polynomial by the given constant, storing it in this
        T constant = constantOne->leadingCoefficient();
        for (size_t i = 0; i <= complexOne->degree(); ++i)
        {
            T coeff = complexOne->getMember(i);
            if (coeff == id_additive<T>::value) continue;

            multiple.setMember(i, coeff * constant);
        }

        // Let the multiplied polynomial be the current one
        *this = multiple;
    }
    // If both polynomials are complex ones, we do the multiplication by hand
    else
    {
        std::vector<T> multi_coefficients;
        multi_coefficients.resize(this->degree() + poly.degree() + 1); // Create space for the coefficients

        for (size_t i = 0; i <= this->degree(); ++i)
        {
            T left = this->getMember(i);
            if (left == id_additive<T>::value) continue; // 0 * anything = 0

            for (size_t j = 0; j <= poly.degree(); ++j)
            {
                T right = poly.getMember(j);
                if (right == id_additive<T>::value) continue; // anything * 0 = 0

                // Basically you have to multiply every member with every member...
                multi_coefficients.at(i+j) += left * right;
            }
        }

        // Set the calculated coefficients
        this->m_coefficients.clear();
        for (size_t i = 0; i < multi_coefficients.size(); ++i)
            this->setMember(i, multi_coefficients.at(i));
    }
}

template<typename T>
bool Polynomial<T>::divide(const Polynomial<T>& divisor, 
                           Polynomial<T>& quotient, 
                           Polynomial<T>& remainder) const
{
    if (divisor.isNull())
        return false;

    if (divisor.degree() > this->degree())
    {
        remainder = *this;
        Polynomial<T> quotient_null;
        quotient = quotient_null;
        return true;
    }
    else
    // Degree of divisor is smaller or equal than dividend
    {
        Polynomial<T> dividend = *this; // Save out the current instance
        // f : g = q
        // f % g = r
        // f: dividend (this), g: divisor, q: quotient, r: remainder
        // The dividend is only initially 'this', it gets consumed as the division happens.

        size_t quotient_max_degree = dividend.degree() - divisor.degree();
        size_t quotient_member_degree = quotient_max_degree;

        // We can only divide further if there is something to divide.
        while (quotient_member_degree >= 0 && 
         quotient_member_degree <= quotient_max_degree && !dividend.isNull())
        {
            std::cout << "-----------------------------" << std::endl;
            std::cout << "We are in a division loop. \n"
                      << "Current degree of the quotient: " 
                      << quotient_member_degree << std::endl;

            std::cout << "(Maximum possible degree of quotient: " 
                      << quotient_max_degree << ")" << std::endl;
            std::cout << "Dividend: " << dividend << std::endl;

            // Subtract the divisor's LC from the dividend's LC
            // (and so the powers), so that we get the LC of the quotient.
            std::cout << "degree of dividend: " << dividend.degree() << std::endl;
            std::cout << "leading coefficient of dividend:  " << dividend.leadingCoefficient() 
                      << std::endl;

            // After subtracting the LCs, we get the quotient's LC.
            // Using polynomial long division, we now multiply the quotient's LC by the divisor (-> getting a polynomial)
            // and thus, we subtract that from the dividend... and loop this whole shit.)

            // If the current member degree underflows, thus reaches a values higher than the maximum degree
            // We terminate the whole cycle because we reached the remainder.
            if (quotient_member_degree > quotient_max_degree)
                break;

            // Set the current member of the quotient to the quotient of the coefficients
            quotient.setMember(quotient_member_degree, dividend.leadingCoefficient() / divisor.leadingCoefficient());

            Polynomial<T> multiplier;
            multiplier.setMember(quotient_member_degree, quotient.getMember(quotient_member_degree));

            std::cout << "Divisor: " << divisor << std::endl;
            Polynomial<T> inner_multiple = divisor * multiplier;

            dividend = dividend - inner_multiple;

            std::cout << "Current quotient (across division operation): " << quotient << std::endl;
            std::cout << "Current multiplier: " << multiplier << std::endl;
            std::cout << "Multiplied dividend: (this will be subtracted) " << inner_multiple << std::endl;
            std::cout << "New dividend (after subtraction): " << dividend << std::endl;

            quotient_member_degree = dividend.degree() - divisor.degree();
            pause();
        }

        // When the loop reaches its terminus, we divided everything we could.
        // Anything that remained in the dividend (after subtraction) 
        // is actually the remainder.

        remainder = dividend; // So assign it into its proper place.

        return true;
    }
}

/* Helper functions to equate two values... needed for floating-point arithmetics. */
template<typename T>
bool equate(const T& a, const T& b)
{
    return a == b;
}

template<typename T>
bool float_equate(const T& a, const T& b, T epsilon)
{
    return fabs(a - b) < epsilon;
}

#ifndef _INC_FLOAT
#include <cfloat>
#endif // _INC_FLOAT

template<> bool equate(const float& a, const float& b) { return float_equate(a, b, FLT_EPSILON); }
template<> bool equate(const double& a, const double& b) { return float_equate(a, b, DBL_EPSILON); }
template<> bool equate(const long double& a, const long double& b) { return float_equate(a, b, (long double)LDBL_EPSILON); }
/* Helper functions over. */

template<typename T>
bool Polynomial<T>::equals(const Polynomial<T>& poly) const
{
    // Two polynomials are equal if their degree is equal and every coefficient is equal for every member.
    /* We could do it this way, but it's safer if we check without assuming the invariant for now...
    return this->m_coefficients.size() == poly.m_coefficients.size() &&
        equal(this->m_coefficients.cbegin(), this->m_coefficients.cend(), poly.m_coefficients.cbegin());
    */

    // Two polynomials equal if their degree is equal and every coefficient is equal for every member
    if (this->degree() != poly.degree())
        return false;

    bool is_equal = true; // Assume that it is equal
    for (size_t i = 0; i <= this->degree() && is_equal == true; ++i) // The sizes are now equal
        is_equal &= equate(this->getMember(i), poly.getMember(i));

    return is_equal;
}

template<typename T>
bool Polynomial<T>::isNull() const
{
    // A polynomial is a null-polynomial if every coefficient is zero.
    // Utilising the invariant, if we are nullpolynomial, 
    // the degree is zero and that 0th coefficient is zero.
    return (this->degree() == 0 && this->getMember(0) == id_additive<T>::value);
}

template<typename T>
bool Polynomial<T>::isConstant() const
{
    // The polynomial is a constant one if there is only a constant (0th power) member
    return (this->degree() == 0 && this->getMember(0) != id_additive<T>::value);
}

template<typename T>
void Polynomial<T>::_performCleanup()
{
    // Cleanup consists of removing the 0 coefficient parts from the map
    std::vector<size_t> removePowers;
    removePowers.reserve(this->m_coefficients.size());

    // And keeping the degree of the object the degree of the polynomial
    size_t maxPower = 0;

    // Get the 0 coefficient elements from the map
    for (typename Polynomial<T>::coefficientsMap::const_iterator cit = this->m_coefficients.cbegin();
        cit != this->m_coefficients.cend(); ++cit)
    {
        if (cit->second == id_additive<T>::value)
            removePowers.push_back(cit->first);
        else
            if (cit->first > maxPower)
                maxPower = cit->first;
    }

    // Remove the elements
    for (std::vector<size_t>::const_iterator cit = removePowers.cbegin(); cit != removePowers.cend(); ++cit)
        this->m_coefficients.erase(*cit);

    this->m_degree = maxPower;
}

template<typename T>
bool Polynomial<T>::hasMember(const size_t index) const
{
    // This is basically the standard way of finding if the given power exists
    // The iterator becomes the end of the collection if the value is not found
    return this->m_coefficients.find(index) != this->m_coefficients.end();
}

template<typename T>
Polynomial<T> operator+(const Polynomial<T>& a, const Polynomial<T>& b)
{
    Polynomial<T> ret = a;
    ret.add(b);

    return ret;
}

template<typename T>
Polynomial<T> operator-(const Polynomial<T>& a, const Polynomial<T>& b)
{
    Polynomial<T> ret = a;
    ret.subtract(b);

    return ret;
}

template<typename T>
Polynomial<T> operator*(const Polynomial<T>& a, const Polynomial<T>& b)
{
    Polynomial<T> ret = a;
    ret.multiply(b);

    return ret;
}

template<typename T>
Polynomial<T> operator/(const Polynomial<T>& a, const Polynomial<T>& b)
{
    Polynomial<T> q; // quotient
    Polynomial<T> r; // remainder

    a.divide(b, q, r);

    return q;
}

template<typename T>
Polynomial<T> operator%(const Polynomial<T>& a, const Polynomial<T>& b)
{
    Polynomial<T> q; // quotient
    Polynomial<T> r; // remainder

    a.divide(b, q, r);

    return r;
}

template<typename T>
std::vector<T> extract_coefficient_vector(const Polynomial<T>& poly)    {
    std::vector<T> coefficient_vector(poly.degree() + 1);
    for (auto&& i : coefficient_vector)  i = T(0);
    for (size_t k = 0; k <= poly.degree(); ++k) {
          coefficient_vector.at(k) = poly.getMember(k);
    }

    return coefficient_vector;
}


std::string Ntab(int n)  {
  std::string s = "";
  for (int i = 0; i < n; ++i)   s += '\t';
  return s;
}



/*   It turns out that size_t is like an unsigned long integer, 
     and I have this size_t j decrementing while j <= 0; but when 
     j decreases by 1 (when current value of j is 0), there is 
     no representation for -1 in data type size_t 
     (from some C-style library <cstddef> ) 
      hence 0 - 1 = 18446744073709551615.   I changed n-1 to n
      and j >= 0 to j > 0, then getMember(j - 1)  
      */
template <typename T>
std::pair< Polynomial<T>, T > 
QR_sDiv(const Polynomial<T>& P, const T& c) 
{
  size_t n(P.degree());
  std::cout << "\nP(x) = " << P << '\n';

  Polynomial<T> R_local(P);  // create Polynomial from P
  // then just initialize elements 0 to (n-1) of R_local.m_coefficients to 0:
  for (size_t i = 0; i < n; i++) { R_local.setMember(i, T(0));  }

   std::cout << "\n\n " << c << "  |" << Ntab(n+1) << "| f(0)\n";
   std::cout << "-----------------------------------------------------";
   std::cout << "\n\t |  " << P.getMember(n) << " \t ";
   for (size_t j = n; j > 0; j--)  {
       std::cout << P.getMember(j - 1) << " \t ";
    }

   std::cout << "\n\n\t |" << Ntab(n+1) << "| f(" << c << ")\n";
   std::cout << "----------------------------------------------------";
  
   T y = P.getMember(n);
   //R_local.setMember(n, y);
   T x = c;
   std::cout << "\n\t | " << y << " \t ";
      for (size_t j = n; j > 0; j--)  {
         y = y*x + P.getMember(j - 1);              // SYNTHETIC DIVISION 
         R_local.setMember(j - 1, y);
         std::cout << R_local.getMember(j - 1) << " \t ";
         //pause();
        // std::cout << "\nj = " << j << '\n';
      }
    std::cout << "\n\n";

   Polynomial<T> Quotient;
   Quotient.setMember(n - 1, P.getMember(n));

    for (size_t j = n-1; j > 0; j--)  {
         Quotient.setMember(j - 1, R_local.getMember(j));
    }

    return std::make_pair(Quotient, R_local.getMember(0));

  }  


// The 'phi' functions of polynomials (in the Euclidean ring order) is their degree
template<typename T>
bool operator<(const Polynomial<T>& a, const Polynomial<T>& b)
{
    return a.degree() < b.degree();
}

template<typename T>
bool operator<=(const Polynomial<T>& a, const Polynomial<T>& b)
{
    return a.degree() <= b.degree();
}

template<typename T>
bool operator>(const Polynomial<T>& a, const Polynomial<T>& b)
{
    return a.degree() > b.degree();
}

template<typename T>
bool operator>=(const Polynomial<T>& a, const Polynomial<T>& b)
{
    return a.degree >= b.degree();
}


#ifdef _ADD_MULT_IDENTITY_H
// The additive and multiplicative identity of a polynomial of type T is
// a constant polynomial with LC of add/mult id of type T.
template<typename T> struct id_multiplicative_exists<Polynomial<T>> : id_multiplicative_exists<T>{};
template<typename T> struct id_multiplicative<Polynomial<T>> { static Polynomial<T> const value; };
template<typename T> Polynomial<T> const id_multiplicative<Polynomial<T>>::value = Polynomial<T>(id_multiplicative<T>::value);

template<typename T> struct id_additive_exists<Polynomial<T>> : id_additive_exists<T>{};
template<typename T> struct id_additive<Polynomial<T>> { static Polynomial<T> const value; };
template<typename T> Polynomial<T> const id_additive<Polynomial<T>>::value = Polynomial<T>(id_additive<T>::value);
#endif // _ADD_MULT_IDENTITY_H

#endif // _POLYNOMIAL_H

/*
void pause()  
{
    std::cin.clear();
    std::cin.ignore(numeric_limits<streamsize>::max(), '\n');
    std::string dummy;
    std::cout << "\t Press ENTER [ <-----||| ] to continue ... \n";
    std::getline(std::cin, dummy);
}
*/