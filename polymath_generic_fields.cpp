
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 *   Algebra++ :=    
 *            [abstract=="modern"]  A N A L Y S I S
 *            with  G E N E R I C  P R O G R A M M I N G  
 *            using Modern C++ 
 * 
 * 
 *   polycalc++:   polymath_generic_fields.cpp
 * 
 *   COMPILE:  g++ -g -Wall polymath_generic_fields.cpp -std=c++14 -o polycalc++
 * 
 *   Format to Print:
 * 
 *   enscript -1rG --line-numbers --portrait -p polymath_generic_fields.ps 
 *            --highlight=cpp -c polymath_generic_fields.cpp
 *           
 *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  * 
 *   This version uses an obscure but novel Polynom-Nom polynomial library 
 *   for its educational value.
 *   https://github.com/whisperity/polynomnom
 */

#include "Polynomial.hpp"
#include "EuclideanAlgorithm.hpp"
#include "Complex.hpp"    // for Polynomial over Complex Field: C
#include "Rational.hpp"   // for Polynomial over Rational Field: Q
#include <complex>      // for private member inner _complex
#include <iostream>
//--------------------------------------------------------------------
// added September 2019: Synthetic Division (see notes below)
//-------------------------------------------------------------------
#include <vector>
#include <utility>   // added by MW Hentrich 2019.02.26 for STL pair
                     // I want division to return quotient and remainder
// I will also add member function template <class K> Polynomial::QR_sDiv(K) 
// to display "synthetic division process"

#include <string>    // for Ntab(int) used to display synthetic division
                     // in QR_sDiv

#include <istream>  // to use cin.get() so as to use function "pause"

using namespace std;

void pause()  
{
    std::cin.clear();
    std::cin.ignore(numeric_limits<streamsize>::max(), '\n');
    std::string dummy;
    std::cout << "\t Press ENTER [ <-----||| ] to continue ... \n";
    std::getline(std::cin, dummy);
}

template<class K>
void getValues(Polynomial<K>& P1, Polynomial<K>& P2);

template<class K>
void getValue(Polynomial<K>& P);

template<class K>
void run(K dummy, char status, char op);

char choose_Field();
char choose_Operation();

template<class K>     // (see notes at end of Polynomial.h)
using Quotient_Remainder = std::pair< Polynomial<K>, K >;

int main() {
  char operationType, F;  // where F represents Field

  do    {
    F = choose_Field();  // Initialize Field 
    if ( (F == 'x') || (F == 'X') ) break;

    operationType = choose_Operation();
    if ( (operationType == 'x') || (operationType == 'X') ) break;
    
    if ( (F == 'Q')  ||  (F == 'q') )   run(Rational(1), F, operationType);
    if ( (F == 'R')  ||  (F == 'r') )   run(1.0, F, operationType);
    if ( (F == 'C')  ||  (F == 'c') )   run(Complex(1.0), F, operationType);

    std::cout << std::endl;
    pause();
    
  } while (operationType != 'x');
  std::cout << std::endl;
}  // end main


char choose_Field()  {
   char field = 'R';
   std::cout << "\n--------------------------------------------------------------\n";
   std::cout << "Polynomials over which field?"<< std::endl;
   std::cout << "\n--------------------------------------------------------------\n";
   std::cout << "   Q (Polynomial over the field of rational numbers x = a/b)\n"
             << "   R (Polynomial over the field of real numbers x = a.)\n"
             << "   C (Polynomial over the field of complex numbers x = a + i b)\n"
             << "\n--------------------------------------------------------------\n"
             << "    x (to exit)\n"
             << "\n--------------------------------------------------------------\n"
             << "This polynomial is over :  ";

 std::cin >> field;
 std::cout << std::endl;
 return field;
}


char choose_Operation()  {

    char op = 'x';

		 std::cout << "\n--------------------------------------------------------\n";
		 std::cout << "Choose the operation you would like to perform:"<< std::endl;
		 std::cout << "\n--------------------------------------------------------\n";
		 std::cout << "	+ (addition)\n"
		 		  << "	- (subtraction)\n"
				  << "	* (multiplication)\n"
				  << "        / (division)\n"
				  << "        % (modulus)\n"
				  << "        s (synthetic division: quotient, remainder)\n"
				  << "	d (compute derivative)\n"
				  << "	~ (compare)\n"
                  << "        = (evaluate)\n"
                  << "        M EXP:::[Modular Arithmetic Sub-Menu]\n"
				  << "\t\t\t x (exit)\n"
				  << "\n-----------------------------------------------------------------------\n"
				  << "Operation to perform: ";
		std::cin >> op;
        return op;
}



template<class K>
void getValues(Polynomial<K>& P1, Polynomial<K>& P2)
{
	size_t degree;
	K coefficient;

	cout << "\nEnter degree of first polynomial A(x): ";
	cin >> degree;
	// Polynomial<K> A;
	for (size_t i = degree; i > 0; i--)  {
		cout << "\nEnter coefficient for x^" << i << " : ";
		cin >> coefficient;
		P1.setMember(i, coefficient);
	}
	cout << "\nEnter constant (coefficient for x^0) : ";
	cin >> coefficient;
	P1.setMember((size_t)(0), coefficient);

	cout << "\n\nEnter degree of second polynomial B(x): ";
	cin >> degree;
	// Polynomial<K> B;
	for (long i = degree; i > 0; i--)  {
		cout << "\nEnter coefficient for x^" << i << " : ";
		cin >> coefficient;
		P2.setMember(i, coefficient);
	}
	cout << "\nEnter constant (coefficient for x^0) : ";
	cin >> coefficient;
	P2.setMember((long)(0), coefficient);
}


template<class K>
void getValue(Polynomial<K>& P)   {
	long degree;
	K coefficient;

	std::cout << "\nEnter degree of polynomial P(x): ";
	std::cin >> degree;

	for (long i = degree; i > 0; i--)  {
		cout << "\nEnter coefficient for x^" << i << " : ";
		cin >> coefficient;
		P.setMember(i, coefficient);
	}
	cout << "\nEnter constant (coefficient for x^0) : ";
		cin >> coefficient;
		P.setMember((long)(0), coefficient);

}


template<class K>
void run(K dummy, char status, char op) {
  Polynomial<K> A, B, result, q, r;	
  Polynomial<Rational> poly_Q_;
  Rational c_Q_;
  Polynomial<Complex> poly_C_;
  Complex c_C_;
  Polynomial<double> poly_R_;
  double c_R_;
  char Field = status;
  char operationType = op;

	  
 if ( (operationType == 's') || (operationType == 'S') 
	        || (operationType == 'd') || (operationType == 'D')
			|| (operationType == '=') )  
  {

    if ( (operationType == 's') || (operationType == 'S') )   {  
   
      if ((Field == 'Q') || (Field == 'q'))   {   
	    getValue(poly_Q_);

      std::cout << "Over Rational Field, c = a/b\n";
      std::cout << "\nInclude division symbol between integers, as in 1/2\n";
      std::cout << "\nEnter c for divisor (x - c): c = :";
      std::cin >> c_Q_;
	   
	  Quotient_Remainder<Rational> QR = QR_sDiv(poly_Q_, c_Q_);

      // Display polynomials:
      std::cout << "\nP(x) = " << poly_Q_;
      std::cout << "\n---------------------------------";
      std::cout << "\nD(x) = (x - [" << c_Q_ << "])";

      std::cout << "\n\nP(x)/D(x)\n_________________________________\n";

      std::cout << "\nQuotient = " << QR.first << "\n";
      std::cout << "\nRemainder = " << QR.second << "\n\n";
    }

	if ((Field == 'R') || (Field == 'r'))   {  
	  getValue(poly_R_);

      std::cout << "Over Real Field, c = [enter any decimal number]\n";
      std::cout << "\nEnter c for divisor (x - c): c = ";
      std::cin >> c_R_;
	   
	  Quotient_Remainder<double> QR = QR_sDiv(poly_R_, c_R_);

      // Display polynomials:
      std::cout << "\nP(x) = " << poly_R_;
      std::cout << "\n---------------------------------";
      std::cout << "\nD(x) = (x - [" << c_R_ << "])";

      std::cout << "\n\nP(x)/D(x)\n_________________________________\n";

      std::cout << "\nQuotient = " << QR.first << "\n";
      std::cout << "\nRemainder = " << QR.second << "\n\n";
    }

    if ((Field == 'C') || (Field == 'c'))   
    { 
	  getValue(poly_C_);
    
      std::cout << "Over Complex Field, c = a + bi\n";
      std::cout << "\nEnter c for divisor (x - c): ";
      std::cin >> c_C_;

      Quotient_Remainder<Complex> 
	  QR = QR_sDiv(poly_C_, c_C_);
  
      // Display polynomials:
      std::cout << "\nP(x) = " << poly_C_;
	  std::cout << "\n---------------------------------";
      std::cout << "\nD(x) = (x - [" << c_C_ << "])";
      std::cout << "\n\nP(x)/D(x)\n_________________________________\n";

      std::cout << "\nQuotient = " << QR.first;
      std::cout << "\n\nRemainder = " << QR.second << "\n\n"; 
     }
	 
    }   // END  (s ----> QR_sDiv) 

	if ((operationType == 'd') || (operationType == 'D'))
	{

     if ((Field == 'Q') || (Field == 'q'))   {  
	      getValue(poly_Q_);

          Polynomial<Rational> dP_Q_ = poly_Q_.derive();

          // Display polynomials:
         std::cout << "\nP(x) = " << poly_Q_;
         std::cout << "\n\nP'(x) = " << dP_Q_ << "\n\n";
      }

	 if ((Field == 'R') || (Field == 'r'))   {  
	     getValue(poly_R_);

         Polynomial<double> dP_R_ = poly_R_.derive();

         // Display polynomials:
         std::cout << "\nP(x) = " << poly_R_;
         std::cout << "\n\nP'(x) = " << dP_R_ << "\n\n";
      }

     if ((Field == 'C') || (Field == 'c'))   
     { 
	   getValue(poly_C_);

       Polynomial<Complex> dP_C_ = poly_C_.derive();

       // Display polynomials:
       std::cout << "\nP(x) = " << poly_C_;
       std::cout << "\n\nP'(x) = " << dP_C_ << "\n\n";
      }

   }  // end Derivative

   if (operationType == '=')  {
      
	  if ((Field == 'Q') || (Field == 'q'))   {  
	   getValue(poly_Q_);

	   std::cout << "\nTo evaluate P(x) = " << poly_Q_;
       std::cout << "\nover Rational Field, x = a/b\n";
       std::cout << "\nInclude division symbol between integers, as in 1/2\n";
       std::cout << "\nEnter value of x to evaluate P(x): x = ";
       std::cin >> c_Q_;

       // Display polynomials:
       std::cout << "\nP(" << c_Q_ << ") = " 
	             << poly_Q_.at(c_Q_) << "\n\n";

      }

	  if ((Field == 'R') || (Field == 'r'))   {  
	   getValue(poly_R_);

	   std::cout << "\nTo evaluate P(x) = " << poly_R_;
       std::cout << "\nover Real Field, x = [decimal number]\n";
       std::cout << "\nEnter value of x to evaluate P(x): x = ";
       std::cin >> c_R_;

       // Display polynomials:
       std::cout << "\nP(" << c_R_ << ") = " 
	             << poly_R_.at(c_R_) << "\n\n";
      }

     if ((Field == 'C') || (Field == 'c'))   
     { 
	   getValue(poly_C_);

	   std::cout << "\nTo evaluate P(x) = " << poly_C_;
       std::cout << "\nover Complex Field, x = a + bi\n";
       std::cout << "\nEnter value of x to evaluate P(x): x = ";
       std::cin >> c_C_;

       // Display polynomials:
       std::cout << "\nP(" << c_C_ << ") = " 
	             << poly_C_.at(c_C_) << "\n\n";
      }

	}   // end = (evaluate at)

  }        // end UNARY operations
  else {
 
     // If NOT Exit, perform Binary operation requiring 2 Polynomials:
	 if  ( (operationType != 'x') && (operationType != 'X') )  {
  
                getValues(A, B);
	 }
    switch (operationType)   {
	
	  case '+':
		result = A + B; 
		std::cout << A <<  " + " << B << " = " << result << std::endl;
		break;

	  case '-':
	  	result = A - B; 
		std::cout << A << " - " << B << " = " << result << std::endl;
		break;

	  case '*':
		result = A*B; 
		std::cout << A << " * " << B << " = " << result << std::endl;
		break;

	  case '/':
		//result = A/B; // ----> q: quotient
		if (A.divide(B, q, r))  {
           std::cout << '\n' << A << " / " << B << " = \n\t" << q;
           //result = A%B;	 // ----> r: remainder	

            if (!r.isNull()) {
			   std::cout << " + [(" << r << ")/(" << B; 
			   std::cout << ")]\n";
		    }
           else std::cout << std::endl;
		}
		else std::cout << "\nYou can't divide by NullPolynomial = 0\n\n";
		break;

	  case '%':
		//result = A%B;
        A.divide(B, q, r);
        std::cout << std::endl;
		std::cout << A << " % " << B << " =\n\t= " << r;

	    break;

	  case '~':
		char compare;
		
		if (A == B) compare = '=';
		else compare = '!';

		std::cout << "Test == :" << A << " == " << B 
		          << " ?      " << (bool)( A == B ) << std::endl;

		std::cout << "Test != :" << A << " != " << B
		          << " ?      " << (bool)(A != B ) << std::endl;
		std::cout << std::endl << "Hence, " << A << " " << compare 
		          << " " << B << std::endl;
		
		std::cout << "Test == :" << A << " == " << B << "?      " 
	              << (bool)(A == B ) << std::endl;
		 
		if (A != B) compare = '!';
		else compare = '=';

		std::cout << "Test != " << A << " != " << B << "?      " 
		          << (bool)(A != B ) << std::endl;

		std::cout << std::endl << "Hence, " << A << " " << compare 
		          << " " << B << std::endl;
		break;

	  case 'x':
	    break;

	  default:
		    std::cout << "Your input was invalid.\n"<< std::endl;
	 }  // end swtich

	std::cout << std::endl;
    
	} // end else (binary operations)

}  // exit Calculator 

 /*
  C++11 added alias declarations, which are a generalization of typedef,
  allowing templates - for example,

 template <size_t N>
 using Vector = Matrix<N, 1>;
 The type Vector<3> is equivalent to Matrix<3, 1>

 CHANGE:
 typedef std::pair< std::vector<double>, double > Quotient_Remainder;

 TO:
  template<class F>
  using Quotient_Remainder = std::pair< std::vector<F>, F>;

 CALL with Quotient_Remainder<F>, where <F> could be <std::complex<double>>

 */
