
1. Four multiple roots
	(fl01 is taken from M.R. Farmer and G. Loizou, algorithm for the total,
		or partial factorization of a polynomial, Math. Proc. Camb. 
		Phil. Soc. Vol 82, pp427-437, 
	 fl02-fl07 are used by Z. Zeng)

	fl01: 	k=1,  	(x-1)^(4k) * (x-2)^(3k) * (x-3)^(2k) * (x-4)^k
	fl02: 	k=2,  	(x-1)^(4k) * (x-2)^(3k) * (x-3)^(2k) * (x-4)^k
	fl03: 	k=3,  	(x-1)^(4k) * (x-2)^(3k) * (x-3)^(2k) * (x-4)^k
	fl04: 	k=4,  	(x-1)^(4k) * (x-2)^(3k) * (x-3)^(2k) * (x-4)^k
	fl05: 	k=5,  	(x-1)^(4k) * (x-2)^(3k) * (x-3)^(2k) * (x-4)^k
	fl06:  	k=6,  	(x-1)^(4k) * (x-2)^(3k) * (x-3)^(2k) * (x-4)^k
	fl07:  	k=7, 	(x-1)^(4k) * (x-2)^(3k) * (x-3)^(2k) * (x-4)^k


2.  Twin multiple roots 0.39 and 0.4, plus the third multiple root -0.2

	twin01		k = 4, 	(x-0.39)^k * (x-0.4)^k * (x+0.2)^k
	twin02		k = 8, 	(x-0.39)^k * (x-0.4)^k * (x+0.2)^k
	twin03		k = 12, 	(x-0.39)^k * (x-0.4)^k * (x+0.2)^k
	twin04		k = 16, 	(x-0.39)^k * (x-0.4)^k * (x+0.2)^k

3.  Frank Uhlig's testing polynomials
	F. Uhlig, General polynomial roots and their multiplicities in O(n)
	memory and O(n^2) time, Linear and Multilinear Algebra, Vol. 46,
	pp 327-359, 1999

	uhlig01		a = 0.01, 	(x-a^4)(x-a)^4
	uhlig02		a = 0.001, 	(x-a^4)(x-a)^4
	uhlig03		(x-3/11)^12 * (x-11/3)^2 * (x-2i/7)^4 * 
				(x-2.5-i/4)^2 * (x-1/4)
	uhlig04		(x-3/11)^12 * (x-11/3)^2 * (x-2i/7)^4 * 
				(x-2.5-i/4)^2 * (x-1/8)
	uhlig05		(x+1)^6 * (x+2)^2;

4.  Jenkins-Traub testing polynomials
	M.A. Jenkins and J.F. Traub, Principles for testing polynomial zerofinding
	programs, ACM Trans Math Software, Vol. 1, pp26-34, 1975

	jt01a 	a = 10^10, 	(x-A)(x+A)(x-1)
	jt01b 	a = 10^(-10), 	(x-A)(x+A)(x-1)
	jt02	(x-1)(x-2)...(x-17)
	jt03	(x-0.1)(x-0.001)...(x-0.00000001)
	jt04	(x-.1)^3 * (x-.5)(x-.6)(x-.7)
	jt05	(x-.1)^4 * (x-.2)^3 * (x-.3)^2 * (x-.4)
	jt06	(x-.1)(x-1.001)(x-.998)(x-1.00002)(x-.99999)
	jt07a	a = 10^(-10)
		(x-.001)(x-.01)(x-.1)(x-.1+a*i)(x-.1-a*i)(x-1)(x-10)
	jt07b	a = 10^(-9)
		(x-.001)(x-.01)(x-.1)(x-.1+a*i)(x-.1-a*i)(x-1)(x-10)
	jt07c	a = 10^(-8)
		(x-.001)(x-.01)(x-.1)(x-.1+a*i)(x-.1-a*i)(x-1)(x-10)
	jt07d	a = 10^(-7)
		(x-.001)(x-.01)(x-.1)(x-.1+a*i)(x-.1-a*i)(x-1)(x-10)
	jt08	(x+1)^5
	jt09	(x^10-10^(-20)) * (x^10+10^20)
	jt10a	a = 10^3,	(x-a)(x-a)(x-1/a)
	jt10b	a = 10^6,	(x-a)(x-a)(x-1/a)
	jt10c	a = 10^9,	(x-a)(x-a)(x-1/a)
	jt11a	m = 15, 	roots:  	exp(i*k*pi/(2*m)), k=1-m:m-1;
				0.9*exp(i*k*pi/(2*m)), k=m:3*m
	jt11a	m = 20, 	roots:  	exp(i*k*pi/(2*m)), k=1-m:m-1;
				0.9*exp(i*k*pi/(2*m)), k=m:3*m
	jt11a	m = 25, 	roots:  	exp(i*k*pi/(2*m)), k=1-m:m-1;
				0.9*exp(i*k*pi/(2*m)), k=m:3*m

5.  Goedecker's testing polynomials
	S. Goedecker, Remark on algorithms to find roots of polynomials,
	SIAM J. Sci. Comput, Vol. 15, pp 1059-1063, 1994

	fib(n)	Fibocacci polynomial of degree n	p = ]-1,1,1,...,1]
	fib05	Fibocacci polynomial, n=5,		p = [-1,1,1,...,1]
	fib10	Fibocacci polynomial, n=10,		p = [-1,1,1,...,1]
	fib15	Fibocacci polynomial, n=15,		p = [-1,1,1,...,1]
	fib20	Fibocacci polynomial, n=20,		p = [-1,1,1,...,1]
	fib30	Fibocacci polynomial, n=30,		p = [-1,1,1,...,1]
	fib50	Fibocacci polynomial, n=10,		p = [-1,1,1,...,1]
	fib100	Fibocacci polynomial, n=100,	p = [-1,1,1,...,1]
	fib150	Fibocacci polynomial, n=150,	p = [-1,1,1,...,1]

	fibsq04	squared Fibocacci polynomial, n=4;
	fibsq08	squared Fibocacci polynomial, n=8;
	fibsq16	squared Fibocacci polynomial, n=16;
	fibsq24	squared Fibocacci polynomial, n=24;
	fibsq32	squared Fibocacci polynomial, n=32;
	fibsq48	squared Fibocacci polynomial, n=48;

	lgd(n)	Legendre polynomial of any degree n
	lgd05	Legendre polynomial, n=5
	lgd10	Legendre polynomial, n=10
	lgd15	Legendre polynomial, n=15
	lgd20	Legendre polynomial, n=20
	lgd24	Legendre polynomial, n=24
	lgd50	Legendre polynomial, n=50
	lgd100	Legendre polynomial, n=100

6. Miyakoda, Iterative methods for multiple zeros of a polynomial by clustering, 
	J. of Computational and Applied Mathematics, Vol. 28, pp315-326, 1989

	miyak00		(x-1.1-1.i*i)^4 * (x-3.2-2.3*i)^2 * (x-2.1-1.5*i)
	miyak02		square of miyak00
	miyak04		square of miyak02
	miyak08		square of miyak04

7. A cluster of three multiple roots (Z. Zeng)

	triple(m,n,k)	(x-0.9)^m * (x-1)^n * (x-1.1)^k
	
	triple01		(m,n,k) = (5,5,5);
	triple02		(m,n,k) = (10,10,10);
	triple03		(m,n,k) = (18,10,16)
	triple04		(m,n,k) = (20,15,10);

8.  Inexact polynomial (Z. Zeng)

	(x-10/11)^5 * (x-20/11)^3 * (x-30/11)^2

	inex01	above polynomial with 10 digits accuracy in each coefficient
	inex02	above polynomial with 9 digits accuracy in each coefficient
	inex03	above polynomial with 8 digits accuracy in each coefficient
	inex04	above polynomial with 7 digits accuracy in each coefficient

	(Multroot requires adjustment of tol)

9. Nearby multiple roots (Z. Zeng)

	(x-1+e)^20 * (x-1)^20 * (x+0.5)^5

	near01	e = 0.1
	near02  	e = 0.01
	near03	e = 0.001

10. A large problem (Z. Zeng)

	roots: 	 1.00000000000247 + 0.30000000000157i
 		-1.00000000000247 - 0.30000000000157i
 		-0.89999999999333 + 0.39999999999963i
		 -0.89999999999333 - 0.39999999999963i
 		-0.70000000000170 + 0.69999999999681i
 		-0.70000000000170 - 0.69999999999681i
		 -0.40000000000029 + 0.90000000000610i
		 -0.40000000000029 - 0.90000000000610i
		 -0.00000000000046 + 1.10000000000425i
		 -0.00000000000046 - 1.10000000000425i
 		 1.19999999999897                    
  		1.00000000000000                    
  		0.89999999999807 + 0.39999999999939i
  		0.89999999999807 - 0.39999999999939i
  		0.59999999999654 + 0.60000000000072i
  		0.59999999999654 - 0.60000000000072i
  		0.39999999999631 + 0.90000000000011i
  		0.39999999999631 - 0.90000000000011i
  		0.00000000000785 + 0.79999999997916i
  		0.00000000000785 - 0.79999999997916i

	large01		polynomial with above roots
	large02		square of large01
	large03		square of large02
	large04		square of large03
	large05		square of large04
	

11. M. Petkovic, Iterative Methods for Simultaneous Inclusion of Polynomial Zeros,
	Lecture Notes in Mathematics, Springer-Verlag, 1989
	Chapter 4

	petk01	(x+1)^2 * (x-3)^3 * (x-i)^4 (x^2-2x+5)^2
	petk02	(x-1)^2 * (x+i)^3 * (x-5i)^2 * (x+5i)^2
	petk03 	70(x^2-2x+3)^2 * (x-1-99i/70) * (x+1);
	petk04	(x-3) * (x+1)^3 * (x-2i)^3 * (x^2+4x+5)^2 * (x^2-4x+5)^2
	petk05	(x-3)^2 * (x^2+2x+5)^2 * (x+1)^3
	petk06	(x-3)^2 * (x^2+2x+5)^2 * (x+1)^4 * (x+i)^2
	petk07	(x-1)^3 * (x^2-4x+5) * (x^2+25)^2
	henrici	(x+4.1)(x+3.8)(x+2.05)(x+1.85)(x-1.95)(x-2.15)(x-3.9)(x-4.05)

	
12. Farmer-Loizou

	farloi01		[(x^2+x+2)(x^2+x+3)]^4
	

13. L. Brugnano, D. Trigiante, Polynomial roots: the ultimate answer? Linear Algebra
	and Its Applications, Vol. 225, pp207-219, 1995

	bt01	(x-1)^6 (x+1)^2 (x+i)^3 (x-i)^3 (x-2)
	bt02	(x-1)^10 (x-2)^2 (x+i) (x-i)
	bt03	(x^2+1)^5 (x-0.5i)^4 (x+0.5i)^4 (x-0.75i)(x+0.75i)
	bt04	(x-1)^3 (x+1)^4 (x-.5-i)^3 (x-.5+i)^3 (x-.5-.5i)^2 (x-.5+.5i)^2

14. A. I. Iliev 
	iliev00 is original, others are extention.

	iliev00	(x+2)^2 (x-1) (x-4)^3
	iliev01 	(x-1)^2 (x-2)^4 (x-3)^6
	iliev02	(x-1)^4 (x-2)^8 (x-3)^12
	iliev03 	(x-1)^8 (x-2)^16 (x-3)^24

15. M. Igarashi and T. Ypma, Relationsips between order and efficiency of a class of
	methods for multiple zeros of polynomials, J. of Computational and Applied
	Mathematics, Vol. 60, pp101-113, 1995

	igyp00	(x-2.35)(x-2.37)(x-2.39)
	igyp01	(x-2.35)^3 (x-2.56)
	igyp02a 	m=8,	(x-10-10i)^m (x+1)^(10-m);
	igyp0ab	m=7,	(x-10-10i)^m (x+1)^(10-m);
	igyp0ab	m=3,	(x-10-10i)^m (x+1)^(10-m);


16. K.C. Toh and L. N. Trefethen, Pseudozeros of polynomials and pseudospectra
	of companion matrices, Numer. Math., Vol. 68, pp403-425, 1994
      cited and modified by
      H. Zhang, Numerical condition of polynomials in different forms, Electronic
	Transactions on Numerical Analysis, Vol. 12, pp. 66-87, 2001

	toh01	monic polynomial with simple roots
		2(k+0.5)/19 + i*sin 2*pi*(k+0.5)/19,   k=-10,-9,...,9
	toh02	sum of  (10x)^k / k!, k=0,1,...,20
	toh03	monic polynomial with roots 10/11-2^(-k), k=1,2,...,20
	toh04	(x-10/11)^20
	toh05	monic polynomial with roots 2^(-k), k=0,1,...,19
	toh06a	1 + x + x^2 + ... + x^20
	toh06b	(1+x+x^2+...+x^10)^2
	toh06c	(1+x+...+x^5)^4
