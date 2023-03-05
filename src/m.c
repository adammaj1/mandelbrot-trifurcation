/*





find parameter c from it's multiplier
for 
fc(z) = z^2+c


 mandelbrot-numerics -- numerical algorithms related to the Mandelbrot set
 Copyright (C) 2015-2017 Claude Heiland-Allen
 License GPL3+ 
 http://www.gnu.org/licenses/gpl.html


 
 Files:
  https://code.mathr.co.uk/mandelbrot-numerics/blob_plain/HEAD:/c/bin/m-interior.c
 https://code.mathr.co.uk/mandelbrot-numerics/blob_plain/HEAD:/c/bin/m-util.h
 
 
 description:  
 https://en.wikibooks.org/wiki/Fractals/mandelbrot-numerics#m-interior
 https://en.wikibooks.org/wiki/Fractals/Mathematics/Derivative#First_derivative_wrt_z
 
 
 
 





c console program

gcc m.c -Wall -Wextra -lm
./a.out > m.txt


*/
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdbool.h>








/*
https://mathr.co.uk/blog/2013-04-01_interior_coordinates_in_the_mandelbrot_set.html

There is an algorithm to find points on the boundary of the Mandelbrot set, 
given a particular hyperbolic component and the desired internal angle. 
It involves Newton's method in two complex variables to solve

F^p(z,c)=z
∂F^p(z,c)/∂z =b

where 
F0(z,c)=z and 
Fq+1(z,c)=Fq(F(z,c)2+c), 
p is the period of the target component, 
and m=e2πiθ with the θ the desired internal angle. 
The resulting c is the coordinates of the point on the boundary. 
It can also be modified to find points in the interior, simply set m=re2πiθ with |r|≤1.

*/



// mandelbrot-numerics/c/include/mandelbrot-numerics.h


enum m_newton { m_failed, m_stepped, m_converged };
typedef enum m_newton m_newton;

// mandelbrot-numerics/c/bin/m-util.h

const double twopi = 6.283185307179586;
// epsilon^2
static const double epsilon2 = 1.9721522630525295e-31;



static inline double cabs2(double complex z) {
  return creal(z) * creal(z) + cimag(z) * cimag(z);
}

static inline bool cisfinite(double complex z) {
  return isfinite(creal(z)) && isfinite(cimag(z));
}


complex double GiveFixed(complex double c){
/* 

Equation defining fixed points : z^2-z+c = 0
	z*2+c = z
	z^2-z+c = 0

coefficients of standard form ax^2+ bx + c  
 	a = 1 , b = -1 , c = c
 
The discriminant d is 

	d = b^2- 4ac 
	d = 1 - 4c

 alfa =  (1-sqrt(d))/2 
*/

	complex double d = 1-4*c;
	complex double z = (1-csqrt(d))/2.0;
	return z;

}


// mandelbrot-numerics/c/lib/m_d_interior.c 
// double precision: m_d_*()  


m_newton m_d_interior_step(double complex *z_out, double complex *c_out, double complex z_guess, double complex c_guess, double complex multiplier, int period) {
  double complex c = c_guess;
  double complex z = z_guess;
  double complex dz = 1;
  double complex dc = 0;
  double complex dzdz = 0;
  double complex dcdz = 0;
  for (int p = 0; p < period; ++p) {
    dcdz = 2 * (z * dcdz + dc * dz);
    dzdz = 2 * (z * dzdz + dz * dz);
    dc = 2 * z * dc + 1;
    dz = 2 * z * dz;
    z = z * z + c;
  }
  double complex det = (dz - 1) * dcdz - dc * dzdz;
  double complex z_new = z_guess - (dcdz * (z - z_guess) - dc * (dz - multiplier)) / det;
  double complex c_new = c_guess - ((dz - 1) * (dz - multiplier) - dzdz * (z - z_guess)) / det;
  if (cisfinite(z_new) && cisfinite(c_new)) {
    *z_out = z_new;
    *c_out = c_new;
    if (cabs2(z_new - z_guess) <= epsilon2 && cabs2(c_new - c_guess) <= epsilon2) {
      return m_converged;
    } else {
      return m_stepped;
    }
  } else {
    *z_out = z_guess;
    *c_out = c_guess;
    return m_failed;
  }
}

/* 
  usage: ./m_d_interior z-guess c-guess interior period maxsteps
  
  input: 
  	z-guess 	=  	the initial guess z_0 for Newton method, use z = 0
  	c-guess 	= 	the initial guess c_0 for Newton method, use nucleus of given hyperbolic component
  	interior	=  multiplier 
	period	= 	Period of given hyperbolic component
	maxsteps	= 
  
  
  
  My c program gives wrong result:
 period = 3 center = -0.1225611668766540+0.7448617666197440*I    multiplier : internal angle = 0.00  internal radius = 1.0000000000000000   c = -1000.0000000000000000+0.0000000000000000*I
original code gives good result:
  ./m-interior double 0 0 -0.1225611668766540 +0.7448617666197440 1 0 3 100
find point c of component with period = 3     multiplier = 1.0000000000000000+0.0000000000000000    located near c=  -0.1225611668766540+0.7448617666197440
-2.5006361141355088e-01 4.3298542796573969e-01 -1.2499999999943506e-01 6.4951905283823752e-01
https://fractalforums.org/programming/11/multiplier-map/5038/msg36763#new
The return policy of m_d_interior should be changed as follows:
  
  
m_newton m_d_interior(double complex *z_out, double complex *c_out, double complex z_guess, double complex c_guess, double complex multiplier, int period, int maxsteps) {
	m_newton result = m_failed;
  	double complex z = z_guess;
  	double complex c = c_guess;
  
  	for (int i = 0; i < maxsteps; ++i) {
    		if (m_stepped != (result = m_d_interior_step(&z, &c, z, c, multiplier, period))) 
    			{ break;  }
  		}
  	// 	
  	*z_out = z;
  	*c_out = c;
  	return result;
}
*/

m_newton m_d_interior(double complex *z_out, double complex *c_out, double complex z_guess, double complex c_guess, double complex multiplier, int period, int maxsteps) {

	m_newton result = m_failed;
  	double complex z = z_guess;
  	double complex c = c_guess;
 
  	for (int i = 0; i < maxsteps; ++i) {
    		if (m_stepped != (result = m_d_interior_step(&z, &c, z, c, multiplier, period)))
    			{ break;  }
  		}
  	// 	
  	*z_out = z;
  	*c_out = c;

	// if the result variable here has value
	// m_stepped, the loop above has not been left
	// by the break statement, hence was successful
	// in all iterations
	// then the function should return
	// m_converged, as the
	// outside routine aproximate_c tests only
	// for m_converged
	
	if (result == m_stepped) return m_converged;
	
  	return result;
}




complex double aproximate_c( const int p, const complex double center, const complex double multiplier){


	complex double c = 0.0;
	complex double z = 0;
	int maxsteps = 100;
	
	m_newton result;
		// m_d_interior(&z, &c, zre + I * zim, cre + I * cim, ir * cexp(I * twopi * it), period, maxsteps);
	result = m_d_interior(&z,  &c, 0.0, center, multiplier, p, maxsteps);
    	if (result != m_converged) 
    		{return -1000;}
	
	
	
	return c;
}




// m_d_interior(&z, &c, zre + I * zim, cre + I * cim, ir * cexp(I * twopi * it), period, maxsteps);

complex double give_c(const int p, const complex double center, const double angle, const double r )
{
	/*
	input:
	Internal Radius = r in [0,1] ; double
  	Internal Angle In Turns = t  or theta in range [0,1) ; double  
  	p = period ; int 
  	
  	output = c = complex point of 2D parameter plane  
  	*/
  	

	complex double m = 0.0; // multiplier
	complex double c = 0.0; // result 
	double t = angle;
	
	t = t*twopi; // from turns to radians
	m = r* cexp(I*t); // point of unit circle = multiplier = ir * cexp(I * twopi * it)
  		
	// map circle to component
	switch (p){
	
		case 1: c = (2.0*m - m*m)/4.0; break;
		case 2: c = (m -4.0)/ 4.0; break;
		default : c = aproximate_c( p, center, m); // for higher periods there is no exact method; use numerical aproximation	 
	}
	
	return c; 
}








// from center to root
complex double GiveInternalRay(const int period, const complex double center, const double InternalAngle, const int iMax)
{
	complex double c;
	// double InternalAngle; in turns = from 0.0 to 1.0
	double InternalRadius ;
	double RadiusStep = 1.0/iMax;
	int i;
	complex double zf;
	
	InternalRadius = 0.0 + RadiusStep; // without center 
	for (i=2; i <= iMax; ++i) { 
		c = give_c(period, center, InternalAngle, InternalRadius);
		
		zf = GiveFixed(c);
		//printf ("%.16f\t %.16f \t %+.16f\n",InternalRadius,  creal(zf), cimag(zf)); 
		printf ("%.16f\t %.16f \t %+.16f\t %.16f \t %+.16f\n", InternalRadius, creal(c), cimag(c), creal(zf), cimag(zf));
				
		InternalRadius = i * RadiusStep;

	}
	
	
	return c;
}


// from root to center
complex double GiveInternalRayInv(const int period, const complex double center, const double InternalAngle, const int iMax)
{
	complex double c;
	// double InternalAngle; in turns = from 0.0 to 1.0
	double InternalRadius;
	double RadiusStep = 1.0/iMax;
	int i;
	complex double zf;

	InternalRadius= 1.0 - RadiusStep; //
	
	for (i=iMax-2; i >= 0; --i) { 
		c = give_c(period, center, InternalAngle, InternalRadius);
		
		zf = GiveFixed(c);
		//printf ("%.16f\t %.16f \t %+.16f\n", InternalRadius, creal(zf), cimag(zf)); 
		printf ("%.16f\t %.16f \t %+.16f\t %.16f \t %+.16f\n", InternalRadius, creal(c), cimag(c), creal(zf), cimag(zf));
		 		
		InternalRadius = i * RadiusStep;

	}
	
	
	return c;
}




m_newton m_d_attractor_step(double _Complex *z, double _Complex z_guess, double _Complex c, int period) {
  double _Complex zz = z_guess;
  double _Complex dzz = 1;
  for (int i = 0; i < period; ++i) {
    dzz = 2 * zz * dzz;
    zz = zz * zz + c;
  }
  if (cabs2(zz - z_guess) <= epsilon2) {
    *z = z_guess;
    return m_converged;
  }
  double _Complex z_new = z_guess - (zz - z_guess) / (dzz - 1);
  double _Complex d = z_new - z_guess;
  if (cabs2(d) <= epsilon2) {
    *z = z_new;
    return m_converged;
  }
  if (cisfinite(d)) {
    *z = z_new;
    return m_stepped;
  } else {
    *z = z_guess;
    return m_failed;
  }
}

m_newton m_d_attractor(double _Complex *z_out, double _Complex z_guess, double _Complex c, int period, int maxsteps) {
  m_newton result = m_failed;
  double _Complex z = z_guess;
  for (int i = 0; i < maxsteps; ++i) {
    if (m_stepped != (result = m_d_attractor_step(&z, z, c, period))) {
      break;
    }
  }
  *z_out = z;
  return result;
}





// *****************************************************

int main (){

	const int child_period = 3;
	const complex double child_center = -0.122561166876654  +0.744861766619744*I; //   period = 3; // = nucleus = center of hyperbolic component of the Mandelbrot set with period p ; there are more then 1 components with periods > 2
	
	
	const int parent_period = 1;
	const complex double parent_center = 0.0; // = nucleus = center of hyperbolic component of the Mandelbrot set with period p ; there are more then 1 components with periods > 2
	
	complex double c;
	complex double zf;
	//complex double z31;
	
	int number_of_points = 10;
	
	
	
	// give_c(const int p, const complex double center, const double angle, const double r )

	printf ("#  radius cx cy  zxf  zyf \n"); // header for gnuplot data file
	// printf ("%.16f\t %.16f \t %+.16f \t %.16f \t %+.16f \t %.16f \t %+.16f\n", 0.0, creal(c), cimag(c), creal(zf), cimag(zf), creal(z31), cimag(z31));
	
	c = parent_center;
	zf = GiveFixed(parent_center);
	//m_d_attractor(&z31, 0.0, c, child_period, 100);
	//printf ("%.16f\t %.16f \t %+.16f \t %.16f \t %+.16f \t %.16f \t %+.16f\n", 0.0, creal(c), cimag(c), creal(zf), cimag(zf), creal(z31), cimag(z31)); 
	printf ("%.16f\t %.16f \t %+.16f\t %.16f \t %+.16f\n", 0.0, creal(c), cimag(c), creal(zf), cimag(zf)); 
	
	// int ray
	GiveInternalRay(parent_period, parent_center, 1.0/child_period, 3*number_of_points);
	
	
	// root is a common point of 2 componnets 
	c = give_c(parent_period, parent_center, 1.0/child_period, 1.0);
	zf = GiveFixed(c);
	printf ("%.16f\t %.16f \t %+.16f\t %.16f \t %+.16f\n", 1.0, creal(c), cimag(c), creal(zf), cimag(zf)); 
	
	
	
	// int ray 
	GiveInternalRayInv(child_period, child_center, 0.0, number_of_points);
	
	
	
	c = child_center; // child center 
	zf = GiveFixed(c);
	//m_d_attractor(&z31, 0.0, c, child_period, 100);
	printf ("%.16f\t %.16f \t %+.16f\t %.16f \t %+.16f\n", 0.0, creal(c), cimag(c), creal(zf), cimag(zf)); 
	//printf ("%.16f\t %.16f \t %+.16f \t %.16f \t %+.16f \t %.16f \t %+.16f\n", 0.0, creal(c), cimag(c), creal(zf), cimag(zf), creal(z31), cimag(z31)); 
	
	return 0;
}

