# include <stdlib.h>		// (to include exit(1)) 
#include <stdio.h>
#include <math.h>

double th_Gamma(double x); 
double th_Hyper1F1GG(double a, double b, double x, double eps); 
double th_Uax9Strip(double aa, double x, double wStart, double Dw);
double th_Hermite(int n, double x);
double th_Power(double x, double p);
double for_large_negative(double a, double x);

#define asinh(x) ( log( (x)+sqrt((x)*(x)+1) ) )

/***************************************************************************
// return log(U(a,x))
// According to Abramowitz/Stegun this function provides a good approximation
// to U(a,x) when x^2+4a is large. 
// Further this function requires a > 0
****************************************************************************/

double for_large_negative(double a, double x)
{
	double d3, d6, v_ax, X, theta, d9, d12;

	const double MY_PI = acos(-1);

	// Abramowitz and Stegun page 690 19.10.13

	d3 = (x*x*x/48.0+a*x/2.0)/a;

	d6 = 3.0*x*x/4.0 - 2.0*a;

	d9 = (-7*x*x*x*x*x*x*x*x*x/5760 - 7*a*x*x*x*x*x*x*x/320 - 49*a*a*x*x*x*x*x/320 
		+ 31*a*a*a*x*x*x/12 - 19*a*a*a*a*x)/(a*a*a);

	d12 = 153*x*x*x*x/8 - 186*a*x*x + 80*a*a;

	// from A/S 19.10.1

	X = sqrt(x*x+4.0*a);

	theta = fabs(x)*X/4.0 + a*asinh(fabs(x)/(2*sqrt(a)));  // notice the fabs()'s I need this to get the correct
																													// theta when x is negative (not completely clear from
	                                                       // A/S I think)

	// A/S 19.10.4

	v_ax = -0.5*log(X) -d3/(X*X*X) + d6/(X*X*X*X*X*X) - d9/(X*X*X*X*X*X*X*X*X)
		+ d12/(X*X*X*X*X*X*X*X*X*X*X*X);


	// from 19.10.3

	return 0.25*log(2*MY_PI) - 0.5*log(th_Gamma(0.5+a)) + theta + v_ax;
}

#undef asinh


/********************************************************************
// Parabolic Cylinder function U(a,x). Fractional accuracy goal is eps for series.
// ParCylU is copied from Thompson but modified here to suit the purpose of the 
// extinction group program.
// E.g. a function which can treat large negative x is encluded "for_large_negative".
// and the function now returns log(exp(x^2/4)U(a,x))
// The parabolic cylinder function in Whittaker notation is U(a,x) = D_{-a-1/2}(x). 
//
// Credits -- Atlas for Computing Mathematical Functions: An Illustrated Guide 
//            for Practitioners in C and Mathematica by William J. Thompson
*********************************************************************/

double ParCylU(double a, double x, double eps)
{
	double aPh,aD2,Dw,integral,epsI,wStart,xS,xSD2,xSD4,
		wMax,DwStart,Uax,epsK;
  //double rem, Ma;
	const double RRt2 = 0.707106781187,RtPi = 1.77245385091;
	int kMax,k;
  //int flMa;
	
	aPh = a + 0.5;
	aD2 = 0.5*a;
	xS = x*x;
	xSD2 = 0.5*xS;
	xSD4 = 0.0; //0.5*xSD2;  
	if ( ( x > 3.0 ) && ( a > 1.0 ) )
	{  /* Use 9-strip integrals */
		wMax = 0.5*x*(sqrt(x*x+4.0*a-2.0) - x);
		wMax = ( wMax > 2.5 ) ? 8.0*wMax : 25.0;
		kMax = 90;  /* is number of blocks */
		DwStart = wMax/( (double) kMax );
		Dw = DwStart/27.0;  /* is step size in a block */
		/* Analytical integral near origin */
		epsI = 0.05;
		integral = pow(epsI,aPh)*(1.0/aPh-epsI/(aPh+1.0)+
			0.5*epsI*epsI*((1.0-1.0/xS)/(aPh+2.0)+
			epsI*(1.0/xS-1.0/3.0)/(aPh+3.0)));
		wStart =  epsI;
		for ( k = 1; k <= kMax; k++ )
		{
			integral += th_Uax9Strip(a,x,wStart,Dw);
			wStart +=  DwStart; 
		}
		Uax = integral/(exp(xSD4+aPh*log(x))*th_Gamma(aPh));
		return  Uax;
	}
	if ( a <= 1.0 )
	{
		printf("The custom made parametric cylinder function is\n");
		printf("is not garantied to work form a <= 1.0\n\n");
		exit(1);
	}

	if ( x < -10 )
	{
		return for_large_negative(a, x) + x*x/4.0;
	}

//	if ( a < 0.0 )
//	{   /* th_Hermite polynomials */
//		Ma = -a;
//		flMa = floor(Ma);
//		rem = Ma - (double) flMa ;
//		if ( fabs(rem - 0.5) < 1.0e-6 ) /* th_Hermite polynomial */
//		{ return  th_Hermite(flMa,RRt2*x)/
//		(pow(2.0,0.5*((double)flMa))*exp(xSD4)); }
//	}
	/* Hypergeometric series */
	epsK = 0.001*eps;
	Uax = RtPi*exp(-xSD4)*
		(pow(2.0,-0.25-aD2)*th_Hyper1F1GG(aD2+0.25,0.5,xSD2,epsK)/
		th_Gamma(aD2+0.75) -
		x*pow(2.0,0.25-aD2)*th_Hyper1F1GG(aD2+0.75,1.5,xSD2,epsK)/
		th_Gamma(aD2+0.25));
	return  log(Uax);
}


double th_Uax9Strip(double aa, double x,
								 double wStart, double Dw)
								 /* Nine-strip integral for parabolic cylinder U
								 with  aa > 0   */
{
	double aMh,R2xS,w,f[28],integral;
	double th_Power(double x, double p);
	const double a[5] = {0.2869754464285714,
		1.581127232142857, 0.1084821428571429,
		1.943035714285714, 0.5803794642857144 };
	int k;
	
	aMh = aa - 0.5;
	R2xS = 1.0/(2.0*x*x);
	w = wStart;
	for ( k = 0; k <= 27; k++ )
	{
		f[k] = exp(-w*(1.0+w*R2xS))*th_Power(w,aMh);
		w += Dw;
	}
	integral = Dw*(a[0]*(f[0]+f[27]+2.0*(f[9]+f[18]))+
		a[1]*(f[1]+f[8]+f[10]+f[17]+f[19]+f[26])+
		a[2]*(f[2]+f[7]+f[11]+f[16]+f[20]+f[25])+
		a[3]*(f[3]+f[6]+f[12]+f[15]+f[21]+f[24])+
		a[4]*(f[4]+f[5]+f[13]+f[14]+f[22]+f[23]));
	return  integral;
}


double th_Hermite(int n, double x)
/* th_Hermite Polynomial of order n */
{
	double HiM2,HiM1,Hi;
	int i;
	
	/* Explicit for  n < 4 */
	if ( n < 1 )	return 	1.0;  /* H0 */
	if ( n < 2 ) 	return 	2.0*x;
	if ( n < 3 )	return 	-2.0+4.0*x*x; /* H2 */
	if ( n < 4 ) 	return 	4.0*x*(-3.0+2.0*x*x);
	
	/* Use iteration for  n >= 4  */
	/* H2  then  H3  */
	HiM2 = -2.0+4.0*x*x;
	HiM1 = 4.0*x*(-3.0+2.0*x*x);
	for ( i = 4; i <= n; i++ )
	{
		Hi = 2.0*(x*HiM1 - ((double)i-1.0)*HiM2);
		HiM2 = HiM1;
		HiM1 = Hi;
	}
	return	Hi;
}


double th_Hyper1F1GG(double a, double b,
									double x, double eps)
									/* Hypergeometric function 1F1
									to fractional accuracy  eps.
									Uses  th_Gamma  rather than  Logth_Gamma   */
{
	double  sum,ratio,aks,bks,kfact,tk,xpow,
		eps10,term,fk,bMaPkM1,OneMaPkM1,termold,termnew,
		aPkM1,aMbPk,Mx,th_Gamma(double x);
	int k;
	
	if ( a == b )  return  exp(x);
	
	if ( b <= 0.0 )
	{
		printf
			("\n\n!!In th_Hyper1F1GG  b<=0; 0 returned");
		return  0.0;
	}
	
	eps10 = 0.1*eps;
	if ( fabs(x) <= 50.0 )
	{
		sum = 1.0;
		ratio = 10.0;
		aks = a;
		bks = b;
		kfact = 1.0;
		tk = 1.0;
		xpow = 1.0;
		while ( ratio > eps10 )
		{
			tk *= aks/(bks*kfact);
			xpow *= x;
			term = tk*xpow;
			sum += term;
			aks += 1.0;
			bks += 1.0;
			kfact += 1.0; 
			ratio = fabs(term/sum);
		}
		return  sum;
	}
	/* Use asymptotic series */
	if ( x > 50.0 )
	{
		if ( a <= 0.0 )
		{
			printf
				("\n\n!!In th_Hyper1F1GG asymptotic a<=0; 0 returned");
			return  0.0;
		}
		k = 0;
		fk = 1.0;
		bMaPkM1 = b - a;
		OneMaPkM1 = 1.0 - a;
		term = 1.0;
		termold = 100.0;
		sum = 1.0;
		ratio = 10.0;
		while ( ratio > eps10 )
		{
			k++;
			term *=  bMaPkM1*OneMaPkM1/(fk*x);
			sum += term;
			bMaPkM1 += 1.0;
			OneMaPkM1 += 1.0;
			fk += 1.0;
			ratio = fabs(term/sum);
			termnew = fabs(term);
			if ( termnew > termold )
			{
				printf
					("\n\n!!In th_Hyper1F1GG diverging after %i terms",k);
				ratio = 0.0;
			}
			else
			{ termold = termnew; }
		}
		return  (th_Gamma(b)*exp(x)/th_Gamma(a))*pow(x,a-b)*sum;
	}
	/*  x < -50  asymptotic series   */
	if ( b <= a )
	{
		printf
			("\n\n!!In asymptotic th_Hyper1F1GG  b<=a; 0 returned");
		return  0.0;
	}
	Mx = -x;
	k = 0;
	fk = 1.0;
	aMbPk = a - b + 1.0;
	aPkM1 = a;
	term = 1.0;
	termold = 100.0;
	sum = 1.0;
	ratio = 10.0;
	while ( ratio > eps10 )
	{
		k++;
		term *= aPkM1*aMbPk/(fk*Mx);
		sum += term;
		aPkM1 += 1.0;
		aMbPk += 1.0;
		fk += 1.0;
		ratio = fabs(term/sum);
		termnew = fabs(term);
		if ( termnew > termold )
		{
			printf
				("\n\n!!In th_Hyper1F1GG diverging after %i terms",k);
			ratio = 0.0;
		}
		else
		{ termold = termnew; }
	}
	return  th_Gamma(b)*sum/(th_Gamma(b-a)*pow(Mx,a));
}


double th_Power(double x, double p)
/* x to th_Power p, including 0 to 0 = 1 */
{
	if ( x == 0.0 && p == 0.0 )  return 1.0;
	return  pow(x,p);
}


double th_Gamma(double x)
/* th_Gamma function; uses asymptotic for log th_Gamma */
{
	double Recxs,sum,term,x10;
	const double EulerG = 0.577215664902,
		hlfLn2Pi = 0.918938533205;
	const double b[] =
	{0.0, 0.0833333333333, -0.00277777777778, 
  0.000793650793651, -0.000595238095238, 
  0.000841750841751};
	const double Pi = 3.14159265359;
	int xInt,k,x11;
	
	if ( x > 10.0 ) /* use asymptotic series */
	{
		Recxs = 1.0/(x*x);
		term = x;
		sum = (x-0.5)*log(x)-x+hlfLn2Pi;
		for ( k = 1; k <= 5; k++ )
		{
			term *= Recxs;
			sum += b[k]*term;
		}
		return  exp(sum);
  }
	if ( x > 0 ) /* recurrence to x > 10 */
	{
		x -= 1.0;
		x11 = (int) 11.0 - x;   // not sure whether this double-to-int conversion here is 100% intensional
		x10 = x + (double) x11;
		xInt = (int) x;  // again, not sure if double-to-int conversion 100% intensional
                     // but hopefully ok. As far as I remember I copied this code from
                     // William J. Thompson's book credited above.
		sum = 0.0;
		for ( k = 1; k <= x11-1; k++ )
		{
			sum -= log(x10 - (double) k );
		}
		return exp(sum)*th_Gamma(x10);
	} 
	if ( x == 0.0 )
	{
		printf("\n!!x=0 in th_Gamma; zero returned");
		return 0.0;
	}
	/* Use reflection formula */;
	return  Pi/(sin(Pi*x)*th_Gamma(1.0-x)); 
}