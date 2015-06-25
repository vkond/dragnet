#include <stdio.h>
#include <math.h>
#include "gsl/gsl_sf_erf.h"
#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_integration.h"

/* Probability density function to be integrated, giving
   cumulative probability */
double pearson4 (double xla, void * mvc_ptr_void)
{
  double *mvc_ptr = (double *)mvc_ptr_void;
  double p4 = exp(*(mvc_ptr+2) - *mvc_ptr*log(1+xla*xla) - *(mvc_ptr+1)*atan(xla));
  return p4;
}

/* Probability density function to be integrated, giving
   cumulative probability */
double pearson6 (double x, void * abc_ptr_void)
{
  double *abc_ptr = (double *)abc_ptr_void;
  double p6 = exp(*(abc_ptr+2) + *abc_ptr*log(x) - *(abc_ptr+1)*log(x+1));
  return p6;
}

int sk_threshold3(int M_int, double s,double Nd, double sk_lims[])
{
  /* This function takes two integers, two doubles and a 2-element
     double array as input, fills the array and gives an integer as
     output */
  /* The function follows PASP papers Nita, Gary & Liu (2007) and Nita
     & Gary (2010) and MNRAS letter Nita & Gary (2010). It calculates
     upper and lower thresholds for the spectral kurtosis (SK) RFI
     estimator. N is the number of FFTs being averaged into power
     spectra (minimum value 1, corresponding to Nyquist-sampled data
     with no averaging) (N is NOT the number of samples per FFT); M is
     number of blocks of N FFTs being used to produce each value of
     the SK estimator (minimum value 3, otherwise equations don't
     work); s is the sigma limit of the thresholds; d gives the type
     of distribution expected from non-RFI (1.0 for standard
     exponential distribution); sk_lims is the array into which the
     calculated threshold values will be placed. */
  
  /* Initial hard-wired variables, which you might want to change */
  int subdiv_lim = 10; /* Maximum number of intervals into which
			  integrations will be subdivided (this can be
			  low for a smoothly varying function such as
			  the probability density function we are
			  integrating) */
  double Ptol = 0.0001; /* Permitted absolute error in cumulative
			    probability values calculated by
			    integration (again, this can be small for
			    our smoothly varying function) */
  int max_its = 1000; /* Maximum number of numerical integrations
			allowed when trying to find an initial or best
			threshold value */
  
  /* Other variables declared */
  int status, ul, fill_lo, fill_hi, n;
  double M = (double)M_int, NN, NN1, M1, MN, MN23, MN45, u1, u2, br, b1, b2, k, r, mvc[3], a, l, u23, u223, rt, alpha, beta, abc[3], delta, P_frac, Ptol_abs,P_thresh, sig, x_thresh, x_lo, x_hi, P, int_abserr;
  gsl_sf_result re_ln_gamma, im_ln_gamma;
  gsl_integration_workspace * int_workspace = gsl_integration_workspace_alloc(subdiv_lim);
  gsl_function p;
  
  /* Constants that depend only on M and N */
  NN = Nd*Nd;
  NN1 = Nd*(Nd+1);
  M1 = M-1;
  MN = M*Nd;
  MN23 = (MN+2)*(MN+3);
  MN45 = (MN+4)*(MN+5);
  /* Constants derived from those above, to be used later */
  u2 = 2*NN1*M*M/M1/MN23;
  br = MN*(Nd+4)-5*Nd-2;
  b1 = 8/NN1/M1*MN23/MN45/MN45*br*br;
  b2 = 3/NN1/M1*MN23/MN45/(MN+6)/(MN+7)*(MN*MN*MN*(Nd+1)+MN*MN*(3*NN+68*Nd+125)-MN*(93*NN+245*Nd+32)+84*NN+48*Nd+24);
  
  /* Decide which probability density function to use */
  k = b1*(b2+3)*(b2+3)/4/(4*b2-3*b1)/(2*b2-3*b1-6);
  
  /* Fail if k<0 */
  if (k<0)
  {
    printf("This function can't calculate thresholds for your chosen parameters! Try increasing M or N.\n");
    sk_lims[0] = 0.0;
    sk_lims[1] = 0.0;
    return 1;
  }
  
  /* Pearson IV if 0<=k<=1 */
  else
  {
    if (k>=0 && k<=1)
    {
      /* More constants derived from those above, dependent ultimately
	 only on M and N */
      u1 = 1.0;
      r = 6*(b2-b1-1)/(2*b2-3*b1-6);
      mvc[0] = (r+2)/2;
      mvc[1] = r*(2-r)*sqrt(b1/(16*(r-1)-b1*(r-2)*(r-2)));
      a = sqrt(u2*(16*(r-1)-b1*(r-2)*(r-2)))/4;
      l = u1 - (r-2)*sqrt(u2*b1)/4;
      delta = 0.0;
      
      /* Constant factor in probability density function, also
	 dependent only on M and N */
      status = gsl_sf_lngamma_complex_e(mvc[0],mvc[1]/2,&re_ln_gamma,&im_ln_gamma);
      mvc[2] = 2*re_ln_gamma.val - gsl_sf_lngamma(2*mvc[0]-1) - log(2)*(2-2*mvc[0]) - log(M_PI); 

      /* Don't try to take mvc[2] out of the integration as a constant
	 factor, because you will have to use exp(mvc[2]), which will
	 be a ridiculously small number */
      
      /* Prepare numerical integrals */
      p.function = &pearson4;
      p.params = &mvc;
    }
    
    /* Pearson VI if k>1 */
    else
    {
      /* More constants derived from those above, dependent only on M
	 and N */
      a = 1.0;
      l = 0.0;
      u23 = M1/M*MN45/4/br;
      u223 = NN1*M/MN23*MN45/2/br;
      rt = 4 + sqrt(16+(4+1/u2)*b1);
      alpha = u23 + u223*(((u223*8-1)*u23+1)*rt+4) - 1; /* alpha in MNRAS letter */
      beta = 3 + 2*rt/b1;
      abc[0] = alpha - 1;
      abc[1] = alpha + beta;
      delta = 1 - alpha/(beta-1);
      
      /* Constant factor in probability density function, also
	 dependent only on M and N */
      abc[2] = -gsl_sf_lnbeta(alpha,beta);
      
      /* Prepare numerical integrals */
      p.function = &pearson6;
      p.params = &abc;
    }
    
    /* Steps common to Pearson IV and VI distributions */
    /* Lower threshold for P that we are aiming for (cumulative probability function) */
    P_frac = gsl_sf_erf(s/sqrt(2));
    /* Absolute tolerance on P-values*/
    Ptol_abs = Ptol*(1-P_frac)/2;
    /* 1-sigma, for Gaussian signal */
    sig = sqrt(u2);
    
    /* Integrations to find upper and lower bounds for SK thresholds by trying out different P thresholds */
    for (ul=-1; ul<=1; ul+=2)
    {
      /* Threshold of cumulative probability distribution that we are aiming for */
      P_thresh = (1-ul*P_frac)/2;
      /* Initial trial values of SK thresholds are symmetric thresholds */
      x_thresh = (1-l)/a - delta + ul*s*sig;
      /* Pearson VI function defined from 0 to infinity, so make sure the trial threshold is not negative in this case */
      if (x_thresh<0 && k>1)
        x_thresh = 0.0;
      fill_lo = 0;
      fill_hi = 0;
      n = 0;
      do
      {
        status = gsl_integration_qagiu(&p,x_thresh,Ptol_abs,0.0,subdiv_lim,int_workspace,&P,&int_abserr);
        if (P>P_thresh)
        {
          x_lo = x_thresh;
          x_thresh += sig; /* Steps of 1-sigma in SK threshold */
          fill_lo = 1;
        }
        else
        {
          x_hi = x_thresh;
          x_thresh -= sig;
          fill_hi = 1;
        }
        n++;
      }
      while ((fill_lo==0 || fill_hi==0) && n<max_its);
      /* Integrations to find lower SK threshold by converging on lower P threshold */
      n = 0;
      do
      {
        x_thresh = (x_lo+x_hi)/2;
        status = gsl_integration_qagiu(&p,x_thresh,Ptol_abs,0.0,subdiv_lim,int_workspace,&P,&int_abserr);
        if (P>P_thresh)
          x_lo = x_thresh;
        else
          x_hi = x_thresh;
        n++;
      }
      while (fabs(P-P_thresh)>Ptol_abs && n<max_its);
      sk_lims[(ul+1)/2] = x_thresh*a + l + delta;
    }
    
    return 0;
  }
  
  gsl_integration_workspace_free (int_workspace);

  return 0;
}

int sk_threshold6(int M_int, float s_float, float d_float, float sk_lims_float[])
{
  /* This function takes two integers, two floats and a 2-element
     float array as input, fills the array and gives an integer as
     output */
  /* The function follows PASP papers Nita, Gary & Liu (2007) and Nita
     & Gary (2010) and MNRAS letter Nita & Gary (2010). It calculates
     upper and lower thresholds for the spectral kurtosis (SK) RFI
     estimator. N is the number of FFTs being averaged into power
     spectra (minimum value 1, corresponding to Nyquist-sampled data
     with no averaging) (N is NOT the number of samples per FFT); M is
     number of blocks of N FFTs being used to produce each value of
     the SK estimator (minimum value 3, otherwise equations don't
     work); s is the sigma limit of the thresholds; d gives the type
     of distribution expected from non-RFI (1.0 for standard
     exponential distribution); sk_lims is the array into which the
     calculated threshold values will be placed. */
  
  /* Initial hard-wired variables, which you might want to change but
     which seem to work well */
  size_t subdiv_lim = 100; /* Maximum number of intervals into which
			      integrations will be subdivided (this
			      can be quite low for a smoothly varying
			      function such as the probability density
			      function we are integrating, but may
			      need to be larger for large M or N) */
  double Ptol = 0.0001; /* Permitted absolute error in cumulative
			   probability values calculated by
			   integration (this can be small for our
			   smoothly varying function) */
  int max_its = 100; /* Maximum number of numerical integrations
			allowed when trying to find an initial or best
			threshold value */
  
  /* Other variables declared */
  int status, ul, fill_lo, fill_hi, n;
  double s=(double)s_float, Nd=(double)d_float, M=(double)M_int, sk_lims[2], NN, NN1, M1, MN, MN23, MN45, u2, u23, u223, B1, B2, B23, br, br1, br2, br3, br32, rt, rt2, sign, k, r, mvc[3], a, l, alpha, beta, abc[3], m11, m21, lamb, Ptol_abs, P_frac, P_thresh, x_thresh, expect, sig, x_lo, x_hi, P, int_abserr;
  gsl_sf_result re_ln_gamma, im_ln_gamma;
  gsl_integration_workspace * int_workspace = gsl_integration_workspace_alloc(subdiv_lim);
  gsl_function p;
  
  if (M_int<2)
    {
      fprintf(stderr, "Thresholds can't be calculated for your chosen parameters! M must be at least 2. Thresholds have been returned as zero.\n");
      sk_lims_float[0] = 0.0;
      sk_lims_float[1] = 0.0;
      gsl_integration_workspace_free(int_workspace);
      return 1;
    }
  
  else
    {    
      /* Constants that depend only on M and N */
      NN = Nd*Nd;
      NN1 = Nd*(Nd+1);
      M1 = M-1;
      MN = M*Nd;
      MN23 = (MN+2)*(MN+3);
      MN45 = (MN+4)*(MN+5);
      /* Constants derived from those above, to be used later */
      u2 = 2*NN1*M*M/M1/MN23;
      br = MN*(Nd+4)-5*Nd-2;
      u23 = M1/M*MN45/4/br;
      B1 = 8/NN1/M1*MN23/MN45/MN45*br*br;
      B2 = 3/NN1/M1*MN23/MN45/(MN+6)/(MN+7)*(MN*MN*MN*(Nd+1)+MN*MN*(3*NN+68*Nd+125)-MN*(93*NN+245*Nd+32)+84*NN+48*Nd+24);
      B23 = B2+3;
      /* Even more constants derived from those above, used by Pearson I and IV */
      br2 = 2*B2-3*B1-6;
      br3 = 10*B2-12*B1-18;
      br32 = br3/br2;
      
      /* Decide which probability density function to use */
      k = B1*B23*B23/4/(4*B2-3*B1)/br2;
      
      /* Pearson IV if 0<=k<=1 */
      if (k>=0 && k<=1)
	{
	  /* More constants derived from those above, dependent ultimately only on M and N */
	  r = br32-2;
	  mvc[0] = br32/2;
	  mvc[1] = r*(2-r)*sqrt(B1/(16*(r-1)-B1*(r-2)*(r-2)));
	  a = 4/sqrt(u2*(16*(r-1)-B1*(r-2)*(r-2))); /* This is 1/a in the paper*/
	  l = ((r-2)/u23/4-1)*a; /* This is -l/a in the paper */
	  
	  /* Constant factor in probability density function, also dependent only on M and N */
	  status = gsl_sf_lngamma_complex_e(mvc[0],mvc[1]/2,&re_ln_gamma,&im_ln_gamma);
	  mvc[2] = 2*re_ln_gamma.val - gsl_sf_lngamma(2*mvc[0]-1) - log(2)*(2-2*mvc[0]) - log(M_PI); /* Don't try to take mvc[2] out of the integration as a constant factor, because you will have to use exp(mvc[2]), which will be a ridiculously small number */
	  
	  /* Prepare numerical integrals */
	  p.function = &pearson4;
	  p.params = &mvc;
	}
      
      /* Pearson VI if k>1 */
      else if (k>1)
	{
	  /* More constants derived from those above, dependent only on M and N */
	  u223 = NN1*M/MN23*MN45/2/br;
	  rt = 4+sqrt(16+(4+1/u2)*B1);
	  alpha = u23+u223*(((u223*8-1)*u23+1)*rt+4)-1; /* alpha in MNRAS letter */
	  beta = 3+2*rt/B1;
	  abc[0] = alpha-1;
	  abc[1] = alpha+beta;
	  a = 1.0;
	  l = alpha/(beta-1)-1; /* This is -delta in the paper */
	  
	  /* Constant factor in probability density function, also dependent only on M and N */
	  abc[2] = -gsl_sf_lnbeta(alpha,beta);
	  
	  /* Prepare numerical integrals */
	  p.function = &pearson6;
	  p.params = &abc;
	}
      
      /* Pearson I if k<0*/
      else
	{
	  br1 = B23/u23;
	  rt2 = fabs(br3)/br3*sqrt(br1*br1+4*u2*(3*B1-4*B2)*br2);
	  sign = fabs(br32)/br32;
	  m11 = 1-br32/2-(1+br32/2)*sign*br1/rt2; /* This is m1+1 in the paper*/
	  m21 = 2-br32-m11; /* This is m2+1 in the paper*/
	  a = sign*br2/rt2;
	  l = (sign*br1/rt2+1)/2;
	  lamb=br1/(B1-B2+1)/3+1;
	}
      
      /* Steps common to Pearson I, IV and VI distributions */
      /* Lower threshold for P that we are aiming for (cumulative probability function) */
      P_frac = gsl_sf_erf(s/sqrt(2));
      /* Absolute tolerance on P-values*/
      Ptol_abs = Ptol*(1-P_frac)/2;
      
      /* Integrations to find upper and lower bounds for SK thresholds by trying out different P thresholds */
      for (ul=-1; ul<=1; ul+=2)
	{
	  /* Threshold of cumulative probability distribution that we are aiming for */
	  P_thresh = (1-ul*P_frac)/2;
	  /* Expectation value of integrating variable, for Gaussian signal */
	  expect = a+l;
	  /* 1-sigma deviation in integrating variable, for Gaussian signal */
	  sig = sqrt(u2)*fabs(a);
	  /* Initial trial values of SK thresholds are symmetric thresholds */
	  x_thresh = expect + ul*s*sig;
	  /* Pearson I and VI functions not defined below 0, so make sure the trial threshold is not negative in these cases */
	  if ((k<0 || k>1) && x_thresh<0)
	    x_thresh = 0.0;
	  /* Pearson I function not defined above 1, so make sure the trial threshold is not above 1 in this case */
	  if (k<0 && x_thresh>1)
	    x_thresh = 1.0;
	  fill_lo = 0;
	  fill_hi = 0;
	  n = 0;
	  do
	    {
	      if (k>=0)
		{
		  status = gsl_integration_qagiu(&p,x_thresh,Ptol_abs,0.0,subdiv_lim,int_workspace,&P,&int_abserr);
		}
	      else
		{
		  P = 1-gsl_sf_beta_inc(m11,m21,x_thresh); /* This incomplete beta function is already normalised by the complete one */
		}
	      if (P>P_thresh)
		{
		  x_lo = x_thresh;
		  x_thresh += sig; /* Steps of 1-sigma in SK threshold */
		  if (k<0 && x_thresh>1)
		    x_thresh = 1.0;
		  fill_lo = 1;
		}
	      else
		{
		  x_hi = x_thresh;
		  x_thresh -= sig;
		  if ((k<0 || k>1) && x_thresh<0)
		    x_thresh = 0.0;
		  fill_hi = 1;
		}
	      n++;
	    }
	  while ((fill_lo==0 || fill_hi==0) && n<max_its);
	  /* Integrations to find SK thresholds by converging on P thresholds */
	  n = 0;
	  do
	    {
	      x_thresh = (x_lo+x_hi)/2;
	      if (k>=0)
		{
		  status = gsl_integration_qagiu(&p,x_thresh,Ptol_abs,0.0,subdiv_lim,int_workspace,&P,&int_abserr);
		}
	      else
		{
		  P = 1-gsl_sf_beta_inc(m11,m21,x_thresh);
		}
	      if (P>P_thresh)
		x_lo = x_thresh;
	      else
		x_hi = x_thresh;
	      n++;
	    }
	  while (fabs(P-P_thresh)>Ptol_abs && n<max_its);
	  sk_lims[(ul+1)/2] = (x_thresh-l)/a+lamb;
	}
      
      sk_lims_float[0] = (float)sk_lims[0];
      sk_lims_float[1] = (float)sk_lims[1];
      
      gsl_integration_workspace_free(int_workspace);
      return 0;
    }
}

/* Compute the spectral kurtosis mask. Input is the dynamic spectrum z
   with nx channels and ny time samples.  Output is the integer mask
   of size nx channels and my intervals, where each interval is m time
   samples. n is the number of averaged spectra (n=6 for
   491.52us=6x81.92us LOTAAS sampling). d is the averaging
   factor. Since two polarizations are added, d=2.0. skmin and skmax
   are the sk limits computed by sk_thresh3.  The mask is filled with
   intervals to mask (value 1) or keep (value 0). The masked values in
   the dynamic spectrum will be replaced by the average of at least 10
   adjacent intervals in frequency (5 on either side).
 */
int compute_sk_mask(float *z,int nx,int ny,int my,int m,float nd,float skmin,float skmax,int *mask)
{
  int i,j,k,l,i0,mact,i1,navg=5,nmask=0;
  double s1,s2,sk;
  float zmax;
  FILE *file;

  // Compute mask
  for (i=0;i<nx;i++) {
    for (j=0;j<my;j++) {
      for (k=0,s1=0.0,s2=0.0,mact=0;k<m;k++) {
	l=i+nx*(j*m+k);
	if (j*m+k>=ny)
	  continue;
	s1+=z[l];
	s2+=z[l]*z[l];
	mact++;
      }
      sk=(mact*nd+1.0)/(mact-1.0)*(mact*s2/(s1*s1)-1.0);

      // Mask block
      if (sk<skmin || sk>skmax || isnan(sk)) {
      	mask[i+nx*j]=1;
	nmask++;
      } else {
      	mask[i+nx*j]=0;
      }
    }
  }

  // Fill masked values
  for (i=0;i<nx;i++) {
    for (j=0;j<my;j++) {
      // Skip unmasked intervals
      if (mask[i+nx*j]==0) 
	continue;

      // Compute average value of adjacent channels
      s1=0.0;
      s2=0.0;

      // Loop over time samples
      for (k=0;k<m;k++) {
	// Skip out of array points
	if (j*m+k>=ny)
	  continue;

	// Right side
	for (i0=i+1,zmax=0.0,i1=0;i1<navg;i0++) {
	  // Skip out of array channels
	  if (i0>=nx)
	    break;
	  // Skip masked intervals
	  if (mask[i0+nx*j]==1) 
	    continue;
	  // Index
	  l=i0+nx*(j*m+k);

	  // Keep track of maximum value
	  if (z[l]>zmax)
	    zmax=z[l];
	  s1+=z[l];
	  s2+=1.0;

	  // Increment counter
	  i1++;
	}

	// Substract maximum value
	s1-=zmax;
	s2-=1.0;

	// Left side
	for (i0=i-1,zmax=0.0,i1=0;i1<navg;i0--) {
	  // Skip out of array channels
	  if (i0<0)
	    break;
	  // Skip masked intervals
	  if (mask[i0+nx*j]==1) 
	    continue;
	  // Index
	  l=i0+nx*(j*m+k);

	  // Keep track of maximum value
	  if (z[l]>zmax)
	    zmax=z[l];
	  s1+=z[l];
	  s2+=1.0;

	  // Increment counter
	  i1++;
	}

	// Substract maximum value
	s1-=zmax;
	s2-=1.0;
      }

      // Fill masked values
      for (k=0;k<m;k++) {
	l=i+nx*(j*m+k);
	// Skip out of array points
	if (j*m+k>=ny)
	  continue;
	z[l]=s1/s2;
      }
    }
  }

  return nmask;
}
