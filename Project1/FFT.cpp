#define _USE_MATH_DEFINES
#include "FFT.h"
#include <stdlib.h>
#include <stdio.h>
#include <cmath>

FFT::FFT(void)
{
	FFT_bit_table = new int *[Max_Fast_Bits];

   int len = 2;
   for (int b = 1; b <= Max_Fast_Bits; b++) {

      FFT_bit_table[b - 1] = new int[len];
      for (int i = 0; i < len; i++)
         FFT_bit_table[b - 1][i] = Reverse_Bits(i, b);
	  

      len <<= 1;
   }
}


FFT::~FFT(void)
{
	if (FFT_bit_table) {
      for (int b = 1; b <= Max_Fast_Bits; b++) {
         delete[] FFT_bit_table[b-1];
      }
      delete[] FFT_bit_table;
   }
}

void FFT::Window_Func(int function_no, int num_samples, float *in)
{
	switch (function_no) {
		default:
		  fprintf(stderr, "FFT::WindowFunc - Invalid window function: %d\n", function_no);
		  break;
		case FFT::Rectangular:
		  break;
		case FFT::Bartlett:
		{
		  // Bartlett (triangular)
		  const int nPairs = (num_samples - 1) / 2; 
		  const float denom = num_samples / 2.0f;
		  in[0] = 0.0f;
		  for (int ii = 1;
			   ii <= nPairs; 
			   ++ii) {
			 const float value = ii / denom;
			 in[ii] *= value;
			 in[num_samples - ii] *= value;
		}
	   }
		case FFT::Hamming:
	   {
		  // Hamming
		  const double multiplier = 2 * M_PI / (num_samples-1);
		  static const double coeff0 = 0.54, coeff1 = -0.46;
		  for (int ii = 0; ii < num_samples; ++ii)
			in[ii] *= coeff0 + coeff1 * cos(ii * multiplier);
				
	   }
		  break;
		  case FFT::Hanning:
	   {
		  // Hanning
		  const double multiplier = 2.0 * M_PI / num_samples;
		  static const double coeff0 = 0.5, coeff1 = -0.5;
		  for (int ii = 0; ii < num_samples; ++ii)
			 in[ii] *= coeff0 + coeff1 * cos(ii * multiplier);
	   }
	}
}


int FFT::Reverse_Bits(int index, int Num_Bits)
{
   int i, rev;

   for (i = rev = 0; i < Num_Bits; i++) {
      rev = (rev << 1) | (index & 1);
      index >>= 1;
   }
   return rev;
}

int FFT::Is_Power_Of_Two(int x)
{
	if (x < 2)
      return false;

   if (x & (x - 1))
      return false;

   return true;
}

int FFT::Number_Of_Bits_Needed(int Power_Of_Two)
{
	int i;

   if (Power_Of_Two < 2) {
      fprintf(stderr, "Error: FFT called with size %d\n", Power_Of_Two);
      exit(1);
   }

   for (i = 0;; i++)
      if (Power_Of_Two & (1 << i))
         return i;
}

inline int FFT::Fast_Reverse_Bits(int i, int Num_Bits)
{
	if (Num_Bits <= Max_Fast_Bits)
      return FFT_bit_table[Num_Bits - 1][i];
   else
      return Reverse_Bits(i, Num_Bits);
}


void FFT::Transform(int num_samples, bool inverse, float *real_in, float *imag_in, float *real_out, float *imag_out)
{
	 int NumBits;                 /* Number of bits needed to store indices */
   int i, j, k, n;
   int BlockSize, BlockEnd;

   double angle_numerator = 2.0 * M_PI;
   double tr, ti;                /* temp real, temp imaginary */

   if (!Is_Power_Of_Two(num_samples)) {
      fprintf(stderr, "%d is not a power of two\n", num_samples);
      exit(1);
   }



   if (!inverse)
      angle_numerator = -angle_numerator;

   NumBits = Number_Of_Bits_Needed(num_samples);

   /*
    **   Do simultaneous data copy and bit-reversal ordering into outputs...
    */

   for (i = 0; i < num_samples; i++) {
      j = Fast_Reverse_Bits(i, NumBits);
	  //printf("[%d]<->[%d]\n",i,j);
      real_out[j] = real_in[i];
      imag_out[j] = (imag_in == NULL) ? 0.0 : imag_in[i];
   }
   //getchar();
   /*
    **   Do the FFT itself...
    */


   BlockEnd = 1;
   for (BlockSize = 2; BlockSize <= num_samples; BlockSize <<= 1) {

      double delta_angle = angle_numerator / (double) BlockSize;

      double sm2 = sin(-2 * delta_angle);
      double sm1 = sin(-delta_angle);
      double cm2 = cos(-2 * delta_angle);
      double cm1 = cos(-delta_angle);
      double w = 2 * cm1;
      double ar0, ar1, ar2, ai0, ai1, ai2;

      for (i = 0; i < num_samples; i += BlockSize) {
         ar2 = cm2;
         ar1 = cm1;

         ai2 = sm2;
         ai1 = sm1;

         for (j = i, n = 0; n < BlockEnd; j++, n++) {
            ar0 = w * ar1 - ar2;
            ar2 = ar1;
            ar1 = ar0;

            ai0 = w * ai1 - ai2;
            ai2 = ai1;
            ai1 = ai0;

            k = j + BlockEnd;
            tr = ar0 * real_out[k] - ai0 * imag_out[k];
            ti = ar0 * imag_out[k] + ai0 * real_out[k];

            real_out[k] = real_out[j] - tr;
            imag_out[k] = imag_out[j] - ti;

            real_out[j] += tr;
            imag_out[j] += ti;

			
         }
		 

      }

      BlockEnd = BlockSize;
   }

   /*
      **   Need to normalize if inverse transform...
    */

   if (inverse) {
      float denom = (float) num_samples;

      for (i = 0; i < num_samples; i++) {
         real_out[i] /= denom;
         imag_out[i] /= denom;
      }
   }
}

void FFT::Power_Spectrum(int num_samples, float *in, float *out)
{
	int Half = num_samples / 2;
   int i;

   float theta = M_PI / Half;

   float *tmp_real = new float[Half];
   float *tmp_imag = new float[Half];
   float *real_out = new float[Half];
   float *imag_out = new float[Half];

   for (i = 0; i < Half; i++) {
      tmp_real[i] = in[2 * i];
      tmp_imag[i] = in[2 * i + 1];
   }

   Transform(Half, 0, tmp_real, tmp_imag, real_out, imag_out);

   float wtemp = float (sin(0.5 * theta));

   float wpr = -2.0 * wtemp * wtemp;
   float wpi = -1.0 * float (sin(theta));
   float wr = 1.0 + wpr;
   float wi = wpi;

   int i3;

   float h1r, h1i, h2r, h2i, rt, it;

   for (i = 1; i < Half / 2; i++) {

      i3 = Half - i;

      h1r = 0.5 * (real_out[i] + real_out[i3]);
      h1i = 0.5 * (imag_out[i] - imag_out[i3]);
      h2r = 0.5 * (imag_out[i] + imag_out[i3]);
      h2i = -0.5 * (real_out[i] - real_out[i3]);

      rt = h1r + wr * h2r - wi * h2i;
      it = h1i + wr * h2i + wi * h2r;

      out[i] = rt * rt + it * it;

      rt = h1r - wr * h2r + wi * h2i;
      it = -h1i + wr * h2i + wi * h2r;

      out[i3] = rt * rt + it * it;

      wr = (wtemp = wr) * wpr - wi * wpi + wr;
      wi = wi * wpr + wtemp * wpi + wi;
   }

   rt = (h1r = real_out[0]) + imag_out[0];
   it = h1r - imag_out[0];
   out[0] = rt * rt + it * it;

   rt = real_out[Half / 2];
   it = imag_out[Half / 2];
   out[Half / 2] = rt * rt + it * it;

   delete[]tmp_real;
   delete[]tmp_imag;
   delete[]real_out;
   delete[]imag_out;

}




