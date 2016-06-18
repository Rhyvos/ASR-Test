#pragma once
class FFT
{
public:
	FFT(void);
	void Transform(int num_samples, bool inverse, float *real_in, float *imag_in, float *real_out, float *imag_out);
	void Power_Spectrum(int num_samples, float *in, float *out);
	~FFT(void);
	void Window_Func(int function_no, int num_samples, float *data);
	enum eWindow_Func
	{
		Rectangular,
		Bartlett,
		Hamming,
		Hanning
	};
private:
	int Reverse_Bits(int index, int Num_Bits);
	int Is_Power_Of_Two(int x);
	int Number_Of_Bits_Needed(int Power_Of_Two);
	int Fast_Reverse_Bits(int i, int Num_Bits);

	int ** FFT_bit_table;
	static const int Max_Fast_Bits = 16;
};










