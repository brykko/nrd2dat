#ifndef _FILTERS_H_INCLUDED_
#define _FILTERS_H_INCLUDED_

#include <cmath>

class FilterButterworth
{
	/// <summary>
	/// rez amount, from sqrt(2) to ~ 0.1
	/// </summary>
public:

	float resonance;
	float frequency;
	int sampleRate;
	float c, a1, a2, a3, b1, b2;

	float inputHistory[2];
	float outputHistory[3];

	const double PI = 3.141592653589793;

public:

	FilterButterworth(float freq, int sampRate, float res)
	{

		/// <summary>
		/// Array of input values, latest are in front
		/// </summary>

		/// <summary>
		/// Array of output values, latest are in front
		/// </summary>

		resonance = res;
		frequency = freq;
		sampleRate = sampRate;

		c = (float)_CMATH_::tan(PI * frequency / (float)sampleRate);
		a1 = 1.0f / (1.0f + resonance * c + c * c);
		a2 = -2.0f * a1;
		a3 = a1;
		b1 = 2.0f * (c * c - 1.0f) * a1;
		b2 = (1.0f - resonance * c + c * c) * a1;

		//cout << "Filter params: c " << c << ", a1 " << a1 << ", a2" << a2 << ", a3 " << a3 << ", b1" << b1 << ", b2" << b2 << endl;

		inputHistory[0] = 0;
		inputHistory[1] = 0;
		outputHistory[0] = 0;
		outputHistory[1] = 0;
		outputHistory[2] = 0;

	}


	void Update(float newInput)
	{
		float newOutput = a1 * newInput + a2 * inputHistory[0] + a3 * inputHistory[1] - b1 * outputHistory[0] - b2 * outputHistory[1];

		///cout << "New value = " << newOutput << endl;

		inputHistory[1] = inputHistory[0];
		inputHistory[0] = newInput;

		outputHistory[2] = outputHistory[1];
		outputHistory[1] = outputHistory[0];
		outputHistory[0] = newOutput;

		//cout << "New value = " << outputHistory[0] << endl;

	}

	float Value()
	{
		return outputHistory[0];
	}

};

#endif