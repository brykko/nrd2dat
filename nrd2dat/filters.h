#ifndef _FILTERS_H_INCLUDED_
#define _FILTERS_H_INCLUDED_

#include <cmath>

#define PI 3.141592653589793

class FilterButterworth
{

public:

	float resonance;
	float frequency;
	int sampleRate;
	double c, a1, a2, a3, b1, b2;

	double inputHistory[2];
	double outputHistory[3];

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

		c = _CMATH_::tan(PI * frequency / (double)sampleRate);
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


	void Update(double newInput)
	{
		double newOutput = a1 * newInput + a2 * inputHistory[0] + a3 * inputHistory[1] - b1 * outputHistory[0] - b2 * outputHistory[1];

		///cout << "New value = " << newOutput << endl;

		inputHistory[1] = inputHistory[0];
		inputHistory[0] = newInput;

		outputHistory[2] = outputHistory[1];
		outputHistory[1] = outputHistory[0];
		outputHistory[0] = newOutput;

		//cout << "New value = " << outputHistory[0] << endl;

	}

	double Value()
	{
		return outputHistory[0];
	}

};


class ExpMovingAverage{

protected:
	float mValue;
	float mFrequency;
	double a, b;

public:

	ExpMovingAverage(float f) {
		init(f, 0);
	}

	ExpMovingAverage(float f, float initVal) {
		init(f, initVal);
	}

	void init(float f, float initVal) {
		mValue = initVal;
		mFrequency = f;
		double alpha = 1.0 - exp(-PI*f);
		b = alpha;      // coefficient of current sample
		a = 1.0 - alpha;  // coefficient of adjacent samples
	}

	void update(float y) {
		mValue = b*(double)y + a*(double)mValue;
	}
	
	float value() {
		return mValue;
	}

	float timeConstant() {
		return 1 / (2 * PI*mFrequency);
	}
};

#endif