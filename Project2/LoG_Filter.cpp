#define _USE_MATH_DEFINES


using namespace std;

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "image_comps.h"
#include "LoG_Filter.h"


float LoG(double s1, double s2, double sigma)
{
	/*   g(s)   */

	double constant = M_PI * pow(sigma, 4);
	constant = 1.0f / constant;
	double sigma_2 = sigma * sigma;
	constant = constant * (((s1*s1 + s2*s2)/(2*sigma_2)) - 1);

	double exponent = exp(-(s1*s1 + s2*s2) / (2*sigma_2));

	double retval = constant * exponent;
	
	return (float)retval;

}

#define FILTER_EXTENT H
#define FILTER_DIM (2*FILTER_EXTENT+1)
#define FILTER_TAPS (FILTER_DIM*FILTER_DIM) 


void my_LoG::init(float sigma, float alpha, int H)
{
	this->H = H;
	this->sigma = sigma;
	this->alpha = alpha;
	int r, c;
	g = new float[FILTER_TAPS];
	float *g_m = g + (FILTER_DIM * FILTER_EXTENT) + FILTER_EXTENT;
	for (r = -FILTER_EXTENT; r <= FILTER_EXTENT; r++)
	{
		for (c = -FILTER_EXTENT; c <= FILTER_EXTENT; c++)
		{
			g_m[r*FILTER_DIM + c] = LoG(r, c, sigma);
		}
	}
#ifdef DEBUG
	printf("Filter coefficients are: \n");
	for (r = -FILTER_EXTENT; r <= FILTER_EXTENT; r++)
	{
		for (c = -FILTER_EXTENT; c <= FILTER_EXTENT; c++)
		{
			//printf("%f ", g_m[r*FILTER_DIM + c]);
			std::cout << g_m[r*FILTER_DIM + c] << " ";
		}
		printf("\n");
	}
#endif


}



void my_LoG::apply_filter(my_image_comp *in, my_image_comp* out) {
	int r, c;
	float *ip, *op, *g_m = g + (FILTER_DIM*FILTER_EXTENT) + FILTER_EXTENT;

	for (r = 0; r < out->height; r++) {
		for (c = 0; c < out->width; c++) {
			// in:out should be 1:1
			ip = in->buf + r * in->stride + c;
			op = out->buf + r * out->stride + c;
			float sum = 0.0F;
			for (int y = -FILTER_EXTENT; y <= FILTER_EXTENT; y++) {
				for (int x = -FILTER_EXTENT; x <= FILTER_EXTENT; x++) {
					//printf("multiply %f and %f\n", ip[y*in->stride + x], g_m[y*FILTER_DIM + x]);
					sum += ip[y*in->stride + x] * g_m[y*FILTER_DIM + x];
				}
			}
			// Now scale and offset the output
			// Samples are naturally centered at 0. Add 128
			sum = sum * alpha;
			sum += 128;
			*op = sum;
			
			//cout << "sum is " << sum << "alpha is " << alpha << endl;
		}
	}
}