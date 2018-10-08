#define _USE_MATH_DEFINES


using namespace std;

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "image_comps.h"
#include "LoG_Filter.h"


float LoG_1(double s, double sigma)
{
	/*   g(s)   */

	double constant = 2 * M_PI * sigma*sigma;
	constant = 1.0f / constant;
	double retval = constant * ((s*s - sigma * sigma) / pow(sigma, 4)) * exp((-1 * s*s) / (2 * sigma*sigma));

	return (float)retval;
}

float LoG_2(double s, double sigma)
{
	return exp((-1 * s*s) / (2 * sigma*sigma));
}

#define FILTER_EXTENT H
#define FILTER_TAPS (2*FILTER_EXTENT+1)

void my_LoG::init(float sigma, float alpha, int H)
{
	this->H = H;
	this->sigma = sigma;
	this->alpha = alpha;
	int r, c;
	g1 = new float[FILTER_TAPS];
	float *g1_m = g1 + FILTER_EXTENT;
	g2 = new float[FILTER_TAPS];
	float *g2_m = g2 + FILTER_EXTENT;

	for (r = -FILTER_EXTENT; r <= FILTER_EXTENT; r++)
	{
		g1_m[r] = LoG_1(r, sigma);
		g2_m[r] = LoG_2(r, sigma);
		
	}
#ifdef DEBUG
	printf("Filter coefficients are: \n");
	for (r = -FILTER_EXTENT; r <= FILTER_EXTENT; r++)
	{
		cout << g1_m[r] << " ";
	}
	printf("\n\n");
	for (r = -FILTER_EXTENT; r <= FILTER_EXTENT; r++)
	{
		cout << g2_m[r] << " ";
	}
	printf("\n");
#endif

}


void my_LoG::apply_dilation(my_image_comp* in, my_image_comp * out,  int * A, int n) {
	int r, c, x, y;
	float *ip, *op;
	int *A_m = A + n * (n / 2) + (n / 2);
	int H1 = n / 2;
	for (r = 0; r < out->height; r++) {
		for (c = 0; c < out->width; c++) {
			op = out->buf + r * out->stride + c;
			ip = in->buf + r * in->stride + c;

			if ( ip[0]*A_m[0] == 256 )
			{
				// Apply dilation using OR

				for (x = -H1; x <= H1; x++) {
					for (y = -H1; y <= H1; y++) {
						unsigned int val = (unsigned int)ip[x*out->stride + y]; // should be either 0 or 1<<8
						unsigned int A_val = A_m[x*n + y] << 8;
						A_val = A_val | val;
						op[x*out->stride + y] = (float)A_val;
					}
				}

			}
		}
	}
}

void my_LoG::apply_filter(my_image_comp *in, my_image_comp* out) {
	int r, c;
	float *ip, *op, *mp, *g1_m = g1 + FILTER_EXTENT, *g2_m = g2 + FILTER_EXTENT;
	my_image_comp intermediate;
	intermediate.init(in->height, in->width, in->border);

	// horizontal filter first
	for (r = 0; r < intermediate.height; r++) {
		for (c = 0; c < intermediate.width; c++) {
			mp = intermediate.buf + r * intermediate.stride + c;
			ip = in->buf + r * in->stride + c;
			float sum = 0.0F;
			for (int x = -FILTER_EXTENT; x <= FILTER_EXTENT; x++)
			{
				sum += ip[x] * g1_m[x];
			}
			*mp = sum;
		}
	}
	intermediate.perform_boundary_extension();

	// vertical filter now
	for (r = 0; r < out->height; r++) {
		for (c = 0; c < out->width; c++) {
			op = out->buf + r * out->stride + c;
			mp = intermediate.buf + r * intermediate.stride + c;
			float sum;
			sum = 0.0F;
			for (int x = -FILTER_EXTENT; x <= FILTER_EXTENT; x++) {
				sum += mp[intermediate.stride*x] * g2_m[x];
			}
			*op = sum;
		}
	}
	// horizontal filter first
	for (r = 0; r < intermediate.height; r++) {
		for (c = 0; c < intermediate.width; c++) {
			mp = intermediate.buf + r * intermediate.stride + c;
			ip = in->buf + r * in->stride + c;
			float sum = 0.0F;
			for (int x = -FILTER_EXTENT; x <= FILTER_EXTENT; x++)
			{
				sum += ip[x] * g2_m[x];
			}
			*mp = sum;
		}
	}
	intermediate.perform_boundary_extension();
	// vertical filter now
	for (r = 0; r < out->height; r++) {
		for (c = 0; c < out->width; c++) {
			op = out->buf + r * out->stride + c;
			mp = intermediate.buf + r * intermediate.stride + c;
			float sum;
			sum = 0.0F;
			for (int x = -FILTER_EXTENT; x <= FILTER_EXTENT; x++) {
				sum += mp[intermediate.stride*x] * g1_m[x];
			}
			*op += sum;
			// Scaling plus thresholding
			*op *= alpha;
			*op *= alpha;
			*op += 128;
			//if (*op > 3 * alpha)
			//	*op = 256;
			//else if (*op < 7 * alpha)
			//	*op = 0;
		}
	}
	
	for (r = 0; r < out->height; r++) {
		for (c = 0; c < out->width; c++) {
			op = out->buf + r * out->stride + c;
			if (*op < 20 && *op > 10 && r > 0 && c > 0 && r < out->height -1 && c < out->width-1)
			{
				// Check surrounding pixels for an edge
				for (int x = -1; x <= 1; x++) {
					for (int y = -1; y <= 1; y++) {
						if (op[x*out->stride + y] == 255) {
							*op = 256;
							break;
						}
					}
				}
			
			}
		}
	}


}
