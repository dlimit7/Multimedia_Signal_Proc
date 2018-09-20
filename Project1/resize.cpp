#define _USE_MATH_DEFINES

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "resize.h"

float windowed_sinc(double s, int extent, double alpha)
{
	double t = extent + 0.5;
	double temp, sinc, window_val;
	if (s == 0.0) return 1.0f;
	else {
		window_val = (0.5)*(1.0 + cos((M_PI*s) / (t)));
		if (alpha < 1) 
		{
			temp = M_PI * s;
			if (temp == 0.0f) return 1;
			sinc = sin(temp) / (temp);
		}
		else // presumable alpha is 5/3
		{
			temp = M_PI * s * (1 / alpha);
			if (temp == 0.0f) return 1;
			sinc = sin(temp) / (temp);
		}	
		return (float) (sinc * window_val);
	}
}

#define FILTER_EXTENT H
#define FILTER_DIM (2*FILTER_EXTENT+1) 
#define FILTER_TAPS (FILTER_DIM*FILTER_DIM)

void my_resizer::init(int H, bool is_expand, bool is_sinc_interp)
{

	this->H = H;
	this->is_expand = is_expand;
	this->is_sinc_interp = is_sinc_interp;
	int i, r, c;
	float *q_m[25], *g_m[5];
	for (i = 0; i < 25; i++) {
		q[i] = new float[FILTER_TAPS];
		q_m[i] = q[i] + (FILTER_DIM*FILTER_EXTENT) + FILTER_EXTENT; // q_m is centered in the middle
	}
	for (i = 0; i < 5; i++) {
		g[i] = new float[FILTER_DIM];
		g_m[i] = g[i] + FILTER_EXTENT;
	}
	double alpha;
	if (is_expand) 
	{
		alpha = 3.0f / 5;
		if (is_sinc_interp)
		{
			// Make the windowed sinc kernels
			float bank[5] = { 0, 0.4f, -0.2f, 0.2f, -0.4f };
			for (i = 0; i < 25; i++) {
				for (r = -FILTER_EXTENT; r <= FILTER_EXTENT; r++) {
					for (c = -FILTER_EXTENT; c <= FILTER_EXTENT; c++) {
						q_m[i][r*FILTER_DIM + c] = (float)windowed_sinc(r + bank[i / 5], H, alpha) * windowed_sinc(c + bank[i % 5], H, alpha);
					}
				}
			}
			for (i = 0; i < 5; i++) {
				for (r = -FILTER_EXTENT; r <= FILTER_EXTENT; r++) {
					g_m[i][r] = windowed_sinc(r + bank[i % 5], H, alpha);
				}
			}
#ifdef DEBUG
			for (i = 0; i < 5; i++) {
				printf("i = %d: ", i);
				for (r = -FILTER_EXTENT; r <= FILTER_EXTENT; r++) {
					printf("%f ", g_m[i][r]);
				}
				printf("\n");
			}
			for (r = -FILTER_EXTENT; r <= FILTER_EXTENT; r++) {
				for (c = -FILTER_EXTENT; c <= FILTER_EXTENT; c++) {
					printf("%f ", q_m[0][r*FILTER_DIM + c]);
				}
				printf("\n");
			}
			printf("\n\n");
#endif
			// Normalising
			float gain;
			for (i = 0; i < 25; i++) { // For each kernel..
				gain = 0;
				// Find their DC Gain
				for (r = -FILTER_EXTENT; r <= FILTER_EXTENT; r++) {
					for (c = -FILTER_EXTENT; c <= FILTER_EXTENT; c++) {
						gain += q_m[i][r*FILTER_DIM + c];
					}
				}
				// Make DC Gain = 1
				gain = 1.0f / gain;
				for (r = -FILTER_EXTENT; r <= FILTER_EXTENT; r++) {
					for (c = -FILTER_EXTENT; c <= FILTER_EXTENT; c++) {
						q_m[i][r*FILTER_DIM + c] = q_m[i][r*FILTER_DIM + c] * gain;
					}
				}
			}
#ifdef DEBUG
			// Debugging
			for (i = 0; i < 25; i++) {
				printf(" i = %d\n", i);
				for (r = -FILTER_EXTENT; r <= FILTER_EXTENT; r++) {
					for (c = -FILTER_EXTENT; c <= FILTER_EXTENT; c++) {
						printf("%f ", q_m[i][r*FILTER_DIM + c]);
					}
					printf("\n");
				}
				printf("\n\n");
			}
			/*
			// Recalculate the gain to check if 1
			for (i = 0; i < 25; i++) { // For each kernel..
				gain = 0;
				// Find their DC Gain
				for (r = -FILTER_EXTENT; r <= FILTER_EXTENT; r++) {
					for (c = -FILTER_EXTENT; c <= FILTER_EXTENT; c++) {
						gain += q_m[i][r*FILTER_DIM + c];
					}
				}
				printf("Gain is %f\n", gain);
			}
			*/
#endif
			
		} else // Bilinear interp
		{
			/*  PSF is
			*	
			*	s1*s2			(1-s2)*s1
			*
			*	s2*(1-s1)	(1-s2)*(1-s1)
			*/
			
		}
	}
	else // Reduction
	{
		alpha = 5.0f / 3;
		float constant = (1.0f / 2) * (1.0f / alpha) * (1.0f / alpha);
		// Make the windowed sinc kernels
		float bank[3] = { 0, 1.0f / 3, -1.0f / 3 };
		for (i = 0; i < 9; i++) {
			for (r = -FILTER_EXTENT; r <= FILTER_EXTENT; r++) {
				for (c = -FILTER_EXTENT; c <= FILTER_EXTENT; c++) {
					q_m[i][r*FILTER_DIM + c] = (float)windowed_sinc(r + bank[i / 3], H, alpha) * windowed_sinc(c + bank[i % 3], H, alpha) * constant;
				}
			}
		}
		// Normalising
		float gain;
		for (i = 0; i < 9; i++) { // For each kernel..
			gain = 0;
			// Find their DC Gain
			for (r = -FILTER_EXTENT; r <= FILTER_EXTENT; r++) {
				for (c = -FILTER_EXTENT; c <= FILTER_EXTENT; c++) {
					gain += q_m[i][r*FILTER_DIM + c];
				}
			}
			// Make DC Gain = 1
			gain = 1.0f / gain;
			for (r = -FILTER_EXTENT; r <= FILTER_EXTENT; r++) {
				for (c = -FILTER_EXTENT; c <= FILTER_EXTENT; c++) {
					q_m[i][r*FILTER_DIM + c] = q_m[i][r*FILTER_DIM + c] * gain;
				}
			}
		}
#ifdef DEBUG
		// Debugging
		for (i = 0; i < 9; i++) {
			printf(" i = %d\n", i);
			for (r = -FILTER_EXTENT; r <= FILTER_EXTENT; r++) {
				for (c = -FILTER_EXTENT; c <= FILTER_EXTENT; c++) {
					printf("%f ", q_m[i][r*FILTER_DIM + c]);
				}
				printf("\n");
			}
			printf("\n\n");
		}
		// Recalculate the gain to check if 1
		for (i = 0; i < 9; i++) { // For each kernel..
			gain = 0;
			// Find their DC Gain
			for (r = -FILTER_EXTENT; r <= FILTER_EXTENT; r++) {
				for (c = -FILTER_EXTENT; c <= FILTER_EXTENT; c++) {
					gain += q_m[i][r*FILTER_DIM + c];
				}
			}
			printf("Gain is %f\n", gain);
		}
#endif
	}

}


void my_resizer::apply_filter(my_image_comp *in, my_image_comp* out) {
	int r, c, i;

	float *q_m[25], *g_m[5];
	for (i = 0; i < 25; i++) {
		q_m[i] = q[i] + (FILTER_DIM*FILTER_EXTENT) + FILTER_EXTENT; // q_m is centered in the middle
	}
	for (i = 0; i < 5; i++) {
		g_m[i] = g[i] + FILTER_EXTENT;
	}

	// Check for consistent dimensions
	//assert(in->border >= FILTER_EXTENT);
	//assert((out->height <= in->height) && (out->width <= in->width));

	if (this->is_expand)
	{
		if (this->is_sinc_interp) {
			int m1, m2, n1, n2, x_r, x_c;
			float *ip, *op;
			// Perform the convolution
#ifdef DEBUG
			printf("Performing convolution: Output image %dx%d\n", out->height, out->width);
#endif
			for (r = 0; r < out->height; r++) {
				for (c = 0; c < out->width; c++) {
					m1 = r / 5; m2 = c / 5;
					n1 = r % 5; n2 = c % 5;
					switch (n1)
					{
					case 0:
						x_r = 3 * m1; break;
					case 1:
					case 2:
						x_r = 3 * m1 + 1; break;
					case 3:
					case 4:
						x_r = 3 * m1 + 2; break;
					default:
						perror(NULL);
					}
					switch (n2)
					{
					case 0:
						x_c = 3 * m2; break;
					case 1:
					case 2:
						x_c = 3 * m2 + 1; break;
					case 3:
					case 4:
						x_c = 3 * m2 + 2; break;
					default:
						perror(NULL);
					}
					// printf("output coords %d %d, input coords %d %d, i = %d\n", r,c, x_r, x_c, i);

					ip = in->buf + x_r * in->stride + x_c;
					op = out->buf + r * out->stride + c;
					float sum;
					i = 5 * (n1)+n2;
					sum = 0.0F;
					for (int y = -FILTER_EXTENT; y <= FILTER_EXTENT; y++) {
						for (int x = -FILTER_EXTENT; x <= FILTER_EXTENT; x++) {
							//printf("multiply %f and %f\n", ip[y*in->stride + x], q_m[i][y*FILTER_DIM + x]);
							sum += ip[y*in->stride + x] * q_m[i][y*FILTER_DIM + x];
						}
					}
					//printf("sum is %f\n", sum);
					*op = sum;
				}
			}
#ifdef DEBUG
			printf("Done processing image\n");
#endif
		
		}
		else { // Bilinear interp 
			/*  PSF is
			*	
			*	s1*s2			(1-s2)*s1
			*
			*	s2*(1-s1)	(1-s2)*(1-s1)
			*
			*/
			float s1, s2, x1, x2;
			int n1,n2;
			float *ip, *op;
			for (r = 0; r < out->height; r++) {
				for (c = 0; c < out->width; c++) {
					x1 = r * (3.0f / 5);
					s1 = ceil(x1) - x1;
					x2 = c * (3.0f / 5);
					s2 = ceil(x2) - x2;
					n1 = (int)ceil(x1);
					n2 = (int)ceil(x2);

					ip = in->buf + n1 * in->stride + n2;
					op = out->buf + r * out->stride + c;

					*op = (1 - s2)*((1 - s1)*ip[0] + s1 * ip[-1 * in->stride]) + s2 * ((1 - s1)*ip[-1] + s1 * ip[-1 * in->stride - 1]);
				}
			}
		}


	}
	else //  Reduction
	{
		int m1, m2, n1, n2, x_r, x_c;
		float *ip, *op;
		// Perform the convolution
#ifdef DEBUG
		printf("Performing convolution: Output image %dx%d\n", out->height, out->width);
#endif
		for (r = 0; r < out->height; r++) {
			for (c = 0; c < out->width; c++) {
				m1 = r / 3; m2 = c / 3;
				n1 = r % 3; n2 = c % 3;
				switch (n1) {
				case 0:
					x_r = 5 * m1; break;
				case 1:
					x_r = 5 * m1 + 2; break;
				case 2:
					x_r = 5 * m1 + 3; break;
				default:
					perror(NULL);
				}
				switch (n2) {
				case 0:
					x_c = 5 * m2; break;
				case 1:
					x_c = 5 * m2 + 2; break;
				case 2:
					x_c = 5 * m2 + 3; break;
				default:
					perror(NULL);
				}
				// printf("output coords %d %d, input coords %d %d, i = %d\n", r,c, x_r, x_c, i);

				ip = in->buf + x_r * in->stride + x_c;
				op = out->buf + r * out->stride + c;
				float sum;
				i = 3 * (n1) + n2;
				sum = 0.0F;
				for (int y = -FILTER_EXTENT; y <= FILTER_EXTENT; y++) {
					for (int x = -FILTER_EXTENT; x <= FILTER_EXTENT; x++) {
						//printf("multiply %f and %f\n", ip[y*in->stride + x], q_m[i][y*FILTER_DIM + x]);
						sum += ip[y*in->stride + x] * q_m[i][y*FILTER_DIM + x];
					}
				}
				//printf("sum is %f\n", sum);
				*op = sum;
			}
		}
#ifdef DEBUG
		printf("Done processing image\n");
#endif

	}
}