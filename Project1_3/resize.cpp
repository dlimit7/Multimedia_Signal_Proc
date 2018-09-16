#define _USE_MATH_DEFINES

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <emmintrin.h> // Include SSE2 processor intrinsic functions
#include "resize.h"
#include "aligned_image_comps.h"

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
#define FILTER_TAPS (2*FILTER_EXTENT+1) 


void my_resizer::init(int H, bool is_expand, bool is_sinc_interp)
{
	this->H = H;
	this->is_expand = is_expand;
	this->is_sinc_interp = is_sinc_interp;
	int i, r, c;
	float *g_m[5];
	__m128 *g_intr_m[5];
	for (i = 0; i < 5; i++) {
		g[i] = new float[FILTER_TAPS]; 
		g_m[i] = g[i] + FILTER_EXTENT;
#ifdef INTRINSICS
		g_intr[i] = new __m128[FILTER_TAPS];
		g_intr_m[i] = g_intr[i] + FILTER_EXTENT;
#endif
	}
	double alpha;
	if (is_expand) 
	{
		alpha = 3.0f / 5;
		if (is_sinc_interp)
		{
			// Make the windowed sinc kernels
			float bank[5] = { 0, 0.4f, -0.2f, 0.2f, -0.4f };
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
#endif
			// Normalising
			float gain;
			for (i = 0; i < 5; i++) { // For each kernel..
				gain = 0;
				// Find their DC Gain
				for (r = -FILTER_EXTENT; r <= FILTER_EXTENT; r++) {
					gain += g_m[i][r];
				}
				// Make DC Gain = 1
				gain = 1.0f / gain;
				for (r = -FILTER_EXTENT; r <= FILTER_EXTENT; r++) {
					g_m[i][r] = g_m[i][r] * gain;
				}
			}			
#ifdef INTRINSICS
			for (i = 0; i < 5; i++) {
				for (int t = -FILTER_EXTENT; t <= FILTER_EXTENT; t++)
					g_intr_m[i][t] = _mm_set1_ps(g_m[i][t]);
			}
#endif
		} else // Bilinear interp
		{
			/*  PSF is
			*	
			*	s1*s2			(1-s2)*s1
			*
			*	s2*(1-s1)	(1-s2)*(1-s1)
			*/
			// Nothing needs to be done
			
		}
	}
	else // Reduction
	{
		alpha = 5.0f / 3;
		float constant = (1.0f / alpha) ;
		// Make the windowed sinc kernels
		float bank[3] = { 0, 1.0f / 3, -1.0f / 3 };
		for (i = 0; i < 3; i++) {
			for (r = -FILTER_EXTENT; r <= FILTER_EXTENT; r++) {
				g_m[i][r] = windowed_sinc(r + bank[i % 3], H, alpha);
			}
		}
		// Normalising
		float gain;
		for (i = 0; i < 3; i++) { // For each kernel..
			gain = 0;
			// Find their DC Gain
			for (r = -FILTER_EXTENT; r <= FILTER_EXTENT; r++) {
				gain += g_m[i][r];
			}
			// Make DC Gain = 1
			gain = 1.0f / gain;
			for (r = -FILTER_EXTENT; r <= FILTER_EXTENT; r++) {
				g_m[i][r] = g_m[i][r] * gain;
			}
		}
#ifdef INTRINSICS
		for (i = 0; i < 3; i++) {
			for (int t = -FILTER_EXTENT; t <= FILTER_EXTENT; t++)
				g_intr_m[i][t] = _mm_set1_ps(g_m[i][t]);
		}
#endif
	}

}


void my_resizer::apply_filter(my_aligned_image_comp *in, my_aligned_image_comp* out) {
	int r, c, i;

	float  *g_m[5];
	__m128 *g_intr_m[5];
	for (i = 0; i < 5; i++) {
		if (g[i] != NULL) g_m[i] = g[i] + FILTER_EXTENT;
#ifdef INTRINSICS
		g_intr_m[i] = g_intr[i] + FILTER_EXTENT;
#endif
	}
	
	if (this->is_expand)
	{
		if (this->is_sinc_interp) {
			int m1, n1, x_r, x_c;
			float *ip, *mp, *op;
			my_aligned_image_comp intermediate;
			intermediate.init(in->height, out->width, in->border);
			// Horizontal filter first
#ifndef INTRINSICS
			for (r = 0; r < intermediate.height; r++) {
				for (c = 0; c < intermediate.width; c++) {
					mp = intermediate.buf + r * intermediate.stride + c;
					m1 = c / 5;
					n1 = c % 5;
					switch (n1)
					{
					case 0:
						x_c = 3 * m1; break;
					case 1:
					case 2:
						x_c = 3 * m1 + 1; break;
					case 3:
					case 4:
						x_c = 3 * m1 + 2; break;
					default:
						perror(NULL);
					}

					ip = in->buf + r * in->stride + x_c;
					float sum;
					i = n1;
					sum = 0.0F;
					for (int x = -FILTER_EXTENT; x <= FILTER_EXTENT; x++) {
						sum += ip[x] * g_m[i][x];
					}
					*mp = sum;
				}
			}
#else
			for (int r = 0; r < intermediate.height; r += 4)
				for (int c = 0; c < intermediate.width; c++)
				{
					mp = intermediate.buf + r * intermediate.stride + c;
					m1 = c / 5;
					n1 = c % 5;
					switch (n1)
					{
					case 0:
						x_c = 3 * m1; break;
					case 1:
					case 2:
						x_c = 3 * m1 + 1; break;
					case 3:
					case 4:
						x_c = 3 * m1 + 2; break;
					default:
						perror(NULL);
					}
					ip = in->buf + r * in->stride + x_c - FILTER_EXTENT;
					__m128 input;
					__m128 sum = _mm_setzero_ps();
					for (int y = -FILTER_EXTENT; y <= FILTER_EXTENT; y++) {
						input = _mm_set_ps(*(ip + 3 * in->stride), *(ip + 2 * in->stride), *(ip + in->stride),*ip);
						sum = _mm_add_ps(sum, _mm_mul_ps(g_intr_m[n1][y], input));
						ip += 1;
					}
					float *tmp = (float*)&sum;
					*mp = *tmp;
					*(mp + intermediate.stride) = *(tmp + 1);
					*(mp + 2 * intermediate.stride) = *(tmp + 2);
					*(mp + 3 * intermediate.stride) = *(tmp + 3);
				}
#endif
			intermediate.perform_boundary_extension();
			// Vertical filter now
#ifndef INTRINSICS
			for (r = 0; r < out->height; r++) {
				for (c = 0; c < out->width; c++) {
					op = out->buf + r * out->stride + c;
					m1 = r / 5;
					n1 = r % 5;
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
					mp = intermediate.buf + x_r * intermediate.stride + c;
					float sum;
					i = n1;
					sum = 0.0F;
					for (int x = -FILTER_EXTENT; x <= FILTER_EXTENT; x++) {
						sum += mp[intermediate.stride*x] * g_m[i][x];
					}
					*op = sum;
				}
			}
#else
			int vec_stride_in = intermediate.stride / 4;
			int vec_stride_out = out->stride / 4;
			int vec_width_out = (out->width + 3) / 4; // Big enough to cover the width
			__m128 *line_out = (__m128 *) out->buf;
			__m128 *line_in;
			for (int r = 0; r < out->height; r++,
				line_out += vec_stride_out) 
			{
				for (int c = 0; c < vec_width_out; c++)
				{
					m1 = r / 5;
					n1 = r % 5;
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
					line_in = (__m128 *)(intermediate.buf + x_r * intermediate.stride);
					__m128 *ip = (line_in + c) - vec_stride_in * FILTER_EXTENT;
					__m128 sum = _mm_setzero_ps();
					for (int y = -FILTER_EXTENT; y <= FILTER_EXTENT; y++) {
						sum = _mm_add_ps(sum, _mm_mul_ps(g_intr_m[n1][y], *ip));
						ip += vec_stride_in; // Next line
					}
					line_out[c] = sum;
				}
			}
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
		int m1, n1, x_r, x_c;
		float *ip, *mp, *op;
		my_aligned_image_comp intermediate;
		intermediate.init(in->height, out->width, in->border);
		// Horizontal filter first
#ifndef INTRINSICS
		for (r = 0; r < intermediate.height; r++) {
			for (c = 0; c < intermediate.width; c++) {
				mp = intermediate.buf + r * intermediate.stride + c;
				m1 = c / 3; 
				n1 = c % 3; 
				switch (n1) {
				case 0:
					x_c = 5 * m1; break;
				case 1:
					x_c = 5 * m1 + 2; break;
				case 2:
					x_c = 5 * m1 + 3; break;
				default:
					perror(NULL);
				}
				ip = in->buf + r * in->stride + x_c;
				float sum;
				i = n1;
				sum = 0.0F;
				for (int x = -FILTER_EXTENT; x <= FILTER_EXTENT; x++) {
					sum += ip[x] * g_m[i][x];
				}
				*mp = sum;
			}
		}
#else
		for (int r = 0; r < intermediate.height; r += 4)
			for (int c = 0; c < intermediate.width; c++)
			{
				mp = intermediate.buf + r * intermediate.stride + c;
				m1 = c / 3;
				n1 = c % 3;
				switch (n1) {
				case 0:
					x_c = 5 * m1; break;
				case 1:
					x_c = 5 * m1 + 2; break;
				case 2:
					x_c = 5 * m1 + 3; break;
				default:
					perror(NULL);
				}
				ip = in->buf + r * in->stride + x_c - FILTER_EXTENT;
				__m128 input;
				__m128 sum = _mm_setzero_ps();
				for (int y = -FILTER_EXTENT; y <= FILTER_EXTENT; y++) {
					input = _mm_set_ps(*(ip + 3 * in->stride), *(ip + 2 * in->stride), *(ip + in->stride), *ip);
					sum = _mm_add_ps(sum, _mm_mul_ps(g_intr_m[n1][y], input)); // Pixels in same column apply the same horizontal filter
					ip += 1;
				}
				float *tmp = (float*)&sum;
				*mp = *tmp;
				*(mp + intermediate.stride) = *(tmp + 1);
				*(mp + 2 * intermediate.stride) = *(tmp + 2);
				*(mp + 3 * intermediate.stride) = *(tmp + 3);
			}
#endif
		intermediate.perform_boundary_extension();
		// Vertical filter now
#ifndef INTRINSICS
		for (r = 0; r < out->height; r++) {
			for (c = 0; c < out->width; c++) {
				op = out->buf + r * out->stride + c;
				m1 = r / 3;
				n1 = r % 3;
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
				mp = intermediate.buf + x_r * intermediate.stride + c;
				float sum;
				i = n1;
				sum = 0.0F;
				for (int x = -FILTER_EXTENT; x <= FILTER_EXTENT; x++) {
					sum += mp[intermediate.stride*x] * g_m[i][x];
				}
				*op = sum;
			}
		}
#else 
		int vec_stride_in = intermediate.stride / 4;
		int vec_stride_out = out->stride / 4;
		int vec_width_out = (out->width + 3) / 4; // Big enough to cover the width
		__m128 *line_out = (__m128 *) out->buf;
		__m128 *line_in;
		for (int r = 0; r < out->height; r++,
			line_out += vec_stride_out)
		{
			for (int c = 0; c < vec_width_out; c++)
			{
				m1 = r / 3;
				n1 = r % 3;
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
				line_in = (__m128 *)(intermediate.buf + x_r * intermediate.stride);
				__m128 *ip = (line_in + c) - vec_stride_in * FILTER_EXTENT;
				__m128 sum = _mm_setzero_ps();
				for (int y = -FILTER_EXTENT; y <= FILTER_EXTENT; y++) {
					//printf("do u get to here%d %d %d\n", r, c, y);
					sum = _mm_add_ps(sum, _mm_mul_ps(g_intr_m[n1][y], *ip));
					ip += vec_stride_in; // Next line
				}
				line_out[c] = sum;
			}
		}
#endif
	}
}