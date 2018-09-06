#include "image_comps.h"

class my_resizer 
{
	public:
		my_resizer()
		{
			H = 0; 
			int i;
			for (i = 0; i < 25; i++) {
				q[i] = NULL;
			}
		}
		~my_resizer()
		{
			cleanup();
		}
		/*	~~~~~~~ Init ~~~~~~~
		*
		*	Initialise the discrete kernels required for expansion or reduction
		*	Normalises the filters to DC Gain 1 and also precalculates and stores their mirror_psfs
		*  Results are stored in the 2D array q. There will be a fixed max of 16 kernels required in this design.
		*
		*	H ---->
		*
		*	is_expand --->
		*
		*	is_sinc_interp --->
		*
		*
		*/
		void init(int H, bool is_expand, bool is_sinc_interp);


		/*	~~~~~~~ Apply Filter ~~~~~~~
		*
		*	~~~ Sinc Interpolation ~~~
		*	~~ Expansion ~~
		*  y[5m+n] can be broken down into:
		*	y[5m1+n1, 5m2+n2]
		*		 r	  ,	c
		*	n1 = r%5, n2 = c%5
		*
		*  if 0 < n1 < 2, pixel of interest is x[3m1+1]
		*  else,      pixel of interest is x[3m1+2]
		*
		*	In the special case when either n1 or n2 = 0,
		*	the output pixel y[5m] = x[3m]
		*
		*	Kernels q_k where -k = [0, 0.4, -0.2, 0.2, -0.4]
		*	
		*	i		Kernels			n1,n2
		*	0		0.4,  0.4		1, 1
		*	1		0.4, -0.2		1, 2
		*	2		0.4,  0.2		1, 3
		*	3		0.4, -0.4		1, 4
		*	4	  -0.2,  0.4		2, 1
		*	5	  -0.2, -0.2		2, 2
		*	6			...			 ...
		*
		*	It should be clear that the kernel index i = 4*(n1-1) + (n2 - 1)
		*
		*
		*	~~ Reduction ~~
		*	y[3m1+n1, 3m2+n2]
		*		 r	  ,	c
		*	n1 = r%3, n2 = c%3
		*
		*	if n == 0, poi is x[5m]
		*	if n == 1, poi is x[5m+2]
		*	if n == 2, poi is x[5m+3]
		*(
		*	Kernels q_k where -k = [0, 1/3, -1/3]
		*
		*	i		Kernels			n1, n2
		*	0		0,  0				0, 0
		*	1		0,	 1/3			0, 1
		*	2		0, -1/3			0, 2
		*	3		1/3, 0			1, 0
		*	4		etc 
		*	5
		*
		*	It should be clear the the kernel index is i = 3*(n1) + n2
		*
		*/
		void apply_filter(my_image_comp * in, my_image_comp * out);

	private:
		void cleanup()
		{
			
			int i;
			if (q != NULL) {
				for (i = 0; i < 25; i++) {
					if (q[i] != NULL) {
						delete[] q[i];
					}
				}
			}
			H = 0;
			for (i = 0; i < 25; i++) {
				q[i] = NULL;
			}
			
		}
		int H;
		bool is_expand;
		bool is_sinc_interp;
		// Discrete 2D Kernels
		float *q[25];
		// 1D Kernels
		float *g[5];
};

