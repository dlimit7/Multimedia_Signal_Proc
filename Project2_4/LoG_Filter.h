#pragma once

#include <stdio.h>
#include "image_comps.h"

#define DEBUG

class my_LoG
{
public:
	my_LoG()
	{
		H = 0;
		alpha = 0;
		sigma = 0;
		g1 = NULL;
		g2 = NULL;
	}
	~my_LoG()
	{
		cleanup();
	}
	
	void init(float sigma, float alpha, int H);

	void apply_filter(my_image_comp* in, my_image_comp* out);

	void apply_dilation(my_image_comp* in, my_image_comp *out, int * A, int n);


private:
	void cleanup()
	{
		H = 0;
		alpha = 0;
		sigma = 0;
		if (g1 != NULL)
			delete[] g1;
		g1 = NULL;
		if (g2 != NULL)
			delete[] g2;
		g2 = NULL;
	}
	// Filter halflength
	int H;
	float alpha;
	float sigma;
	// Laplacian of Gaussian Filter
	float *g1;
	float *g2;
};