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
		g = NULL;
	}
	~my_LoG()
	{
		cleanup();
	}
	
	void init(float sigma, float alpha, int H);

	void apply_filter(my_image_comp* in, my_image_comp* out);


private:
	void cleanup()
	{
		H = 0;
		alpha = 0;
		sigma = 0;
		if (g != NULL)
			delete[] g;
		g = NULL;
	}
	// Filter halflength
	int H;
	float alpha;
	float sigma;
	// Laplacian of Gaussian Filter
	float *g;
};