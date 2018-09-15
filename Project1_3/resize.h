#include "aligned_image_comps.h"

class my_resizer 
{
	public:
		my_resizer()
		{
			H = 0; 
			int i;
			for (i = 0; i < 5; i++) {
				g[i] = NULL;
			}
		}
		~my_resizer()
		{
			cleanup();
		}

		void init(int H, bool is_expand, bool is_sinc_interp);

		void apply_filter(my_aligned_image_comp * in, my_aligned_image_comp * out); 

	private:
		void cleanup()
		{
			
			int i;
			H = 0;
			for (i = 0; i < 5; i++) {
				g[i] = NULL;
			}
			
		}
		int H;
		bool is_expand;
		bool is_sinc_interp;
		// Discrete 2D Kernels
		//float *q[25];
		// 1D Kernels
		float *g[5];
};

