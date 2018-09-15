/*****************************************************************************/
// File: filtering_main.cpp
// Author: David Taubman
// Last Revised: 13 August, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#include "io_bmp.h"
#include "image_comps.h"
#include "resize.h"
#include <cmath>
#include <stdlib.h>


/* ========================================================================= */
/*                 Implementation of `my_image_comp' functions               */
/* ========================================================================= */

/*****************************************************************************/
/*                  my_image_comp::perform_boundary_extension                */
/*****************************************************************************/

void my_image_comp::perform_boundary_extension()
{
	/* Z-order Hold Boundary Extension Method*/
	int r, c;

	// First extend upwards
	float *first_line = buf;
	for (r=1; r <= border; r++)
		for (c=0; c < width; c++)
			first_line[-r*stride+c] = first_line[c];

	// Now extend downwards
	float *last_line = buf+(height-1)*stride;
	for (r=1; r <= border; r++)
		for (c=0; c < width; c++)
			last_line[r*stride+c] = last_line[c];

	// Now extend all rows to the left and to the right
	float *left_edge = buf-border*stride;
	float *right_edge = left_edge + width - 1;
	for (r=height+2*border; r > 0; r--, left_edge+=stride, right_edge+=stride)
		for (c=1; c <= border; c++)
		{
			left_edge[-c] = left_edge[0];
			right_edge[c] = right_edge[0];
		}
		

	/*
		// Zero Padding Method
		  int r, c;

	  // First extend upwards
	  float *first_line = buf;
	  for (r=1; r <= border; r++)
		  for (c=0; c < width; c++)
			  first_line[-r*stride+c] = 0;

	  // Now extend downwards
	  float *last_line = buf+(height-1)*stride;
	  for (r=1; r <= border; r++)
		  for (c=0; c < width; c++)
			  last_line[r*stride+c] = 0;

	  // Now extend all rows to the left and to the right
	  float *left_edge = buf-border*stride;
	  float *right_edge = left_edge + width - 1;
	  for (r = height + 2 * border; r > 0; r--, left_edge += stride, right_edge += stride) {
		  for (c = 1; c <= border; c++)
		  {
			  left_edge[-c] = 0;
			  right_edge[c] = 0;
		  }
	  }
	*/  
	/*
	  // Symmetric Boundary Extension Method 
	int r, c;

	// Top and Bottom, back and forth to account for if the border is larger than the image itself 
	float *first_line = buf;
	float *last_line = buf + (height - 1)*stride;
	for (r = 1; r <= border; r++) {
		for (c = 0; c < width; c++) {
			first_line[-r * stride + c] = first_line[(r)*stride + c]; //  Y X O | X Y 
			last_line[r * stride + c] = last_line[-(r)*stride + c];
		}
	}
	// Left and Right
	float *left_edge = handle + border;
	float *right_edge = left_edge + (width - 1);
	for (r = height + 2 * border; r > 0; r--, left_edge += stride, right_edge += stride) {
		for (c = 1; c <= border; c++) {
			left_edge[-c] = left_edge[c];
			right_edge[c] = right_edge[-c];
		}
	}
	/*
	printf("Left edge boundary check\n");
	left_edge = buf;
	right_edge = left_edge + width - 1;
	for (r = -border; r <= border; r++) {
		for (c = -border; c <= border; c++) {
			if (c == 0) {
				printf("| ");
			}
			printf("%lf  ", left_edge[r*stride + c]);
		}
		printf("\n");
	}
	
	printf("\n\nRight edge now\n");
	for (r = -border; r <= border; r++) {
		for (c = -border; c <= border; c++) {
			if (c == 0) {
				printf("| ");
			}
			printf("%lf  ", right_edge[r*stride + c]);
		}
		printf("\n");
	}
	printf("\n\n");
	*/
}

/* ========================================================================= */
/*                              Global Functions                             */
/* ========================================================================= */

/*****************************************************************************/
/*                                apply_filter                               */
/*****************************************************************************/



float filter_coeff(float *arr, int n) {
	int i;
	//printf("size of array in %d\n", n);
	float sum = 0.0;
	for (i = 0; i < n; i++) {
		sum += arr[i];
	}
	return 1.0F/sum;
}
/*
void apply_filter(my_image_comp *in, my_image_comp *out, bool is_expand, bool is_sinc_interp)
{

	printf("Setting up kernels\n");
	my_resizer resizer;
	resizer.init(5, 1, 1);
	resizer.apply_filter(in, out);

}*/

/*****************************************************************************/
/*                                    main                                   */
/*****************************************************************************/

int
  main(int argc, char *argv[])
{
  if (argc != 4)
    {
      fprintf(stderr,"Usage: %s <in bmp file> <in bmp file> <out bmp file> \n",argv[0]);
      return -1;
    }

  int err_code=0;
  try {
      // Read the input image #1 
      bmp_in in;
      if ((err_code = bmp_in__open(&in,argv[1])) != 0)
        throw err_code;

      int width1 = in.cols, height1 = in.rows;
      int n, num_comps1 = in.num_components;
      my_image_comp *input1_comps = new my_image_comp[num_comps1];
      for (n=0; n < num_comps1; n++)
        input1_comps[n].init(height1,width1, 0); // Leave a border of 4
      
      int r; // Declare row index
      io_byte *line = new io_byte[width1*num_comps1];
      for (r=height1-1; r >= 0; r--)
        { // "r" holds the true row index we are reading, since the image is
          // stored upside down in the BMP file.
          if ((err_code = bmp_in__get_line(&in,line)) != 0)
            throw err_code;
          for (n=0; n < num_comps1; n++)
            {
              io_byte *src = line+n; // Points to first sample of component n
              float *dst = input1_comps[n].buf + r * input1_comps[n].stride;
				  for (int c = 0; c < width1; c++, src += num_comps1)
                dst[c] = (float) *src; // The cast to type "float" is not
                      // strictly required here, since bytes can always be
                      // converted to floats without any loss of information.
            }
        }
      bmp_in__close(&in);
		delete[] line;

		// Read the input image #2
		if ((err_code = bmp_in__open(&in, argv[2])) != 0)
			throw err_code;

		int width2 = in.cols, height2 = in.rows;
		int num_comps2 = in.num_components;
		my_image_comp *input2_comps = new my_image_comp[num_comps2];
		for (n = 0; n < num_comps2; n++)
			input2_comps[n].init(height2, width2, 0); // Leave a border of 0

		line = new io_byte[width2*num_comps2];
		for (r = height2 - 1; r >= 0; r--)
		{ // "r" holds the true row index we are reading, since the image is
		  // stored upside down in the BMP file.
			if ((err_code = bmp_in__get_line(&in, line)) != 0)
				throw err_code;
			for (n = 0; n < num_comps2; n++)
			{
				io_byte *src = line + n; // Points to first sample of component n
				float *dst = input2_comps[n].buf + r * input2_comps[n].stride;
				for (int c = 0; c < width2; c++, src += num_comps2)
					dst[c] = (float)*src; // The cast to type "float" is not
												 // strictly required here, since bytes can always be
												 // converted to floats without any loss of information.
			}
		}
		bmp_in__close(&in);
		delete[] line;

      // Allocate storage for the difference image
		int height, width, num_comps;
		height = (height1 > height2) ? height2 : height1; // get the smaller dimensions
		width = (width1 > width2) ? width2 : width1;
		num_comps = (num_comps1 > num_comps2) ? num_comps2 : num_comps1;

		printf("New height %d, new width %d\n\n", height, width);

      my_image_comp *output_comps = new my_image_comp[num_comps];
      for (n=0; n < num_comps; n++)
        output_comps[n].init(height,width,0); // Don't need a border for output

      // Process the input images --> generate difference images, Mean Error, MSE, PSNR per component
		printf("Process Image: num comps %d\n", num_comps);
		float *meanError = new float[num_comps];
		float *MSE = new float[num_comps];
		float* PSNR = new float[num_comps];
		for (n = 0; n < num_comps; n++) {
			meanError[n] = 0.0f;
			MSE[n] = 0.0f;
			PSNR[n] = 0.0f;
		}
		float diff, constant = 1.0f / (width*height);
		int row, col;
		for (n = 0; n < num_comps; n++) {
			for (row = 0; row < height; row++) {
				float * x = input1_comps[n].buf + row * input1_comps[n].stride;
				float * y = input2_comps[n].buf + row * input2_comps[n].stride;
				float * out = output_comps[n].buf + row * output_comps[n].stride;
				for (col = 0; col < width; col++) {
					diff = x[col] - y[col];
					meanError[n] += diff;
					MSE[n] += diff * diff;
					out[col] = 128.0f + (0.5f)*diff;
				}
			}
			meanError[n] = meanError[n] * constant;
			MSE[n] = MSE[n] * constant;
			PSNR[n] = 10 * log10(255 * 255 * (1 / MSE[n]));
		}
		for (n = 0; n < num_comps; n++) {
			printf("For image component #%d,\n Mean Error is %f\n MSE is %f\n PSNR is %f\n", n, meanError[n], MSE[n], PSNR[n]);
		}

		printf("Write out the image\n");
      // Write the image back out again
      bmp_out out;
		line = new io_byte[width*num_comps];
      if ((err_code = bmp_out__open(&out,argv[4],width,height,num_comps)) != 0)
        throw err_code;
      for (r=height-1; r >= 0; r--)
		{ // "r" holds the true row index we are writing, since the image is
				// written upside down in BMP files.
			//printf("row %d\n", r);
         for (n=0; n < num_comps; n++)
         {
            io_byte *dst = line+n; // Points to first sample of component n
            float *src = output_comps[n].buf + r * output_comps[n].stride;
				for (int c = 0; c < width; c++, dst += num_comps) {
					float t = floor(src[c]+0.5F);
					t = (t > 255.0f) ? 255.0f : t;
					t = (t < 0.0f) ? 0.0f : t;
					*dst = (io_byte)t; // The cast to type "io_byte" is
							// required here, since floats cannot generally be
							// converted to bytes without loss of information.  The
							// compiler will warn you of this if you remove the cast.
							// There is in fact not the best way to do the
							// conversion.  You should fix it up in the lab.
				}
         }
         bmp_out__put_line(&out,line);
		}
      bmp_out__close(&out);
      delete[] line;
      delete[] input1_comps;
		delete[] input2_comps;
      delete[] output_comps;
		printf("DONNNEEEEE\n");
  }
  catch (int exc) {
      if (exc == IO_ERR_NO_FILE)
        fprintf(stderr,"Cannot open supplied input or output file.\n");
      else if (exc == IO_ERR_FILE_HEADER)
        fprintf(stderr,"Error encountered while parsing BMP file header.\n");
      else if (exc == IO_ERR_UNSUPPORTED)
        fprintf(stderr,"Input uses an unsupported BMP file format.\n  Current "
                "simple example supports only 8-bit and 24-bit data.\n");
      else if (exc == IO_ERR_FILE_TRUNC)
        fprintf(stderr,"Input or output file truncated unexpectedly.\n");
      else if (exc == IO_ERR_FILE_NOT_OPEN)
        fprintf(stderr,"Trying to access a file which is not open!(?)\n");
      return -1;
  }
  return 0;
}
