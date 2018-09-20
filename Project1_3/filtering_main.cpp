/*****************************************************************************/
// File: filtering_main.cpp
// Author: David Taubman
// Last Revised: 13 August, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#include "io_bmp.h"
#include "aligned_image_comps.h"
#include "resize.h"
#include <cmath>
#include <stdlib.h>


/* ========================================================================= */
/*                 Implementation of `my_image_comp' functions               */
/* ========================================================================= */

/*****************************************************************************/
/*                  my_image_comp::perform_boundary_extension                */
/*****************************************************************************/

void my_aligned_image_comp::perform_boundary_extension()
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


/*****************************************************************************/
/*                                    main                                   */
/*****************************************************************************/

int
  main(int argc, char *argv[])
{
  if (argc != 6)
    {
      fprintf(stderr,"Usage: %s <in bmp file> <out bmp file> <H> <is_expand> <is_sinc_interp>\n",argv[0]);
      return -1;
    }

  int err_code=0;
  try {
      // Read the input image
      bmp_in in;
      if ((err_code = bmp_in__open(&in,argv[1])) != 0)
        throw err_code;

      int width = in.cols, height = in.rows;
      int n, num_comps = in.num_components;
      my_aligned_image_comp *input_comps = new my_aligned_image_comp[num_comps];
      for (n=0; n < num_comps; n++)
        input_comps[n].init(height,width, 16); // Leave a border of 16
      
      int r; // Declare row index
      io_byte *line = new io_byte[width*num_comps];
      for (r=height-1; r >= 0; r--)
        { // "r" holds the true row index we are reading, since the image is
          // stored upside down in the BMP file.
          if ((err_code = bmp_in__get_line(&in,line)) != 0)
            throw err_code;
          for (n=0; n < num_comps; n++)
            {
              io_byte *src = line+n; // Points to first sample of component n
              float *dst = input_comps[n].buf + r * input_comps[n].stride;
				  for (int c = 0; c < width; c++, src += num_comps)
                dst[c] = (float) *src; // The cast to type "float" is not
                      // strictly required here, since bytes can always be
                      // converted to floats without any loss of information.
            }
        }
      bmp_in__close(&in);
		delete[] line;

      // Allocate storage for the filtered output
#ifdef DEBUG
		printf("original height %d, width %d\n", height, width);
#endif
		int H = (atoi(argv[3]));
		bool is_expand = (atoi(argv[4])) ? true : false;
		bool is_sinc_interp = (atoi(argv[5])) ? true : false;
		double temp;
		if (is_expand) {
			temp = (5.0 / 3) * height; height = ceil(temp);
			temp = (5.0 / 3) * width; width = ceil(temp);
		}
		else {
			temp = (3.0 / 5) * height; height = ceil(temp);
			temp = (3.0 / 5) * width; width = ceil(temp);
		}
#ifdef DEBUG
		printf("New height %d, new width %d\n\n", height, width);
		printf("Setting up kernels\n");
#endif
		my_resizer resizer;
		resizer.init(H, is_expand, is_sinc_interp);
      my_aligned_image_comp *output_comps = new my_aligned_image_comp[num_comps];
      for (n=0; n < num_comps; n++)
        output_comps[n].init(height,width,0); // Don't need a border for output

      // Process the image, all in floating point (easy)

#ifdef INTRINSICS
		printf("Intrinsics used\n");
#endif
#ifdef DEBUG
		printf("Process Image: num comps %d\n", num_comps);
		if (is_expand) {
			if (is_sinc_interp) {
				printf("Expansion using Sinc interpolation\n");
			}
			else {
				printf("Expansion using Bilinear Interpolation\n");
			}
		}
		else {
			printf("Reduction by Sinc Interpolation\n");
		}
#endif
      for (n=0; n < num_comps; n++)
        input_comps[n].perform_boundary_extension();
		for (n = 0; n < num_comps; n++) {
			resizer.apply_filter(input_comps + n, output_comps + n);
		}
#ifdef DEBUG
		printf("Write out the image\n");
#endif
      // Write the image back out again
      bmp_out out;
		line = new io_byte[width*num_comps];
      if ((err_code = bmp_out__open(&out,argv[2],width,height,num_comps)) != 0)
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
      delete[] input_comps;
      delete[] output_comps;
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
