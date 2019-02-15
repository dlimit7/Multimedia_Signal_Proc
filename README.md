# ELEC4622 Labs and Projects

          Lab 1 - Reading & Writing BMP Images
          Lab 2 - Filter Implementation & Vectorisation, Boundary extension methods
          Lab 3 - DFT & Separability
          Lab 4 - Block Motion Estimation

          Project 1 - Image Reduction & Expansion, via Sinc or Bilinear Interpolation Kernels
                    - Includes optimised methods using separable filters as well as SSE2 Intrinsics
                    - Difference Image Generator + Statistics
          
          Project 2 - Laplacian of Gaussian Filter for Edge Detection on Greyscale Images and MFLEX Videos
                    - Implemented using Floating Point Arithmetic, with separable filters
                    - Simple Morphological Dilation on foreground using a 3x3 Cross Structuring Set
          
          Project 3 - Global Motion Estimation based on Local Block Motion Vectors
                    - Generates a Motion Compensated Image.
                    - Point-Symmetric Boundary Extension and Bilinear Interpolation used
                    - Motion vectors made to minimise MSE rather than SAD over the block
                    TODO
                    - Modify the global motion estimation strategy to compute a reliability weight for each local motion block
                      based on the figure of merit that is used in the Harris corner detector. 
                      This method requires you to perform Gaussian filtering, with scale parameter σ, take local horizontal and vertical
                      differences, compute the 2 × 2 matrices Γσ [n] and aggregate these over each block to form Γ¯σ,b, 
                      after which the ratio of determinant to rank can be computed and used as the weight wb. 
                      The should then find the weighted average of the block motion vectors vb, using this as the global vector v.
