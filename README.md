# 2d-wavelets

This project implements the 2-dimensional discrete wavelet decomposition of an image, using the 1-dimensional discrete Daubechies wavelet.

1) An original NxN image is high-pass filtered and gives:

i) an N/2 x N/2 image called *approximation* image, which is the result of a _low_-pass filtering of the rows and a subsequent _low_-pass filtering of the columns

ii) an N/2 x N/2 image called *horizontal detail*, which is the result of a _low_-pass filtering of the rows and a subsequent _high_-pass filtering of the columns

iii) an N/2 x N/2 image called *vertical detail*, which is the result of a _high_-pass filtering of the rows and a subsequent _low_-pass filtering of the columns

iv) an N/2 x N/2 image called *diagonal detail*, which is the result of a _high_-pass filtering of the rows and a subsequent _high_-pass filtering of the columns

2) One can apply the above described procedure again, on the approximation image. This would yield four N/4 x N/4 images. In general, one can apply this procedure *l* times, which would yield a N/(2^l) x N/(2^l) approximation image, where *l* the level of the wavelet decomposition.

3) The low-pass filter is implemented via an *o*th-order Daubechies wavelet, called the _mother_ wavelet. The high-pass filter is the quadratic mirror filter of the mother wavelet, and is called the _father_wavelet

In order to compile the program, run on a command window that supports CUDA:
nvcc 2d-daubechies.cu -o 2d-daubechies.out

which will create the executable file. Then run
./2d-daubechies.out 4 1024 2

which will tell the program to decompose both serially and parallely a random 1024x1024 image using a 4th order Daubechies wavelet, and decompose the image in 2 levels, yielding a 256x256 approximation image.
