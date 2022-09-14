This program implements an optical model to calculate the irradiance pattern and
peak retinal irradiance in an aberrated eye from incident monochromatic,
coherent plane waves. The model embeds aberration data from two published
studies that was clinically measured and adapted to a range of pupil sizes. The
studies are the Indiana Aberration Study (IAS) from L.N.Thibos et.al, and the
Texas Investigation of Normal and Cataract Optics (TINCO) from R.A.Applegate
et.al. These studies tabulate aberration data in accordance with ANSI Z80.28
using Zernike coefficients.

Since, in the face of aberrations, the peak irradiance need not occur at the
focal plane. The model is configured to sweep through a small volume around the
focal point of an unaberrated eye lens. The sweep can be by shift of the image
plane or by adjustment of the 4th Zernike coefficent (focus).

This program is for demonstration purposes, takes a brute force approach and
presents a heavy computational load. It computes a mesh of phase adjusted
Huygens sources in the pupil plane and then sums the contribution of each source
at each point in each image plane comprising the swept volume. It uses the
Rayleigh-Sommerfeld Diffraction Formula with the Fresnel approximation. By
default the computation in an image plane is parallelized across 64 threads
using the pthreads library.

Each image plane irradiance result is saved in a normalized 16-bit plain PGM
format image so that the peak irradiance has a pixel value of 65535 and includes
an embedded comment recording the peak irradiance in nominal, but standardized,
units to allow comparison between different image plane irradiance results.

Image file names embed the parameters used in the creation of that image.

The name of the open source project, IOM, stands for Irradiance Optical Model.
The program consists of a single .c source file named “airy”, which is an
allusion to the [Airy disk](https://en.wikipedia.org/wiki/Airy_disk) produced in
the focal plane when imaging a plane wave through a perfect lens with a circular
aperture. The .c source file includes .inc files containing clinical aberration
measurements. It requires pthreads and the standard math library. To compile on
linux with gcc use a command line such as:

```
gcc -o airy -O5 airy.c -lm -pthread
```
