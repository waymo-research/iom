// Copyright 2022 Waymo LLC
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//     http://www.apache.org/licenses/LICENSE-2.0
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// To recreate the README.inc file:
//   awk '{print "\""$0"\\n\""}' README.md > README.inc
const char *overview =
#include "README.inc"
    ;

/* ADDITIONAL NOTES:
// * Internal spatial dimensions are in meters (m) except where interfacing with
//   clincal data pupil radii which use millimeters (mm) and ANSI Z80.28 Zernike
//   coefficients which use microns (um).
// Limitations
// * Model computes relative irradiances in arbitrary units. Objective is to
//   compare irradiance patterns from different aberrations and pupil sizes.
// * Non-dispersive model, same refractive index applied to all wavelengths.
//   Since a volume around the focal plane is swept impact of a non-dispersive
//   model is spatial scaling of computed irradiance in xy plane.
// * Scattering at the retina is not modeled, no veiling corona.
// * Sweep by defocus and sweep by focal plane shift are imperfect analogs of
//   each other. They yield comparable but not identical results.
// * Model is brute force in the spatial domain. If Fourier techniques can be
//   applied it could be much much faster.
// * Model divides the computation across patches of the image plane using
//   pthreads. On a many core processor this gives a near linear speed-up in
//   the number of cores applied. A compilation option can disable use of
//   pthreads (define NO_PTHREADS).
// * Under both gcc and clang compilation option -Ofast gives attractive
//   execution time gains at the cost of slightly different numerical results.
//   This suggests that some refactor of code could yield the same gains at
//   less extreme levels of optimization.
// * Beware data integrity! Data was provided in xlsx sheets. Depending whether
//   the cells are formatted as general or scientific notation the data may
//   suffer from a loss of precision when saved as .csv files.
// * Thin lens model is used (https://en.wikipedia.org/wiki/Thin_lens).
// * Some power falls outside the image, particularly at large pupil sizes for
//   uncorrected eyes.
// * Acceptable source mesh coarseness is limited by highly aberrated lenses
//   where the difference between neighboring sample points is a significant
//   fraction of a wavelength. Some mitigation is achieved by treating sources
//   as oriented square patches. An obliquity factor of sinc(phase_delta/2)
//   where phase_delta is the phase change across the patch.
//   Note that in the absence of aberrations each source patch is oriented
//   normal to the image forming location by the thin lens construction and the
//   image forming location is 51.2 um wide compared to focal length of 17 mm.
//   An additional gradient compensation is applied to account for the distance
//   on the image forming location from the nominal focal point. Empirically the
//   error due to source mesh coarseness is less than 1% for kSDim==300.
*/

#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <pthread.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

// Clinical study data. --------------------------------------------------------

typedef struct {
  enum { kUNSPECIFIED = 0, kIAS = 1, kTINCO_best = 2, kTINCO_unaid = 3 } id;
  int patient;
  int left;
  double age;
  double lambda_nm;  // Lambda used for the physical aberration measurement.
                     // Ignored by this model.
  double pupil_rad_mm;
  double z[66];
} StudyRow;

const StudyRow ias[] = {
#include "third_party_data/IAS_dataset.inc"
};

const StudyRow app_best[] = {
#include "third_party_data/TINCO_best_dataset.inc"
};

const StudyRow app_unaid[] = {
#include "third_party_data/TINCO_unaid_dataset.inc"
};

typedef struct {
  int id;
  int n_rows;
  const StudyRow *data;
  const char *name;
  const char *long_name;
} AberrationStudy;

const AberrationStudy data_sets[] = {
    {kIAS, sizeof ias / sizeof ias[0], ias, "IAS", "Indiana Aberration Study"},
    {kTINCO_best, sizeof app_best / sizeof app_best[0], app_best, "TINCO_best",
     "Texas Investigation of Normal and Cataract Optics - best corrected"},
    {kTINCO_unaid, sizeof app_unaid / sizeof app_unaid[0], app_unaid,
     "TINCO_unaid",
     "Texas Investigation of Normal and Cataract Optics - unaided optics"},
};

// END: Clinical study data. ---------------------------------------------------

void Fatal(const char *message) {
  perror(message);
  exit(1);
}

// Zernike machinery. ----------------------------------------------------------

int ZernikeJ2n(int j) { return (int)(-0.5 + sqrt(0.25 + 2 * j)); }

int ZernikeJn2m(int j, int n) { return 2 * j - n * (n + 2); }

double ZernikeJ2Nmn(int n, int m) { return sqrt((2 * n + 2) / (1 + (m == 0))); }

int Factorial(int n) {
  int f = 1;

  while (n > 1) {
    f *= n--;
  }
  return f;
}

int Choose(int n, int k) {
  return Factorial(n) / (Factorial(k) * Factorial(n - k));
}

double IntegerPower(double x, unsigned n) {
  double r = 1.0;
  while (n--) r *= x;
  return r;
}

double ZernikeJRho2Rmn(int n, int m, double rho) {
  double r = 0;
  for (int s = 0; s <= (n - abs(m)) / 2; s++) {
    // coeff is an integer, hence can use integer division:
    // https://en.wikipedia.org/wiki/Zernike_polynomials#Other_representations
    const int coeff =
        Factorial(n - s) / (Factorial(s) * Factorial((n + abs(m)) / 2 - s) *
                            Factorial((n - abs(m)) / 2 - s));
    r += IntegerPower(rho, n - 2 * s) * ((1 - 2 * (s & 1)) * coeff);
  }

  return r;
}

double ZernikeJxy2Rmn(int n, int m, double x, double y) {
  double r = 0;
  double rsq = x * x + y * y;
  for (int s = 0; s <= (n - abs(m)) / 2; s++) {
    // coeff is an integer, hence can use integer division:
    // https://en.wikipedia.org/wiki/Zernike_polynomials#Other_representations
    const int coeff =
        Factorial(n - s) / (Factorial(s) * Factorial((n + abs(m)) / 2 - s) *
                            Factorial((n - abs(m)) / 2 - s));
    r += IntegerPower(rsq, (n - 2 * s - abs(m)) / 2) *
         ((1 - 2 * (s & 1)) * coeff);
  }

  return r;
}

// Polar Zernike (not used).
double Z(int j, double rho, double theta) {
  const int n = ZernikeJ2n(j);
  const int m = ZernikeJn2m(j, n);
  double N = ZernikeJ2Nmn(n, m);
  double R = ZernikeJRho2Rmn(n, m, rho);
  double M = (m == 0 ? 1 : (m > 0) ? cos(theta * m) : sin(theta * -m));
  return N * R * M;
}

// Cartesian Zernike.
double Zxy(int j, double x, double y) {
  const int n = ZernikeJ2n(j);
  const int m = ZernikeJn2m(j, n);
  const double N = ZernikeJ2Nmn(n, m);
  const double R = ZernikeJxy2Rmn(n, m, x, y);
  double M = 0;
  int s;
  for (int k = (m < 0), s = 0; k <= abs(m); k += 2, s++) {
    M += (1 - 2 * (s & 1)) * Choose(abs(m), k) * IntegerPower(x, abs(m) - k) *
         IntegerPower(y, k);
  }
  return N * R * M;
}

// Special case Cartesian Z4 computation (heavily used in focus sweep).
double Z4xy(double x, double y) { return sqrt(3) * (2 * (x * x + y * y) - 1); }

// END: Zernike machinery. -----------------------------------------------------

// Model parameters. -----------------------------------------------------------

const double kLambda = 555e-9;   // Default wavelength.
const double kDilatedD = 0.007;  // Dilated pupil diameter.
const double kF = 0.017;         // Nominal effective focal length of eye.
enum {
  kSDim = 300,       // Mesh points across a nominal dilated pupil.
  kFineSteps = 100,  // Fine step for mesh cells on the pupil boundary.
  kIDim = 512,       // Output image pixels in x and y.
  kZSweep = 35,      // Default volume sweep bounds (microns).
  kZStep = 5         // Default volume sweep step (microns)
};
const double kImagePitch = 1e-7;        // 100 nm.
const double kGradientBaseline = 1e-5;  // 10 um for wavefront gradient calc.

// Representation of a source in the pupil plane. ------------------------------

typedef struct {
  double o[2];  // Source origin coordinates in pupil plane.
  double w[2];  // Source patch spatial extent.
  double g[2];  // Source patch wavefront gradient referenced to spherical wave.
  double amplitude;
  double opd;  // Optical path difference induced by thin lens + aberrations.
} HuygensSource;

// Usage message and commandline argument support. -----------------------------

static const char *program_name = NULL;

void Usage(int status, const char *format, ...) {
  FILE *f = status ? stderr : stdout;
  fprintf(f, "Usage:\n");
  fprintf(f, " %s <-h|--help>\n", program_name);
  fprintf(f, " %s [-d <data-set>] [-m <subject-mask>]\n", program_name);
  fprintf(f,
          " %s [-c <coeff-mask>] [-d <data-set>] [-f <use-focus>]\n"
          " %*s [-l <lambda>] [-m <subject-mask>] [-z <a:b:c>] <tag>\n",
          program_name, (int)strlen(program_name), "");
  fprintf(f, "Where:\n");
  fprintf(f,
          " <coeff-mask> is a hex mask forcing corresponding Zernike "
          "coefficients to zero,\n"
          "    a value of \"3F\" can approximate correction of "
          "tip/tilt/defocus/astigmatism.\n");
  fprintf(f,
          " <data-set> is a known data set. Run the program without\n"
          "            arguments for a list of data sets.");
  fprintf(f,
          " <use-focus> if 0, sweep volume with varying focal plane distance\n"
          "             if 1, sweep volume with varying Z4.\n");
  fprintf(f, " <lambda> is wavelength in m (default 555e-9).\n");
  fprintf(f, " <subject-mask> select subset of subjects in data-set.\n");
  fprintf(f,
          " <a:b:c> set start:step:end of volume sweep "
          "(default %d:%d:%d).\n",
          -kZSweep, kZStep, kZSweep);
  fprintf(f,
          " <tag> is a token added to output file names. "
          "If omitted the model is not run,\n"
          "       instead a summary of the data-sets is produced.\n");
  fprintf(f, "Examples:\n");
  fprintf(f,
          " %s -d TINCO_unaid -l 555e-9 -m 1 -z 0:1:0 rc0\n"
          "    model the first subject from the TINCO_unaid data-set\n"
          "    with lambda = 555 nm, at the focal plane.\n",
          program_name);
  fprintf(f,
          " %s -d IAS -c 3F rc1\n"
          "    model all subjects from the IAS data-set with lambda = 555 nm"
          " (default),\n"
          "    forcing the first 6 Zernike coeffs to zero to approximate "
          "corrected vision\n"
          "    and sweeping the default modeled volume.\n",
          program_name);
  fprintf(f,
          " %s -d TINCO_best -l 905-9 rc2\n"
          "    model all subjects from the TINCO_best data-set with lambda = "
          "905 nm,\n"
          "    and sweeping the default modeled volume.\n",
          program_name);
  fprintf(f,
          "Warning: the last two examples produce tens of gigabytes of output\n"
          "and take several days to complete.\n");
  if (format) {
    va_list args;
    fprintf(f, "\nError: ");
    va_start(args, format);
    vfprintf(f, format, args);
    va_end(args);
    fprintf(f, "\n");
  }
  exit(status);
}

const AberrationStudy *FindDataSetByName(const char *name) {
  for (int i = 0; i < sizeof data_sets / sizeof data_sets[0]; i++) {
    if (strcmp(data_sets[i].name, name) == 0) {
      return data_sets + i;
    }
  }
  Usage(1,
        "Unknown data_set name %s\n"
        "run the program without arguments for a list of data sets",
        name);
  return NULL;  // Not reached.
}

const AberrationStudy *FindDataSetByID(int id) {
  for (int i = 0; i < sizeof data_sets / sizeof data_sets[0]; i++) {
    if (data_sets[i].id == id) return data_sets + i;
  }
  assert(0);
  return NULL;  // Not reached.
}

int SelectFromHexMaskString(int i, const char *mask) {
  if (mask == NULL) return 1;
  const int cindex = i / 4;
  const int bindex = i % 4;
  if (cindex >= strlen(mask)) return 0;
  const int mc = toupper(mask[strlen(mask) - cindex - 1]);
  // Mask is a hex string, for readability of sparse masks '.' is permitted as
  // an alias for '0'.
  if (mc != '.' && !isxdigit(mc)) {
    Usage(1, "Malformed mask character %c in %s.",
          mask[strlen(mask) - cindex - 1], mask);
  }
  const int m = mc == '.' ? 0 : isdigit(mc) ? mc - '0' : mc - ('A' - 10);
  return (m >> bindex) & 1;
}

// Core model machinery. -------------------------------------------------------

int Sign(double v) { return v >= 0 ? 1 : -1; }

double Aberration(const StudyRow *study, double x, double y) {
  double aberration = 0;
  const int num_zernikes = sizeof(study->z) / sizeof(study->z[0]);
  ;
  const double xr = x * 1e3 / study->pupil_rad_mm;
  const double yr = y * 1e3 / study->pupil_rad_mm;
  for (int j = 0; j < num_zernikes; j++) {
    if (study->z[j]) {
      aberration += study->z[j] * 1e-6 * Zxy(j, xr, yr);
    }
  }
  return aberration;
}

void SetGradient(const StudyRow *study, HuygensSource *s) {
  for (int j = 0; j < 2; j++) {
    double a[2] = {s->o[0], s->o[1]}, b[2] = {s->o[0], s->o[1]};
    a[j] += kGradientBaseline;
    b[j] -= kGradientBaseline;
    double da[2] = {Aberration(study, a[0], a[1]),
                    Aberration(study, b[0], b[1])};
    s->g[j] = (da[0] - da[1]) / (2 * kGradientBaseline);
  }
}

// Set up Huygens sources in pupil plane.
int MakeSource(double lambda, const StudyRow *study, HuygensSource *source,
               int max_n_source) {
  const double kD = study->pupil_rad_mm * 2e-3;
  const double r_pupil = study->pupil_rad_mm * 1e-3;
  const double sourceStep = (kDilatedD / kSDim);
  assert(sourceStep < kD / sqrt(2));  // Ensure some cells fall in disk.
  const double halfSourceStep = sourceStep / 2;
  const double r_pupil_sqrd = r_pupil * r_pupil;
  const double r_spherical_wave = sqrt(kF * kF + r_pupil_sqrd);
  int n = 0;
  memset(source, 0, max_n_source);
  // Set up sources at pupil.
  // Add every source rect with inner corner inside circle
  for (int nx = 0; nx <= kSDim; nx++) {
    double gridx = (nx * 2 - kSDim) * halfSourceStep;
    for (int ny = 0; ny <= kSDim; ny++) {
      double gridy = (ny * 2 - kSDim) * halfSourceStep;
      // Is inner corner inside pupil, if yes, it contributes.
      const double inx = gridx - halfSourceStep * Sign(gridx);
      const double iny = gridy - halfSourceStep * Sign(gridy);
      if (inx * inx + iny * iny <= r_pupil_sqrd) {
        HuygensSource *s = source + (n++);
        assert(n <= max_n_source);
        s->o[0] = gridx;
        s->o[1] = gridy;
        s->w[0] = sourceStep;
        s->w[1] = sourceStep;
        s->amplitude = s->w[0] * s->w[1];
        // Calculate optical path difference due to ideal thin lens. This is the
        // difference between the spherical wavefront and the point in the pupil
        // plane.
        s->opd = r_spherical_wave -
                 sqrt(s->o[0] * s->o[0] + s->o[1] * s->o[1] + kF * kF);
        // Add aberration to optical path difference of ideal lens.
        s->opd += Aberration(study, gridx, gridy);
      }
    }
  }

  // Fine adjust sources overlapping the edge of the pupil.
  // (performing adjustment of pupil edge sources as a separate step causes some
  // duplicate computation but opens the possibility of an adaptive source mesh
  // in a future version of the code).
  for (int i = 0; i < n; i++) {
    HuygensSource *s = source + i;
    double outer[2];  // Outer corner of source patch.
    for (int k = 0; k < 2; k++) {
      outer[k] = s->o[k] + s->w[k] / 2 * Sign(s->o[k]);
    }
    // Is outer corner inside pupil, if not, partial contrib.
    if (outer[0] * outer[0] + outer[1] * outer[1] >= r_pupil_sqrd) {
      int n_tot = 0;
      int n_in = 0;
      double point_sum[2] = {0, 0};
      const double sx = s->w[0] / kFineSteps;
      const double sy = s->w[1] / kFineSteps;
      for (double fx = s->o[0] - (sx / 2.0 * (kFineSteps - 1));
           fx < s->o[0] + s->w[0] / 2.0; fx += sx) {
        for (double fy = s->o[1] - (sy / 2.0 * (kFineSteps - 1));
             fy < s->o[1] + s->w[1] / 2.0; fy += sy) {
          n_tot++;
          if (fx * fx + fy * fy <= r_pupil_sqrd) {
            n_in++;
            point_sum[0] += fx;
            point_sum[1] += fy;
          }
        }
      }
      s->amplitude *= ((double)n_in) / n_tot;
      if (s->amplitude > 0) {
        for (int j = 0; j < 2; j++) {
          double point_avg = point_sum[j] / n_in;
          s->w[j] -= fabs(s->o[j] - point_avg) * 2;
          s->o[j] = point_avg;
        }
        s->opd = r_spherical_wave -
                 sqrt(s->o[0] * s->o[0] + s->o[1] * s->o[1] + kF * kF);
        s->opd += Aberration(study, s->o[0], s->o[1]);
      }
    }
    SetGradient(study, s);
  }
  return n;
}

// Structure to support multi-threaded image formation. Logically this is the
// function signature of ImagePatch().
struct PatchArgs {
  int is, ie, js, je;
  int n_sources;
  const HuygensSource *source;
  double lambda;
  double imzdelta;
  double z4;
  double r_pupil;
  double (*image)[kIDim][2];
};

// Compute the complex amplitude on a patch in the image plane delineated by the
// bounds in args. Argument type is void to support pthreads API, but the actual
// argument must always be a struct PatchArgs *.
// Bulk consumer of CPU time.
void *ImagePatch(void *args) {
  struct PatchArgs *p = (struct PatchArgs *)args;
  int n_sources = p->n_sources;
  const HuygensSource *source = p->source;
  // Patch iteration bounds.
  int is = p->is, ie = p->ie, js = p->js, je = p->je;
  double(*image)[kIDim][2] = p->image;
  const double pi2_div_lambda = 2 * M_PI / p->lambda;
  const double z = kF + p->imzdelta;
  const double z_sqrd = z * z;
  const double z_div_lambda = z / p->lambda;
  for (int n = 0; n < n_sources; n++) {
    const HuygensSource *s = source + n;
    const double opd =
        s->opd         // Thin lens optical path difference.
        + (!p->z4 ? 0  // Option of adjustment by focus (4th Zernike).
                  : p->z4 * 1e-6 *
                        Z4xy(s->o[0] / p->r_pupil, s->o[1] / p->r_pupil));
    // x_sub_i,y_sub_i are the spatial coordinates in meters of a point in the
    // image plane.
    for (int i = is; i < ie; i++) {
      const double x_sub_i = (2 * i - kIDim + 1) * kImagePitch / 2;
      const double xd = s->o[0] - x_sub_i;
      const double x_sqrd = xd * xd;
      for (int j = js; j < je; j++) {
        const double y_sub_i = (2 * j - kIDim + 1) * kImagePitch / 2;
        const double yd = s->o[1] - y_sub_i;
        const double r_sqrd = x_sqrd + yd * yd + z_sqrd;
        const double opl_m = sqrt(r_sqrd) + opd;  // Optical path length (m).
        const double opl_pi2_lambda = opl_m * pi2_div_lambda;
        // scale is -j*z0/(lambda*r*r)  Milster eqn 5.33.
        const double scale = z_div_lambda / r_sqrd;
        const double ampC = sin(opl_pi2_lambda) * scale;
        const double ampS = cos(opl_pi2_lambda) * scale;
        double amplitude = s->amplitude;
        // Approximate separable compensation for non-normal source patch.
        // Attenuate by sinc of half the phase difference of the patch limits,
        // in x and y dimensions.
        for (int k = 0; k < 2; k++) {
          double gd = (!k ? x_sub_i : y_sub_i) / kF;
          // For small angles gradients are additive (small angle approximation)
          double phase = 2 * M_PI * (s->g[k] - gd) * s->w[k] / p->lambda;
          amplitude *= fabs(phase) < 1e-10 ? 1 : sin(phase / 2) / (phase / 2);
        }
        // Add to complex amplitude.
        image[i][j][0] += ampS * amplitude;
        image[i][j][1] += ampC * amplitude;
      }
    }
  }
  return NULL;
}

// Calculate irradiance of light field in image plane modulated by imzdelta
// which shifts the image plane in z relative to the nominal focal plane and/or
// z4 which adjusts the 4th Zernike coefficient to apply a parabolic focus
// adjustment. The complex amplitude is first calculated by ImagePatch, then
// converted to irradiance and phase, and finally save to a file.
void FormImage(double lambda, const StudyRow *study, FILE *f,
               const HuygensSource *source, int n, double imzdelta, double z4) {
  typedef double ImageType[kIDim][kIDim][2];
  // "image" stores complex amplitude during accumulation. This is then
  // converted in place to irradiance and phase.
  ImageType *image = malloc(sizeof *image);
  if (image == NULL) Fatal("malloc(image)");
  const int kPixelRange = 65535;
  memset(image, 0, sizeof *image);

  struct PatchArgs patch_args_base = {.n_sources = n,
                                      .imzdelta = imzdelta,
                                      .z4 = z4,
                                      .r_pupil = study->pupil_rad_mm * 1e-3,
                                      .source = source,
                                      .is = 0,
                                      .ie = kIDim,
                                      .js = 0,
                                      .je = kIDim,
                                      .lambda = lambda,
                                      .image = *image};

#ifdef NO_PTHREADS
  ImagePatch(&patch_args_base);
#else
  /* Divide the work among worker threads. Each thread processes a
   * contiguous slice in the x dimension.
   */
  const int kMaxWorkers = 64;
  pthread_t pt[kMaxWorkers];
  struct PatchArgs patch_args[kMaxWorkers];
  const int workers = getenv("WORKERS") ? atoi(getenv("WORKERS")) : kMaxWorkers;
  assert(workers > 0 && workers <= kMaxWorkers);
  for (int t = 0; t < workers; t++) {
    patch_args[t] = patch_args_base;
    patch_args[t].is = (t) * (kIDim / workers);
    // Handle case that workers is not an integer divisor of kIDim.
    patch_args[t].ie = (t == workers - 1) ? kIDim : (t + 1) * (kIDim / workers);
    if (pthread_create(pt + t, NULL, &ImagePatch, patch_args + t)) {
      Fatal("pthread_create");
    }
  }
  for (int t = 0; t < workers; t++) {
    if (pthread_join(pt[t], NULL)) {
      Fatal("pthread_join");
    }
  }
#endif

  // Compute power and phase from complex amplitude.
  double maxi = (*image)[0][0][0] * (*image)[0][0][0] +
                (*image)[0][0][1] * (*image)[0][0][1];
  double mini = maxi;
  double total = 0;
  for (int i = 0; i < kIDim; i++) {
    for (int j = 0; j < kIDim; j++) {
      const double arg = atan2((*image)[i][j][1], (*image)[i][j][0]);
      (*image)[i][j][0] = (*image)[i][j][0] * (*image)[i][j][0] +
                          (*image)[i][j][1] * (*image)[i][j][1];
      // Preserve arg. Not currently used.
      (*image)[i][j][1] = arg;
      total += (*image)[i][j][0];
      if ((*image)[i][j][0] < mini) mini = (*image)[i][j][0];
      if ((*image)[i][j][0] > maxi) maxi = (*image)[i][j][0];
    }
  }

  // Save result.
  fprintf(f, "P%d\n%d %d\n%d\n", 2, kIDim, kIDim, kPixelRange);
  fprintf(f, "# lambda %g imZ %g Z4 %g\n", lambda, imzdelta, z4);
  // divide by (*image) pixel area for irradiance.
  const double kPixelArea = kImagePitch * kImagePitch;
  fprintf(f, "#m %g %d\n", maxi / kPixelArea, n);  // Normalize to irradiance.
  fprintf(f, "#t %g\n", total);  // Total power captured by model. Used to
                                 // gauge how much power falls outside image.
  double range;
  if (mini < 0) {
    fprintf(stderr, "Anomalous minimum irradiance! %g\n", mini);
  }
  range = maxi;
  if (range == 0) range = 1;
  for (int i = 0; i < kIDim; i++) {
    for (int j = 0; j < kIDim; j++) {
      int pixel = lrint(kPixelRange * (*image)[i][j][0] / range);
      if (fprintf(f, "%d\n", pixel) <= 0) {
        Fatal("fprintf");
      }
    }
  }
  free(image);
}

FILE *OpenCanonicallyNamedResultFile(const char *rev, double lambda,
                                     const StudyRow *study,
                                     int sweep_by_defocus, int imz,
                                     const char *coeff_mask) {
  FILE *f;
  char fname[4096];
  const int imz_bias = 1 << 11;  // Render imz as 3 hex digits biased by 0x800.
  const char *study_name = FindDataSetByID(study->id)->name;
  int len = snprintf(
      fname, sizeof fname, "%s_%03d_%s_%1.1f_%03d_O%c_%02d%s%s_%03X%c.pgm",
      study_name, (int)(lambda * 1e9), rev, study->pupil_rad_mm * 2,
      study->patient, "DS"[study->left], (int)study -> age,
      (coeff_mask ? "_c" : ""), (coeff_mask ? coeff_mask : ""), imz + imz_bias,
      sweep_by_defocus ? 'f' : 'z');
  if (len < 0 || len >= sizeof fname) Fatal("snprintf");
  f = fopen(fname, "w");
  if (f == NULL) Fatal("fopen");
  return f;
}

// Create pupil plane sources and run model over specified volume.
void SweepVolume(const char *rev, double lambda, const StudyRow *study,
                 int z_min_um, int z_step, int z_max_um, int sweep_by_defocus,
                 const char *coeff_mask) {
  // Focus scale is empirically determined and approximate.
  const double focus_scale = -(study->pupil_rad_mm * study->pupil_rad_mm) / 2;
  // imz is limited to integer micron values. This limitation simplifies
  // consistent repeatable result file naming.

  const int max_n_source = (kSDim + 1) * (kSDim + 1);
  HuygensSource *source = malloc(max_n_source * sizeof *source);
  if (source == NULL) Fatal("malloc(source)");
  int n = MakeSource(lambda, study, source, max_n_source);
  for (int imz = z_min_um; imz <= z_max_um; imz += z_step) {
    FILE *f = OpenCanonicallyNamedResultFile(rev, lambda, study,
                                             sweep_by_defocus, imz, coeff_mask);
    FormImage(lambda, study, f, source, n, sweep_by_defocus ? 0 : 1e-6 * imz,
              sweep_by_defocus ? focus_scale * 0.001 * imz : 0);
    fclose(f);
  }
  free(source);
}

void TincoIntegrityCheck() {
  /*
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2083284/
  table 7 contains summary statistics of RMS high order wavefront errors for
  Z6 to Z14.
  The following table permit the cross check of table 7 to the raw TINCO data.
  */
  static struct ZernikeToHighOrderWavefrontErrorMapping {
    int z_coeff[2];
    const char *name;
  } z_to_ho_wfe_rms_map[] = {
      {{6, 9}, "Trefoil"},
      {{7, 8}, "Coma"},
      {{10, 14}, "Tetrafoil"},
      {{11, 13}, "Secondary Astigmatism"},
      {{12, 0}, "Spherical"}  // Other RMS results combine 2 Zernikes.
                              // Spherical just needs z12, exploit that z0=0.
  };

  enum { kNumStats = 5 };
  static struct ApplegateTable7Row {
    int age_lo, age_hi;
    double pupil_mm;
    // HO WFE (um): {z6+z9},{z7+z8},{z10+z14},{z11+z13},{z12} X {mean,SD}
    double ho_wfe[kNumStats][2];
  } summary_data[] = {
#include "third_party_data/TINCO_table7.inc"
  };

  const AberrationStudy *tinco_best = FindDataSetByID(kTINCO_best);
  const AberrationStudy *tinco_unaid = FindDataSetByID(kTINCO_unaid);

  // Check that best and unaid match for Zernikes above z5.
  assert(tinco_best->n_rows == tinco_unaid->n_rows);
  for (int i = 0; i < tinco_best->n_rows; i++) {
    const StudyRow *row_best = tinco_best->data + i;
    const StudyRow *row_unaid = tinco_unaid->data + i;
    for (int j = 6; j < sizeof row_best->z / sizeof row_best->z[0]; j++) {
      assert(row_best->z[j] == row_unaid->z[j]);
    }
  }

  // Cross check published Applegate table 7 summary data against raw data.
  for (int k = 0; k < sizeof summary_data / sizeof summary_data[0]; k++) {
    const struct ApplegateTable7Row *s = summary_data + k;
    struct statistics {
      double sumsq;
      double sum;
      int count;
    } stats[5];
    memset(stats, 0, sizeof stats);
    // Accumulate statistics.
    for (int i = 0; i < tinco_best->n_rows; i++) {
      const StudyRow *study = tinco_best->data + i;
      if (study->pupil_rad_mm == s->pupil_mm / 2 && study->age >= s->age_lo &&
          study->age < s->age_hi) {
        for (int j = 0; j < kNumStats; j++) {
          const double r0 = study->z[z_to_ho_wfe_rms_map[j].z_coeff[0]];
          const double r1 = study->z[z_to_ho_wfe_rms_map[j].z_coeff[1]];
          const double mean_square = r0 * r0 + r1 * r1;
          stats[j].sumsq += mean_square;
          stats[j].sum += sqrt(mean_square);  // RMS.
          stats[j].count += 1;
        }
      }
    }
    // Cross check statistics against Applegate table 7.
    for (int j = 0; j < kNumStats; j++) {
      const double *summary = s->ho_wfe[j];
      if (stats[j].count == 0) {
        fprintf(stderr, "TincoIntegrityCheck: count 0 %d %g %s\n", s->age_lo,
                s->pupil_mm, z_to_ho_wfe_rms_map[j].name);
        break;
      }
      const double mean = stats[j].sum / stats[j].count;
      if (nearbyint(mean * 1000) != summary[0] * 1000) {
        fprintf(stderr, "TincoIntegrityCheck: bad mean %d %g %s %g %g\n",
                s->age_lo, s->pupil_mm, z_to_ho_wfe_rms_map[j].name,
                nearbyint(mean * 1000), summary[0] * 1000);
      }
      if (stats[j].count == 1) return;
      double std = sqrt((double)stats[j].count / (stats[j].count - 1)) *
                   sqrt(stats[j].sumsq / stats[j].count - mean * mean);
      if (nearbyint(std * 1000) != summary[1] * 1000) {
        fprintf(stderr, "TincoIntegrityCheck: bad std %d %g %s %g %g %g %g\n",
                s->age_lo, s->pupil_mm, z_to_ho_wfe_rms_map[j].name, std,
                summary[1], nearbyint(std * 1000), summary[1] * 1000);
      }
    }
  }
}

int main(int argc, const char **argv) {
  const char *study_mask = NULL;
  const char *coeff_mask = NULL;
  const AberrationStudy *data_set = NULL;
  double lambda = kLambda;
  int zmin = -kZSweep, zstep = kZStep, zmax = kZSweep;
  int sweep_by_defocus = 0;
  program_name = argv[0];
  argc--;
  argv++;
  for (int i = 0; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
      fprintf(stderr, "%s\n", overview);
      Usage(0, NULL);
    }
  }
  while (argc > 1) {
    if (argc >= 2 && strcmp(argv[0], "-c") == 0) {
      coeff_mask = argv[1];
    } else if (argc >= 2 && strcmp(argv[0], "-d") == 0) {
      data_set = FindDataSetByName(argv[1]);
    } else if (argc >= 2 && strcmp(argv[0], "-f") == 0) {
      sweep_by_defocus = atoi(argv[1]) != 0;
    } else if (argc >= 2 && strcmp(argv[0], "-l") == 0) {
      lambda = atof(argv[1]);
      if (lambda < 1e-9 || lambda > 1e-3) {
        Usage(1, "Lambda = %g m, require, 1 nm <= lambda <= 1 mm", lambda);
      }
    } else if (argc >= 2 && strcmp(argv[0], "-m") == 0) {
      study_mask = argv[1];
    } else if (argc >= 2 && strcmp(argv[0], "-z") == 0) {
      if (sscanf(argv[1], "%d:%d:%d", &zmin, &zstep, &zmax) != 3 ||
          zstep <= 0 || zmin > zmax) {
        Usage(1, "-z a:b:c where a,b,c integers, for a to c step b");
      }
    } else {
      Usage(1, "Unknown flag %s", argv[0]);
    }
    argc -= 2;
    argv += 2;
  }
  const char *rev = argc > 0 ? argv[0] : "";
  if (strpbrk(rev, "/ \t\n") != NULL) {
    Usage(1, "The tag string %s contains characters best avoided in file names",
          rev);
  }
  if (argc == 0) {  // Print summary of available data sets.
    TincoIntegrityCheck();
    if (data_set == NULL) {
      printf("Available data sets:\n");
      for (int i = 0; i < sizeof data_sets / sizeof data_sets[0]; i++) {
        printf(" Data set \"%s\" id=%d subjects=%d\n", data_sets[i].name,
               data_sets[i].id, data_sets[i].n_rows);
      }
      printf(
          "Select a data-set with the -d flag for "
          "details on that data-set.\n");
      printf("Type:\n %s --help\nfor help\n", program_name);
    } else {
      if (coeff_mask) {
        Usage(1, "coeff-mask %s not yet supported in this context", coeff_mask);
      }
      printf(
          "#data-set,patient,OS/OD,age,ref-lambda,pupil-radius-mm,"
          "RMSLOA,RMSHOA,RMSTot\n");
      for (int j = 0; j < data_set->n_rows; j++) {
        if (SelectFromHexMaskString(j, study_mask)) {
          const StudyRow *s = data_set->data + j;
          double wfe_sq[2] = {0, 0};
          for (int k = 0; k < sizeof s->z / sizeof s->z[0]; k++) {
            // Collect low order aberrations in [0] and high order in [1].
            wfe_sq[k >= 6] += s->z[k] * s->z[k];
          }
          printf("%s, %d, O%c, %g, %g, %g, ... LOA %g HOA %g Tot %g\n",
                 FindDataSetByID(s->id)->name, s->patient, "DS"[s->left],
                 s -> age, s -> lambda_nm, s -> pupil_rad_mm, sqrt(wfe_sq[0]),
                 sqrt(wfe_sq[1]), sqrt(wfe_sq[0] + wfe_sq[1]));
        }
      }
    }
  } else {  // Calculate irradiance maps with selected parameters.
    if (data_set == NULL) Usage(1, "No data_set specified, use -d option.");
    for (int i = 0; i < data_set->n_rows; i++) {
      if (SelectFromHexMaskString(i, study_mask)) {
        StudyRow study = data_set->data[i];
        if (coeff_mask) {
          for (int j = 0; j < sizeof study.z / sizeof study.z[0]; j++) {
            if (SelectFromHexMaskString(j, coeff_mask)) {
              memset(study.z + j, 0, sizeof study.z[j]);
            }
          }
        }
        SweepVolume(rev, lambda, &study, zmin, zstep, zmax, sweep_by_defocus,
                    coeff_mask);
      }
    }
  }
  exit(0);
}
