#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "gmath.h"
#include "matrix.h"
#include "ml6.h"

/*============================================
  IMPORANT NOTE
  Ambient light is represeneted by a color value
  Point light sources are 2D arrays of doubles.
       - The fist index (LOCATION) represents the vector to the light.
       - The second index (COLOR) represents the color.
  Reflection constants (ka, kd, ks) are represened as arrays of
  doubles (red, green, blue)
  ============================================*/

// Lighting functions
color get_lighting( double * normal, double * view, color alight, double light[2][3], double * areflect, double * dreflect, double * sreflect) {
  color i;
  return i;
}

color calculate_ambient(color alight, double * areflect ) {
  color a;
  return a;
}

color calculate_diffuse(double light[2][3], double * dreflect, double * normal ) {
  color d;
  return d;
}

color calculate_specular(double light[2][3], double * sreflect, double * view, double * normal ) {

  color s;
  return s;
}

// Limit each component of c to a max of 255
void limit_color( color * c ) {
}


// Vector functions
// Normalize vetor, should modify the parameter
void normalize(double * vector) {
	double dp = dot_product(vector, vector);
	double magnitude = sqrt(dp);

	for (int i = 0; i < 3; i++) {
		vector[i] = vector[i] / magnitude;
	}
}

// Return the dot product of a . b
double dot_product(double * a, double * b) {
	return (a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2]);
}


// Calculate the surface normal for the triangle whose first
// point is located at index i in polygons
double * calculate_normal(struct matrix * polygons, int i) {
	double ** matrix = polygons -> m;
	double * norm = malloc(3 * sizeof(double));
  	double a[3];
  	double b[3];

  	a[0] = matrix[0][i + 1] - matrix[0][i];
  	a[1] = matrix[1][i + 1] - matrix[1][i];
  	a[2] = matrix[2][i + 1] - matrix[2][i];

  	b[0] = matrix[0][i + 2] - matrix[0][i];
  	b[1] = matrix[1][i + 2] - matrix[1][i];
  	b[2] = matrix[2][i + 2] - matrix[2][i];

  	norm[0] = (a[1] * b[2]) - (a[2] * b[1]);
  	norm[1] = (a[2] * b[0]) - (a[0] * b[2]);
  	norm[2] = (a[0] * b[1]) - (a[1] * b[0]);

  	return norm;
}