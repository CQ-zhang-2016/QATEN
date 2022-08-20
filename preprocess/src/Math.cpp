#include <cmath>
#include "Math.h"

namespace Math {

double dot(double x[], double y[]) {

	return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];

}

void cross(double x[], double y[], double res[]) {

	res[0] = x[1] * y[2] - x[2] * y[1];
	res[1] = x[2] * y[0] - x[0] * y[2];
	res[2] = x[0] * y[1] - x[1] * y[0];

}

} // namespace Math