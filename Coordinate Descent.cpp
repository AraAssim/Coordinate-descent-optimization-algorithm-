#define NMAX 50  // Defining the maximum number of iterations. 

#include<iostream>
#include<cmath>
#include<cstdio>

using namespace std;

// The function is written here below:
double f(double x, double y) {
	return (11 + pow(x, 2) + pow(y, 2) - 2);
}

// The first derivative of the function with respect to x is
double f_dx(double x, double y) {
	return 2 * x;

}

// The first derivative of the function with respect to y is
double f_dy(double x, double y) {
	return 2 * y;
}

// The following functions below are defined to do the coordinate descent search process
double g(double x, double y, double alpha) {
	return f(x - alpha * f_dx(x, y), y - alpha * f_dy(x, y));
}

double Function(double x, double y) {
	return sqrt((f_dx(x, y)) * (f_dx(x, y)) + (f_dy(x, y)) * (f_dy(x, y)));
}

// Half division method to find the minimum in coordinate descent method 
double Dihotomia(double a0, double b0, double epsilon, double x, double y)
{
	//Number of steps
	int k;
	//Deviation 
	double lk, mk;
	//The value by which we deviate from the middle of the segment
	double delta = 0.5 * epsilon;
	//Minimum point 
	double x_;
	//Minimum localization segment
	double ak = a0, bk = b0;
	k = 1;
	//While the length of the segment is greater than the specified accuracy
	do {
		//We take the middle point 
		lk = (ak + bk - delta) / 2;
		mk = (ak + bk + delta) / 2;

		k++;
		//Checking where the minimum point falls to the left of the partition or to the right and select corresponding point 
		if (g(x, y, lk) <= g(x, y, mk)) {
			// The right boundary of the localization segment is mk


			bk = mk;
		}
		else {
			//The left boundary of the localization segment is mk
			ak = lk;
		}
	} while ((bk - ak) >= epsilon);

	x_ = (ak + bk) / 2; //minimum point

	return x_;
}


// Coordinate descent method
double CoordinateDescent(double bx, double by, double epsilon)
{
	double x[NMAX];
	double y[NMAX];
	double alpha[NMAX];
	int k;

	// Initial approximations
	x[0] = bx;
	y[0] = by;

	cout << "Results:" << endl << "x(" << 0 << "): (" << x[0] << ", " << y[0] << ")" << endl;

	for (k = 0; ; k++)
	{
		//We find alpha_k of at least function g on the interval -10000,100000
		alpha[k] = Dihotomia(-10000, 100000, epsilon, x[k], y[k]);

		//We calculate u [k]
		x[k + 1] = x[k] - alpha[k] * f_dx(x[k], y[k]);
		y[k + 1] = y[k] - alpha[k] * f_dy(x[k], y[k]);

		cout << "x(" << k + 1 << "): " << "(" << x[k + 1] << ", " << y[k + 1] << ")" << endl
			<< " f(" << x[k + 1] << ", " << y[k + 1] << ") = " << f(x[k + 1], y[k + 1]) << endl;
		if (k > 1) {
			//Checking the stop condition 
			if (Function(x[k + 1] - x[k], y[k + 1] - y[k]) < epsilon) {
				break;
			}
		}
	}

	cout << " \n That is the minimum point for (epsilon=" << epsilon << ") : \n" <<
		"f(x,y)= (" << x[k + 1] << ", " << y[k + 1] << ") = " << f(x[k + 1], y[k + 1]) << endl;
	return f(x[k + 1], y[k + 1]);

}
int main() {
	double x, y, epsilon;
	cout << "Please enter initial x value\n"; cin >> x;
	cout << "Please enter initial y value\n"; cin >> y;
	cout << "Please enter the error rate/epsilon\n"; cin >> epsilon;
	CoordinateDescent(x, y, epsilon);

	system("pause");
	return 0;
}
