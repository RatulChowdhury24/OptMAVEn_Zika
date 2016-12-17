#ifndef UTILITIES_H
#define UTILITIES_H

#include <cmath>
#include <iostream>
using namespace std;


// Compute the norm of a coordinate vector of length 3.
float norm3(float x[3], bool squared=false)
{
	float ss = 0.0;
	for (int i = 0; i < 3; i++)
		ss += x[i] * x[i];
	if (squared)
		return ss;
	else
		return sqrt(ss);
}


// Normalize a coordinate vector of length 3, mutating it.
float normalize3(float coords[3])
{
	float norm = norm3(coords);
	if (norm > 0)
	{
		for (int i = 0; i < 3; i++)
			coords[i] /= norm;
	}
	else
	{
		for (int i = 0; i < 3; i++)
			coords[i] = 0.0;
	}
	return norm;
}


// Normalize a coordinate vector of length 3 and return a new vector.
float normalize3(float coords[3], float normalized[3])
{
	float norm = norm3(coords);
	if (norm > 0)
	{
		for (int i = 0; i < 3; i++)
			normalized[i] = coords[i] / norm;
	}
	else
	{
		for (int i = 0; i < 3; i++)
			normalized[i] = 0.0;
	}
	return norm;
}


// Return the sign of a number.
int sign(float x)
{
	if (x > 0)
		return 1;
	else if (x < 0)
		return -1;
	else
		return 0;
}


// Compute the 3x3 rotation matrix needed to rotate by the given angle along the given axis.
void rotation_matrix3(float matrix[3][3], float axis[3], float angle, bool degrees=true)
{
	int i, j;  // Declare index variables.
	float pi = 3.14159265359;
	// Normalize the axis.
	float unit_axis[3];
	normalize3(axis, unit_axis);
	// If the angle is in degrees, convert to radians.
	if (degrees)
		angle = angle * pi / 180.0;
	// Compute sine and cosine of the angle.
	float sine = sin(angle); float cosine = cos(angle); float sub1cos = 1.0 - cosine;
	// Compute the rotation matrix. This formula comes from https://en.wikipedia.org/wiki/Rotation_matrix.
	float identity_cosine, cross_product, tensor_product;
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			if (i == j)
			{
				identity_cosine = cosine;
				cross_product = 0.0;
			}
			else
			{
				identity_cosine = 0.0;
				cross_product = unit_axis[3 - (i + j)];
				if (i + j == 2)
					cross_product *= sign(j - i);
				else
					cross_product *= sign(i - j);
			}		
			tensor_product = unit_axis[i] * unit_axis[j];
			matrix[i][j] = identity_cosine + cross_product * sine + tensor_product * sub1cos;
		}
	} 
}

#endif
