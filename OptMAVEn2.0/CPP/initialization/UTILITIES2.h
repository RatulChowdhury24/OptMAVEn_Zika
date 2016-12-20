#ifndef UTILITIES_H
#define UTILITIES_H

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <vector>
using namespace std;


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


// Compute the norm of a coordinate vector of length 3.
inline float norm3(array<float, 3> x, bool squared=false)
{
	float ss = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];  // Hard-code for speed.
	if (squared)
		return ss;
	else
		return sqrt(ss);
}


// Compute the distance between two coordinates in 3D.
inline float dist3(array<float, 3> a1, array<float, 3> a2, bool squared=false)
{
	// Hard-code for speed.
	float d0 = a2[0] - a1[0]; float d1 = a2[1] - a1[1]; float d2 = a2[2] - a1[2];
	float ss = d0 * d0 + d1 * d1 + d2 * d2;
	if (squared)
		return ss;
	else
		return sqrt(ss);
}


// Normalize a coordinate vector of length 3 and return a new vector.
array<float, 3> normalize3(array<float, 3> coords)
{
	float norm = norm3(coords);  // Compute the norm of the input vector.
	array<float, 3> normalized;  // Declare a normalized vector.
	if (norm > 0)
	{
		for (int i = 0; i < 3; i++)
			normalized[i] = coords[i] / norm;
	}
	else
	{
		// Cannot normalize a zero vector.
		throw "Cannot normalize a zero vector.";
	}
	return normalized;
}


// Compute the 3x3 rotation matrix needed to rotate by the given angle along the given axis.
array<array<float, 3>, 3> rotation_matrix3(array<float, 3> axis, float angle, bool degrees=true)
{
	int i, j;  // Declare index variables.
	array<array<float, 3>, 3> matrix;  // Declare the rotation matrix.
	// Normalize the axis.
	array<float, 3> unit_axis = normalize3(axis);
	// If the angle is in degrees, convert to radians.
	float pi = 3.14159265359;
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
	return matrix;
}

// Perform an argument sort on a vector of floats: return a vector of the order of the elements in the given vector.
vector<int> argsort(vector<float> val)
{
	// Create a vector of ascending indices.
	vector<int> idx(val.size());
	int i = 0;
	generate(idx.begin(), idx.end(), [&]{return i++;});
	// Sort the vector based on the order of the elements in the first vector.
	sort(idx.begin(), idx.end(), [&](int i1, int i2) {return val[i1] < val[i2];});
	return idx;
}

// Perform an argument sort on a vector of floats: return a vector of the order of the elements in the given vector.
vector<int> argsort(vector<float>* val)
{
        // Create a vector of ascending indices.
        vector<int> idx(val->size());
        int i = 0;
        generate(idx.begin(), idx.end(), [&]{return i++;});
        // Sort the vector based on the order of the elements in the first vector.
        sort(idx.begin(), idx.end(), [&](int i1, int i2) {return val->at(i1) < val->at(i2);});
        return idx;
}

// Perform an argument sort on a vector of ints: return a vector of the order of the elements in the given vector.
vector<int> argsort(vector<int> val)
{
        // Create a vector of ascending indices.
        vector<int> idx(val.size());
        int i = 0;
        generate(idx.begin(), idx.end(), [&]{return i++;});
        // Sort the vector based on the order of the elements in the first vector.
        sort(idx.begin(), idx.end(), [&](int i1, int i2) {return val[i1] < val[i2];});
        return idx;
}

// Perform an argument sort on a vector of ints: return a vector of the order of the elements in the given vector.
vector<int> argsort(vector<int>* val)
{
        // Create a vector of ascending indices.
        vector<int> idx(val->size());
        int i = 0;
        generate(idx.begin(), idx.end(), [&]{return i++;});
        // Sort the vector based on the order of the elements in the first vector.
        sort(idx.begin(), idx.end(), [&](int i1, int i2) {return val->at(i1) < val->at(i2);});
        return idx;
}

#endif
