#ifndef OVERLOADS_H
#define OVERLOADS_H

// Vector - Vector
template <typename T>
inline std::vector<T> operator-(const std::vector<T> &LVec, const std::vector<T> &RVec) {

	// Check vectors are the same size
	if (LVec.size() != RVec.size()) {
		std::cout << "Vectors not right size for subtracting..." << std::endl;
		exit(-1);
	}

	// Get number of rows
	size_t rows = LVec.size();

	std::vector<T> ResVec;
	ResVec.reserve(rows);

	// Add left and right vectors
	for (size_t i = 0; i < rows; i++) {
		ResVec.push_back(LVec[i] - RVec[i]);
	}

	return ResVec;
}

// Scalar * Vector
template <typename T>
inline std::vector<T> operator*(const T LScalar, const std::vector<T> &RVec) {

	// Get number of rows
	size_t rows = RVec.size();

	std::vector<T> ResVec;
	ResVec.reserve(rows);

	// Multiply  left and right matrices
	for (size_t i = 0; i < rows; i++) {
		ResVec.push_back(LScalar * RVec[i]);
	}

	return ResVec;
}

// Vector * Vector
template <typename T>
inline T operator*(const std::vector<T> &LVec, const std::vector<T> &RVec) {

	// Check vectors are the same size
	if (LVec.size() != RVec.size()) {
		std::cout << "Vectors not right size for multiplying..." << std::endl;
		exit(-1);
	}

	// Get number of rows
	size_t rows = LVec.size();

	T res = 0.0;

	// Multiply  left and right matrices
	for (size_t i = 0; i < rows; i++) {
		res += LVec[i] * RVec[i];
	}

	return res;
}

#endif
