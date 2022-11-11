#pragma once
#include <string>
#include <math.h>
#include <ostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#define PI 3.14159265358979323846
class Math
{
public:
  static constexpr double Clamp(double, double = 0.0, double = 1.0);

  // Minimum and maximum
  static constexpr double Min(double, double);
  static constexpr double Max(double, double);
  static constexpr double Min(double, double, double);
  static constexpr double Max(double, double, double);

  static constexpr double DegreeToRadian(double);
  static constexpr double RadianToDegree(double);
};

/*!
\brief Clamp a double value between two bounds.
\param x Input value.
\param a, b Lower and upper bounds.
*/
inline constexpr double Math::Clamp(double x, double a, double b)
{
  return (x < a ? a : (x > b ? b : x));
}

/*!
\brief Minimum of two reals.
\param a, b Real values.
*/
inline constexpr double Math::Min(double a, double b)
{
  return (a < b ? a : b);
}

/*!
\brief Maximum of two reals.
\param a, b Real values.
*/
inline constexpr double Math::Max(double a, double b)
{
  return (a > b ? a : b);
}

/*!
\brief Maximum of three reals.
\param a, b, c Real values.
*/
inline constexpr double Math::Max(double a, double b, double c)
{
  return Math::Max(Math::Max(a, b), c);
}

/*!
\brief Minimum of three reals.
\param a, b, c Real values.
*/
inline constexpr double Math::Min(double a, double b, double c)
{
  return Math::Min(Math::Min(a, b), c);
}

/*!
\brief Convert degrees to randians.
\param a Angle in degrees.
*/
inline constexpr double Math::DegreeToRadian(double a)
{
  return a * 3.14159265358979323846 / 180.0;
}

/*!
\brief Convert radian to degrees.
\param a Angle in radian.
*/
inline constexpr double Math::RadianToDegree(double a)
{
  return a * 180.0 / 3.14159265358979323846;
}

// Class
class Vector
{
protected:
  double c[3]; //!< Components.
public:
  //! Empty 
  Vector() {}

  explicit Vector(double);
  explicit Vector(double, double, double);

  // Access members
  double& operator[] (int);
  double operator[] (int) const;

  // Unary operators
  Vector operator+ () const;
  Vector operator- () const;

  // Assignment operators
  Vector& operator+= (const Vector&);
  Vector& operator-= (const Vector&);
  Vector& operator*= (const Vector&);
  Vector& operator/= (const Vector&);
  Vector& operator*= (double);
  Vector& operator/= (double);

  // Binary operators
  friend int operator> (const Vector&, const Vector&);
  friend int operator< (const Vector&, const Vector&);

  friend int operator>= (const Vector&, const Vector&);
  friend int operator<= (const Vector&, const Vector&);

  // Binary operators
  friend Vector operator+ (const Vector&, const Vector&);
  friend Vector operator- (const Vector&, const Vector&);

  friend constexpr double operator* (const Vector&, const Vector&);

  friend Vector operator* (const Vector&, double);
  friend Vector operator* (double, const Vector&);
  friend Vector operator/ (const Vector&, double);

  friend Vector operator/ (const Vector&, const Vector&);

  // Boolean functions
  friend int operator==(const Vector&, const Vector&);
  friend int operator!=(const Vector&, const Vector&);

  // Norm
  friend double Norm(const Vector&);
  friend double SquaredNorm(const Vector&);
  friend double Distance(const Vector& v1, const Vector& v2);

  friend void Normalize(Vector&);
  friend Vector Normalized(const Vector&);

  // Compare functions
  static Vector Min(const Vector&, const Vector&);
  static Vector Max(const Vector&, const Vector&);

  // Abs
  friend Vector Abs(const Vector&);

  // Orthogonal and orthonormal vectors
  Vector Orthogonal() const;
  void Orthonormal(Vector&, Vector&) const;

  friend Vector Lerp(const Vector&, const Vector&, double);
  static Vector Bilinear(const Vector&, const Vector&, const Vector&, const Vector&, double, double);

  // Scale
  Vector Scaled(const Vector&) const;
  Vector Inverse() const;

  std::string vectorToString()const;

  friend std::ostream& operator<<(std::ostream&, const Vector&);

public:
  static const Vector Null; //!< Null vector.
  static const Vector X; //!< Vector(1,0,0).
  static const Vector Y; //!< Vector(0,1,0).
  static const Vector Z; //!< Vector(0,0,1).
};

/*!
\brief Create a vector with the same coordinates.
\param a Real.
*/
inline Vector::Vector(double a)
{
  c[0] = c[1] = c[2] = a;
}

/*!
\brief Create a vector with argument coordinates.
\param a,b,c Coordinates.
*/
inline Vector::Vector(double a, double b, double c)
{
  Vector::c[0] = a;
  Vector::c[1] = b;
  Vector::c[2] = c;
}

//! Gets the i-th coordinate of vector.
inline double& Vector::operator[] (int i)
{
  return c[i];
}

//! Returns the i-th coordinate of vector.
inline double Vector::operator[] (int i) const
{
  return c[i];
}

// Unary operators

//! Overloaded.
inline Vector Vector::operator+ () const
{
  return *this;
}

//! Overloaded.
inline Vector Vector::operator- () const
{
  return Vector(-c[0], -c[1], -c[2]);
}

// Assignment unary operators

//! Destructive addition.
inline Vector& Vector::operator+= (const Vector& u)
{
  c[0] += u.c[0]; c[1] += u.c[1]; c[2] += u.c[2];
  return *this;
}

//! Destructive subtraction.
inline Vector& Vector::operator-= (const Vector& u)
{
  c[0] -= u.c[0]; c[1] -= u.c[1]; c[2] -= u.c[2];
  return *this;
}

//! Destructive scalar multiply.
inline Vector& Vector::operator*= (double a)
{
  c[0] *= a; c[1] *= a; c[2] *= a;
  return *this;
}

/*!
\brief Scale a vector.
\param a Scaling vector.
*/
inline Vector Vector::Scaled(const Vector& a) const
{
  return Vector(c[0] * a[0], c[1] * a[1], c[2] * a[2]);
}

/*!
\brief Inverse of a vector.

This function inverses the components of the vector. This is the same as:
\code
Vector v=Vector(1.0/u[0],1.0/u[1],1.0/u[2]);
\endcode
*/
inline Vector Vector::Inverse() const
{
  return Vector(1.0 / c[0], 1.0 / c[1], 1.0 / c[2]);
}

//! Destructive division by a scalar.
inline Vector& Vector::operator/= (double a)
{
  c[0] /= a; c[1] /= a; c[2] /= a;
  return *this;
}

/*!
\brief Destructively scale a vector by another vector.

This is the same as Scale:
\code
Vector u(2.0,-1.0,1.0);
u=u.Scaled(Vector(3.0,1.0,2.0)); // u*=Vector(3.0,1.0,2.0);
\endcode
*/
inline Vector& Vector::operator*= (const Vector& u)
{
  c[0] *= u.c[0]; c[1] *= u.c[1]; c[2] *= u.c[2];
  return *this;
}

//! Destructively divide the components of a vector by another vector.
inline Vector& Vector::operator/= (const Vector& u)
{
  c[0] /= u.c[0]; c[1] /= u.c[1]; c[2] /= u.c[2];
  return *this;
}

//! Compare two vectors.
inline int operator> (const Vector& u, const Vector& v)
{
  return ((u.c[0] > v.c[0]) && (u.c[1] > v.c[1]) && (u.c[2] > v.c[2]));
}

//! Compare two vectors.
inline int operator< (const Vector& u, const Vector& v)
{
  return ((u.c[0] < v.c[0]) && (u.c[1] < v.c[1]) && (u.c[2] < v.c[2]));
}

//! Overloaded
inline int operator>= (const Vector& u, const Vector& v)
{
  return ((u.c[0] >= v.c[0]) && (u.c[1] >= v.c[1]) && (u.c[2] >= v.c[2]));
}

//! Overloaded
inline int operator<= (const Vector& u, const Vector& v)
{
  return ((u.c[0] <= v.c[0]) && (u.c[1] <= v.c[1]) && (u.c[2] <= v.c[2]));
}

//! Adds up two vectors.
inline Vector operator+ (const Vector& u, const Vector& v)
{
  return Vector(u.c[0] + v.c[0], u.c[1] + v.c[1], u.c[2] + v.c[2]);
}

//! Difference between two vectors.
inline Vector operator- (const Vector& u, const Vector& v)
{
  return Vector(u.c[0] - v.c[0], u.c[1] - v.c[1], u.c[2] - v.c[2]);
}

//! Scalar product.
inline constexpr double operator* (const Vector& u, const Vector& v)
{
  return (u.c[0] * v.c[0] + u.c[1] * v.c[1] + u.c[2] * v.c[2]);
}

//! Right multiply by a scalar.
inline Vector operator* (const Vector& u, double a)
{
  return Vector(u.c[0] * a, u.c[1] * a, u.c[2] * a);
}

//! Left multiply by a scalar.
inline Vector operator* (double a, const Vector& v)
{
  return v * a;
}

//! Cross product.
inline Vector operator/ (const Vector& u, const Vector& v)
{
  return Vector(u.c[1] * v.c[2] - u.c[2] * v.c[1], u.c[2] * v.c[0] - u.c[0] * v.c[2], u.c[0] * v.c[1] - u.c[1] * v.c[0]);
}

//! Left multiply by a scalar
inline Vector operator/ (const Vector& u, double a)
{
  return Vector(u.c[0] / a, u.c[1] / a, u.c[2] / a);
}

// Boolean functions

//! Strong equality test.
inline int operator== (const Vector& u, const Vector& v)
{
  return ((u.c[0] == v.c[0]) && (u.c[1] == v.c[1]) && (u.c[2] == v.c[2]));
}

//! Strong difference test.
inline int operator!= (const Vector& u, const Vector& v)
{
  return (!(u == v));
}

/*!
\brief Compute the Euclidean norm of a vector.

This function involves a square root computation, it is in general more efficient to rely on
the squared norm of a vector instead.
\param u %Vector.
\sa SquaredNorm
*/
inline double Norm(const Vector& u)
{
  return sqrt(u.c[0] * u.c[0] + u.c[1] * u.c[1] + u.c[2] * u.c[2]);
}

/*!
\brief Compute the squared Euclidean norm of a vector.
\param u %Vector.
\sa Norm
*/
inline double SquaredNorm(const Vector& u)
{
  return (u.c[0] * u.c[0] + u.c[1] * u.c[1] + u.c[2] * u.c[2]);
}


inline double Distance(const Vector& v1, const Vector& v2) {
	double res = (v2[0] - v1[0]) * (v2[0] - v1[0]);
	res += (v2[1] - v1[1]) * (v2[1] - v1[1]);
	res += (v2[2] - v1[2]) * (v2[2] - v1[2]);
	return sqrt(res);
}
/*!
\brief Return a normalized vector.

Compute the inverse of its norm and scale the components.

This function does not check if the vector is null.
\param u %Vector.
*/
inline Vector Normalized(const Vector& u)
{
  return u * (1.0 / Norm(u));
}

/*!
\brief Computes the absolute value of a vector.
\param u %Vector.
*/
inline Vector Abs(const Vector& u)
{
  return Vector(u[0] > 0.0 ? u[0] : -u[0], u[1] > 0.0 ? u[1] : -u[1], u[2] > 0.0 ? u[2] : -u[2]);
}

/*!
\brief Return a vector with coordinates set to the minimum coordinates
of the two argument vectors.
*/
inline Vector Vector::Min(const Vector& a, const Vector& b)
{
  return Vector(a[0] < b[0] ? a[0] : b[0], a[1] < b[1] ? a[1] : b[1], a[2] < b[2] ? a[2] : b[2]);
}

/*!
\brief Return a vector with coordinates set to the maximum coordinates
of the two argument vectors.
*/
inline Vector Vector::Max(const Vector& a, const Vector& b)
{
  return Vector(a[0] > b[0] ? a[0] : b[0], a[1] > b[1] ? a[1] : b[1], a[2] > b[2] ? a[2] : b[2]);
}

/*!
\brief Linear interpolation between two vectors.
\param a,b Interpolated points.
\param t Interpolant.
*/
inline Vector Lerp(const Vector& a, const Vector& b, double t)
{
  return a + t * (b - a);
}

/*!
\brief Bi-linear interpolation between four vectors.

The values are given in trigonometric order.

\param a00,a10,a11,a01 Interpolated vectors.
\param u,v Interpolation coefficients.

\sa Math::Bilinear
*/
inline Vector Vector::Bilinear(const Vector& a00, const Vector& a10, const Vector& a11, const Vector& a01, double u, double v)
{
  return (1 - u) * (1 - v) * a00 + (1 - u) * (v)*a01 + (u) * (1 - v) * a10 + (u) * (v)*a11;
}

inline std::string Vector::vectorToString()const {
	return std::to_string(this-> c[0]) + ";" + std::to_string(this->c[1]) +";"+ std::to_string(this->c[2]);
}



//! Matrix class
class Matrix {
private:
	std::vector<double> values; //! values of the matrix
	int nbColumns;          //! Number of columns.
	int nbLines;         //!  Number of lines.

	/*!
	  \brief Initialize values of the matrix 
	  \param values is a table with values to initialize the matrix. If values is nullptr, initialize with zeros.
	*/
	void initMatrix(double* values);
public:
	/*!
		\brief Creat a 3x3 matrix with only zeros.
	*/
	Matrix() { //defaut : matrice 3x3 avec des zeros
		nbColumns = 3;
		nbLines = 3;
		initMatrix(nullptr);
	}
	/*!
	  \brief Creat a matrix with only zeros.
	  \param lines is the number of lines
	  \param cols is the number of colums.
	*/
	Matrix(int lines, int cols) {
		this->nbColumns = cols;
		this->nbLines = lines;
		initMatrix(nullptr);
	}
	/*!
	  \brief Creat a matrix with values given.
	  \param lines is the number of lines
	  \param cols is the number of colums.
	  \param values to initialize with initMatrix function.
	*/
	Matrix(int lines, int cols, double* values) {
		this->nbColumns = cols;
		this->nbLines = lines;
		initMatrix(values);
	}
	/*!
	  \brief Copy constructor
	  \param other is the matrix to copy.
	*/
	Matrix(const Matrix& other) {
		this->nbColumns = other.nbColumns;
		this->nbLines = other.nbLines;
		this->values = other.values;
	}
	/*!
	  \brief Creat a 3x1 Matrix with a Vector
	  \param v is the vector to copy.
	*/
	Matrix(const Vector& v) {
		this->nbColumns = 1;
		this->nbLines = 3;
		values.push_back(v[0]);
		values.push_back(v[1]);
		values.push_back(v[2]);
	}
	/*!
	  \brief Destroy the matrix
	*/
	virtual void destroy() {
		values.clear();
		nbColumns = 0;
		nbLines = 0;
	}
	/*!
	  \brief Destructor
	*/
	virtual ~Matrix() {
		destroy();
	}

	/*!
	  \brief Cast a matrix into a std::vector
	  \return a copy of the member values.
	*/
	std::vector<double> MatrixToVector()const {
		return values;
	}
	
	/*!
	  \brief Cast a matrix into a table.
	  \return a table with values of the matrix.
	*/
	double* MatrixToArray()const;

	/*!
	  \brief Getter for nblines
	*/
	int getNbLines()const {
		return nbLines;
	}

	/*!
	  \brief Getter for nbcols
	*/
	int getNbCols()const {
		return nbColumns;
	}


	/*!
	  \brief set a value of the matrix
	  \param i line number
	  \param j column number
	  \param value is the new value
	*/
	void setValue(int i, int j, double value) {
		if(i< nbLines && j < nbColumns)
			values[i * nbColumns + j] = value;
	}

	/*!
	  \brief set all values of the matrix
	  \param values is the table with new values
	  \param size is the number of values
	*/
	void setValues(double* values, int size);

	/*!
	  \brief set all values of the matrix
	  \param v is the std::vector with new values
	*/
	void setValues(std::vector<double> v) {
		this->values = v;
	}

	/*!
	  \brief getter for a value of the matrix
	  \param i,j coordinates of the value
	  \return the value
	*/
	double getValue(int i, int j)const {
		if (i < nbLines && j < nbColumns)
			return values[i * nbColumns + j];
		else {
			throw "Matrix::getValue : out of range";
		}
	}

	/*!
	  \brief copy an other matrix
	  \param other is the matrix to copy
	*/
	void copy(const Matrix& other) {
		destroy();
		nbColumns = other.nbColumns;
		nbLines = other.nbLines;
		values = other.values;
	}
	/*!
	  \brief Creat a clone 
	  \return the clone of this matrix
	*/
	Matrix * clone()const {
		Matrix* m = new Matrix(*this);
		return m;
	}

	//methodes

	//MATRICES DE TRANSFORMATIONS


	/*!
	  \brief this matrix became a dilation matrix
	  \param v is the vector of dilation
	*/
	void homothetieMatrix(const Vector& v) {
		double values[9] = {v[0],0,0,0,v[1],0,0,0,v[2]};
		destroy();
		nbColumns = 3;
		nbLines = 3;
		initMatrix(values);
	}
	/*!
	  \brief this matrix became a rotation matrix
	  \param p is the rotation axis
	  \param angle is the rotation angle
	*/
	void rotateMatrix(const Vector& p, double angle) {
		double values[9];

		values[0] = p[0] * p[0] * (1 - cos(angle)) + cos(angle); //Ux² * (1-c) + c
		values[1] = p[0] * p[1] * (1 - cos(angle)) - p[2] * sin(angle);
		values[2] = p[0] * p[2] * (1 - cos(angle)) + p[1] * sin(angle);
		values[3] = p[0] * p[1] * (1 - cos(angle)) + p[2] * sin(angle); //Ux * Uy * (1-c) + Uz * s
		values[4] = p[1] * p[1] * (1 - cos(angle)) + cos(angle);
		values[5] = p[1] * p[2] * (1 - cos(angle)) - p[0] * sin(angle);
		values[6] = p[0] * p[2] * (1 - cos(angle)) - p[1] * sin(angle);
		values[7] = p[1] * p[2] * (1 - cos(angle)) + p[0] * sin(angle);
		values[8] = p[2] * p[2] * (1 - cos(angle)) + cos(angle);

 		destroy();
		nbColumns = 3;
		nbLines = 3;
		initMatrix(values);
	}

	/*!
	  \brief determinant of the matrix
	  \return the determinant of the matrix
	*/
	double getDeterminant()const;

	/*!
	  \brief cofactor of the matrix
	  \return the cofactor of the matrix
	*/
	Matrix* getCofactor()const;

	/*!
	  \brief transpose the matrix
	*/
	void transposition();
	
	/*!
	  \brief inverse matrix
	  \return return the inverse matrix of this matrix
	*/
	Matrix* getInverse()const;

	//! std::string for the matrix
	std::string toString()const;
	//operator

	//! + between this matrix and an other
	Matrix* operator +(const Matrix& other)const;

	//! - between this matrix and an other
	Matrix* operator -(const Matrix& other)const;

	//! + between this matrix and a double
	Matrix* operator *(double s)const;

	//! * between this matrix and an other
	Matrix* operator *(const Matrix& other)const;
	
	//! * between this matrix and a vector
	Vector operator*(const Vector& v)const;

	//! = with a matrix
	Matrix& operator =(const Matrix& other) {
		if (this != &other) {
			this->copy(other);
		}
		return *this;
	}

	//! equality between this matrix and an other
	bool operator ==(const Matrix& other)const;

	//! inequality between this matrix and an other
	bool operator !=(const Matrix& other)const {
		bool res = *this == other;
		return !res;
	}

	//! return a line of the matrix
	std::vector<double> operator [](int i)const;
};

//! print a matrix 

inline std::ostream& operator<<(std::ostream& s, const Matrix& m)
{
	s << m.toString();
	return s;
}