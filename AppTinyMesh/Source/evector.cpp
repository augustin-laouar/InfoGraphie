// Vector  

// Self include
#include "mathematics.h"

#include <iostream>

/*!
\class Vector mathematics.h
\brief Vectors in three dimensions.

Most binary operators have been overloaded as expected,
destructive operators, such as addition and subtraction
have been implemented and behave as one could expect.

<P><I>How do I compute the cross product of two vectors?</I>
<BR>Simply use the overloaded Vector::operator/, for instance
\code
Vector c=a/b; // Cross product
\endcode
computes the cross product of a and b.
<P><I>How can I get access to the x, y and z components of a vector?</I>
<BR>Use v[0], v[1] and v[2] to get access to the x, y and z components of a vector v respectively.
<P><I>How do I compute the normal of a triangle?</I>
<BR>Let a,b,c the vertices of the triangle, simply compute the cross product
\code
Vector n=(a-b)/(a-c);  // Cross product
\endcode
or use the member function of the Triangle class:
\code
Vector n=Triangle(a,b,c).Normal(); // Compute the normal
\endcode
<P><I>How can I sort the three elements in a vector?</I>
<BR>Use Vector::Sort() as follows:
\code
Vector s=Vector(2.0,3.0,-1.0).Sort(); // Sort components in ascending order
\endcode

<P><I>How do I perform bi-linear interpolation on vectors?</I>
<BR>Use Vector::Bilinear() with four vectors and bilinear coordinates.
Alternatively, some geometric classes implement bilinear interpolation,
such as Quadrangle::Vertex().
*/

const Vector Vector::Null = Vector(0.0, 0.0, 0.0);
const Vector Vector::X = Vector(1.0, 0.0, 0.0);
const Vector Vector::Y = Vector(0.0, 1.0, 0.0);
const Vector Vector::Z = Vector(0.0, 0.0, 1.0);

/*!
\brief Normalize a vector, computing the inverse of its norm and scaling
the components.

This function does not check if the vector is null,
which might resulting in errors.
*/
void Normalize(Vector& u)
{
  u *= 1.0 / Norm(u);
}

/*!
\brief Returns a vector orthogonal to the argument vector.

The returned orthogonal vector is not computed randomly.
First, we find the two coordinates of the argument vector with
maximum absolute value. The orthogonal vector is defined by
swapping those two coordinates and changing one sign, whereas
the third coordinate is set to 0.

The returned orthogonal vector lies in the plane orthogonal
to the first vector.
*/
Vector Vector::Orthogonal() const
{
  Vector a = Abs(*this);
  int i = 0;
  int j = 1;
  if (a[0] > a[1])
  {
    if (a[2] > a[1])
    {
      j = 2;
    }
  }
  else
  {
    i = 1;
    j = 2;
    if (a[0] > a[2])
    {
      j = 0;
    }
  }
  a = Vector::Null;
  a[i] = c[j];
  a[j] = -c[i];
  return a;
}

/*!
\brief Overloaded output-stream operator.
\param u Vector.
\param s Stream.
*/
std::ostream& operator<<(std::ostream& s, const Vector& u)
{
  s << "Vector(" << u.c[0] << ',' << u.c[1] << ',' << u.c[2] << ')';
  return s;
}

/*!
\brief Given a vector, creates two vectors xand y that form an orthogonal basis.

This algorithm pickes the minor axis in order to reduce numerical instability
\param x, y Returned vectors such that (x,y,n) form an orthonormal basis (provided n is normalized).
*/
void Vector::Orthonormal(Vector& x, Vector& y) const
{
  x = Normalized(Orthogonal());
  y = Normalized(*this / x);
}



//MATRIX

void Matrix::initMatrix(double* values) { //initialise avec les valeurs. Si nullptr, init avec des 0
    for (int i = 0; i < nbColumns * nbLines; i++) {
        if (values != nullptr)
            this->values.push_back(values[i]);
        else
            this->values.push_back(0);
    }
}

double* Matrix::MatrixToArray()const {
    double* res = new double[values.size()];
    for (int i = 0; i < values.size(); i++) {
        res[i] = values[i];
    }
    return res;
}

void Matrix::setValues(double* values, int size) {
    this->values.clear();
    for (int i = 0; i < size; i++) {
        this->values.push_back(values[i]);
    }
}

double Matrix::getDeterminant()const
{
	if (this->nbColumns != this->nbColumns)
	{
		throw "Matrix not quadratic";
	}
	int dimension = this->nbColumns;
	if (dimension == 0)
		return 1;
	if (dimension == 1)
		return getValue(0, 0);

	if (dimension == 2)
	{
		return getValue(0, 0) * getValue(1, 1) - getValue(0, 1) * getValue(1, 0);
	}
	double result = 0;
	int sign = 1;
	for (int i = 0; i < dimension; i++) {

		//Submatrix

		Matrix subVect(dimension - 1, dimension - 1);
		for (int m = 1; m < dimension; m++) {
			int z = 0;
			for (int n = 0; n < dimension; n++) {
				if (n != i) {
					subVect.setValue(m - 1, z, getValue(m, n));
					z++;
				}
			}
		}



		result = result + sign * getValue(0, i) * subVect.getDeterminant();
		sign = -sign;
	}

	return result;

}

Matrix* Matrix::getCofactor()const {

	if (this->nbColumns == this->nbLines) {
		int size = this->nbColumns;

		Matrix* solution = new Matrix(size, size);

		Matrix* subVect = new Matrix(size - 1, size - 1);

		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {

				int p = 0;
				for (int x = 0; x < size; x++) {
					if (x == i) {
						continue;
					}
					int q = 0;

					for (int y = 0; y < size; y++) {
						if (y == j) {
							continue;
						}
						subVect->setValue(p, q, this->getValue(x, y));
						q++;
					}
					p++;
				}
				solution->setValue(i, j, pow(-1, i + j) * subVect->getDeterminant());
			}
		}
		return solution;
	}
}


void Matrix::transposition()
{
	Matrix temp(this->nbLines, this->nbColumns);
	for (int i = 0; i < this->nbLines; i++)
	{
		for (int j = 0; j < this->nbColumns; j++)
		{
			temp.setValue(i, j, getValue(j, i));
		}
	}
	*this = temp;
}

Matrix* Matrix::getInverse()const {
	if (getDeterminant() == 0) {
		throw std::runtime_error("Determinant is 0");
	}
	double d = 1.0 / getDeterminant();
	Matrix* solution = new Matrix(this->nbLines, this->nbColumns);
	for (int i = 0; i < solution->nbLines; i++) {
		for (int j = 0; j < solution->nbColumns; j++) {
			solution->setValue(i, j, getValue(i, j));
		}
	}
	solution = solution->getCofactor();
	solution->transposition();

	for (int i = 0; i < nbLines; i++) {
		for (int j = 0; j < nbColumns; j++) {
			solution->setValue(i, j, solution->getValue(i, j) * d);
		}
	}
	return solution;
}

std::string Matrix::toString()const {
	std::string res = "Matrix " + std::to_string(nbColumns) + "x" + std::to_string(nbLines) + " : \n";
	for (int i = 0; i < nbLines; i++) {
		res += "| ";
		for (int j = 0; j < nbColumns; j++) {
			res += std::to_string(getValue(i, j)) + " ";
		}
		res += "|\n";
	}
	return res;
}

Matrix* Matrix::operator +(const Matrix& other)const {
	if (nbColumns != other.nbColumns || nbLines != other.nbLines)
		throw "Matrix::operator + : lines or column are inequals";
	Matrix* res = new Matrix(*this);
	for (int i = 0; i < nbLines; i++) {
		for (int j = 0; j < nbColumns; j++) {
			res->setValue(i, j, res->getValue(i, j) + other.getValue(i, j));
		}
	}
	return res;
}

Matrix* Matrix::operator -(const Matrix& other)const {
	if (nbColumns != other.nbColumns || nbLines != other.nbLines)
		throw "Matrix::operator - : lines or column are inequals";
	Matrix* res = new Matrix(*this);
	for (int i = 0; i < nbLines; i++) {
		for (int j = 0; j < nbColumns; j++) {
			res->setValue(i, j, res->getValue(i, j) - other.getValue(i, j));
		}
	}
	return res;
}

Matrix* Matrix::operator *(double s)const {
	Matrix* res = new Matrix(*this);
	for (int i = 0; i < nbLines; i++) {
		for (int j = 0; j < nbColumns; j++) {
			res->setValue(i, j, res->getValue(i, j) * s);
		}
	}
	return res;
}

Matrix* Matrix::operator *(const Matrix& other)const {
	if (nbColumns != other.nbLines)
		throw "Matrix::operator * with matrix : lines and column are inequals";
	Matrix* res = new Matrix(nbLines, other.nbColumns);
	for (int i = 0; i < res->nbLines; i++) {
		for (int j = 0; j < res->nbColumns; j++) {
			double c = 0;
			for (int k = 0; k < nbColumns; k++) {
				c += this->getValue(i, k) * other.getValue(k, j);
			}
			res->setValue(i, j, c);
		}
	}
	return res;
}

Vector Matrix::operator*(const Vector& v)const {
	if (nbColumns != 3) {
		throw "Matrix::operator * with vector : matrix has not 3 cols";
	}
	Matrix temp(v);
	temp = *(*this * temp);
	Vector res(temp.getValue(0, 0), temp.getValue(1, 0), temp.getValue(2, 0));
	return res;
}



bool Matrix::operator ==(const Matrix& other)const {
	if (nbColumns != other.nbColumns || nbLines != other.nbLines) {
		return false;
	}
	for (int i = 0; i < this->nbLines; i++) {
		for (int j = 0; j < this->nbColumns; j++) {
			if (getValue(i, j) != other.getValue(i, j))
				return false;
		}
	}
	return true;
}

std::vector<double> Matrix::operator [](int i)const {
	if (i > nbColumns) {
		throw "Matrix::operator []  : out of range";
	}
	std::vector<double> res;
	for (int j = 0; j < nbColumns; j++) {
		res.push_back(getValue(i, j));
	}
	return res;
}