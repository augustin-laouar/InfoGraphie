#pragma once
#include "mathematics.h"
#include <QtGui/QImage> 
#include <random>
//! class Sphere
class Sphere
{
protected:
	//! Center of the sphere
	Vector center;
	//! Radius of the sphere
	double radius;

public:
	//! Constructor with the center and the radius
	Sphere(Vector c, double r) {
		center = c;
		radius = r;
	}

	//! Constructor with radius only. The center is (0,0,0)
	Sphere(double r) {
		center = Vector(0, 0, 0);
		radius = r;
	}
	//!Defaut constructor. Center is (0,0,0) and radius is 1.0 .
	Sphere() {
		center = Vector(0, 0, 0);
		radius = 1.0;
	}
	//! Copy constructor with an other sphere
	Sphere(const Sphere& other) {
		center = other.center;
		radius = other.radius;
	}
	//! Destructor. Do nothing.
	virtual ~Sphere() {}

	/*!
	  \brief Setter for the center of the sphere
	  \param c is the new center 
	*/
	void setCenter(const Vector& c) {
		center = c;
	}
	/*!
	  \brief Setter for the center of the sphere
	  \param x,y,z coordinates for the new center
	*/
	void setCenter(double x, double y, double z) {
		center = Vector(x, y, z);
	}
	//! Setter for the radius of the sphere
	void setRadius(double r) {
		radius = r;
	}
	//! Getter for the center of the sphere
	Vector getCenter()const {
		return center;
	}
	//! Getter for the radius of the sphere
	double getRadius()const {
		return radius;
	}
};

//! Class disk
class Disk
{
protected:
	//! Center of the disk
	Vector center;
	//! Radius of the disk
	double radius;

public:
	//! Constructor with the center and the radius
	Disk(Vector c, double r) {
		center = c;
		radius = r;
	}
	//! Constructor with the radius only. Center is (0,0,0)
	Disk(double r) {
		center = Vector(0, 0, 0);
		radius = r;
	}
	//! Defaut constructor. Center is (0,0,0) and radius is 1.0 .
	Disk() {
		center = Vector(0, 0, 0);
		radius = 1.0;
	}
	//! Copy constructor with an other disk.
	Disk(const Disk& other) {
		center = other.center;
		radius = other.radius;
	}
	//! Destructor. Do nothing.
	virtual ~Disk() {}


	/*!
		\brief Setter for the center of the disk
		\param c is the new center
		*/
	void setCenter(const Vector& c) {
		center = c;
	}
	/*!
	  \brief Setter for the center of the disk
	  \param x,y,z coordinates for the new center
	*/
	void setCenter(double x, double y, double z) {
		center = Vector(x, y, z);
	}
	//! Setter for the radius of the disk
	void setRadius(double r) {
		radius = r;
	}
	//! Getter for the center of the disk
	Vector getCenter()const {
		return center;
	}
	//! Getter for the radius of the disk
	double getRadius()const {
		return radius;
	}
};

//! Class cylinder
class Cylinder
{
protected:
	//!Radius of the cylinder
	double radius;
	//! Height of the cylinder
	double height;
	//! Center of the cylinder.
	Vector center;
public:
	//! Constructor with the center, the radius and the height
	Cylinder(Vector center,double r, double h) : radius(r), height(h), center(center) {}

	//! Default constructor, center is (0,0,0), radius is 1.0 and height is 1.0
	Cylinder() {
		center = Vector(0, 0, 0);
		radius = 1.0;
		height = 1.0;
	}

	//! Detructor. Do nothing.
	virtual ~Cylinder() {}

	//getters

	//! Getter for the radius
	double getRadius()const {
		return this->radius;
	}

	//! Getter for the height
	double getHeight()const {
		return this->height;
	}

	//! Getter for the center
	Vector getCenter()const {
		return center;
	}

	//setters
	//! Setter for the radius
	void setRadius(double r) {
		this->radius = r;
	}

	//! Setter for the height
	void setHeight(double height)
	{
		this->height = height;
	}

	//! Setter of the center
	void setCenter(const Vector& c)
	{
		this->center = c;
	}

};

//pour capsule : un cylindre et deux spheres
//! class Capsule
class Capsule {
protected:
	//! Radius of the capsule
	double radius;
	//! Height of the capsule
	double height;
	//! Center of the capsule. 
	Vector center;
public:

	//! Constructor of the capsule with the center, the radius and the height
	Capsule(Vector center ,double r, double h) : radius(r), height(h), center(center) {}

	//! Default constructor, center is (0,0,0), radius is 1.0 and height is 1.0
	Capsule() {
		center = Vector(0, 0, 0);
		radius = 1.0;
		height = 1.0;
	}
	//! Detructor. Do nothing.
	virtual ~Capsule() {}

	//getters

	//! Getter for the radius
	double getRadius()const {
		return this->radius;
	}

	//! Getter for the height
	double getHeight()const {
		return this->height;
	}

	//! Getter for the center
	Vector getCenter()const {
		return center;
	}

	//setters

	//! Setter for the radius
	void setRadius(double r) {
		this->radius = r;
	}

	//! setter for the height
	void setHeight(double height)
	{
		this->height = height;
	}

	//! Setter for the center
	void setCenter(const Vector& c)
	{
		this->center = c;
	}

};

//! Class Tore
class Tore {
protected:
	//!Radius of the ring
	double innerRadius;
	//! Total radius
	double outerRadius;
	//! Center of the tore
	Vector center;

public:
	//! Constructor with the radius of the ring, the total radius and the center
	Tore(Vector center,double innerRad,double outerRad): innerRadius(innerRad), outerRadius(outerRad),center(center) {}

	//! Detructor. Do nothing.
	virtual ~Tore() {}
	//! Getter for the ring's radius
	double getInnerRadius()const {
		return this->innerRadius;
	}
	//! Getter for the total radius
	double getOuterRadius()const {
		return this->outerRadius;
	}
	//! Getter for the center
	Vector getCenter()const {
		return center;
	}

	//setters

	//! Setter for the ring's radius
	void setInnerRadius(double r) {
		this->innerRadius = r;
	}
	//! Setter for the total radius
	void setOuterRadius(double r) {
		this->outerRadius = r;
	}
	//! Setter for the center
	void setCenter(const Vector& c)
	{
		this->center = c;
	}


};

//! Class HeightField
class HeightField {
private:
	//! a is the bottom-left vertex and b is the top-right vertex of the field in 2D .
	Vector a, b;
	//! nx and ny are the number of vertices ( and the number of pixels ) on a line and a column respectivly.
	int nx, ny; // nombre de points
	//! List of the height of the associated point.
	std::vector<double> z;
	//! Ratio to get the real X position of a vertex.
	double scaleX;
	//! Ratio to get the real Y position of a vertex.
	double scaleY;
public:
	/*!
	  \brief Constructor for a heightField with parameters.
	  \param qi is the QImage to generate the field
	  \param a is the bottom left point we want for the field
	  \param b is the top right vertex we want for the field
	  \param scaleZ is the ratio to get the real height of our vertices ( because grey level are only from 0 to 255 ).
	*/
	HeightField(const QImage& qi, const Vector& a, const Vector& b, double scaleZ);
	//! Copy constructor with an other heightfield
	HeightField(const HeightField& h) : a(h.a), b(h.b), nx(h.nx), ny(h.ny), z(h.z), scaleX(h.scaleX), scaleY(h.scaleY) {}

	//! Destructor. Destroy the std::vector.
	virtual ~HeightField() { z.clear(); }

	//! Return the index of the vertex (i,j) in the std::vector .
	int getIndex(int i, int j)const {// on stock ligne par ligne 
		return nx * i + j;
	}

	//! Getter for the bottom left vertex 
	Vector getA()const {
		return a;
	}

	//! Getter for the top right vertex 
	Vector getB()const {
		return b;
	}

	//! Getter for number of vertices on a line
	int getNx()const {
		return nx;
	}

	//! Getter for the number of vertices on a column
	int getNy()const {
		return ny;
	}

	//! Return the vertex in 3D at the point (i,j) (in 2D)
	Vector getPoint(int i, int j)const {
		return Vector(i * scaleX + a[0], j * scaleY + a[1], z.at(getIndex(i, j)));
	}
	//! Getter for the total number of points
	int getNbPoints()const {
		return nx * ny;
	}

	//! Return the height of the vertex at the point (i,j).
	double height(int i, int j)const {
		if (i < 0 || j < 0 || i >= nx || j >= ny)
			return 0;
		return getPoint(i, j)[2];
	}

	//! Return the gradient vector of the vertex (i,j) in the plane.
	Vector gradient(int i, int j)const {
		return Vector(height(i + 1, j) - height(i - 1, j), height(i, j + 1) - height(i, j - 1), 0);
	}

	//! Return the unitary vector of the vertex (i,j) in the plane.
	Vector N(int i, int j)const { // vecteur unitaire
		Vector g3 = -gradient(i, j);
		g3 += Vector(0, 0, 1);
		Vector n = g3 / Norm(g3);
		return n;
	}

	/*!
	  \brief Creat relief on the field.
	  \param c is the center of the disk where we apply the relief
	  \param r is the disk's radius
	  \param h is the height of our deformation. Can be nagative.
	*/
	void createRelief(const Vector& c, double r, double h);

	/*!
	  \brief Creat a bump. Their is no smoothing around.
	  \param c is the center of the disk where we apply the relief
	  \param r is the disk's radius
	  \param h is the height of our bump. Can be nagative.
	*/
	void bump(const Vector& c, double d, double h);
	/*!
	  \brief Terrace at the desired height with a smoothing
	  \param c is the center of the disk we terrace
	  \param r is the disk's radius where it will be set at the height we want
	  \param r2 is a bigger disk around the first disk, to do the smoothing.
	  \param h is the height for the earthworks.
	*/
	void terrasser(const Vector& c, double r, double r2, double h);

	/*!
	  \brief Similiar to terrasser function. It will terrace at the desired height with a somoothing.
	  \param c is the center of the disk we terrace
	  \param size is the disk's radius where it will be set at the height we want
	  \param h is the height for the earthworks. 0 by default.
	*/
	void brush(const Vector& c, double size, double h = 0) {
		terrasser(c, size, size + size * 0.5, h);
	}

	/*!
	  \brief Generate a random field by using the function "creatRelief" randomly.
	  \param nb is the number of deformations wanted.
	  \param minH, maxH are the minimum and maximum height we want for deformations
	  \param minR, maxR are the minimum and maximum radius we want for deformations.
	*/
	void randomField(int nb, double minH, double minR, double maxH, double maxR);
	/*!
	  \brief Generate a random field by using the function "creatRelief" randomly, into a defined perimeter.
	  \param nb is the number of deformations wanted.
	  \param minH, maxH are the minimum and maximum height we want for deformations
	  \param minR, maxR are the minimum and maximum radius we want for deformations.
	  \param minX, maxX, minY, maxY are used to define the perimeter wanted.
	*/
	void randomFieldInSquare(int nb, double minH, double minD, double maxH, double maxD, int minX, int maxX, int minY, int maxY);
	
};