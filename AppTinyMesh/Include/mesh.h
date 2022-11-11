#pragma once

#include "box.h"
#include "ray.h"
#include "mathematics.h"
#include "basic_shapes.h"



// Triangle
class Triangle
{
protected:
  Vector p[3] = {Vector(0.0,0.0,0.0),Vector(1.0,0.0,0.0), Vector(0.0,1.0,0.0), }; //!< Array of vertices.
public:
  //! Empty.
  Triangle() {}
  explicit Triangle(const Vector&, const Vector&, const Vector&);

  //! Empty.
  ~Triangle() {}

  Vector operator[] (int) const;

  // Point in triangle
  Vector Vertex(double, double) const;

  // Intersection
  bool Intersect(const Ray&, double&, double&, double&) const;

  void Translate(const Vector&);

  // Geometry
  Vector Normal() const;
  Vector AreaNormal() const;
  Vector Center() const;

  double Area() const;
  double Aspect() const;
  Box GetBox() const;

  // Stream
  friend std::ostream& operator<<(std::ostream&, const Triangle&);

  double InscribedRadius() const;
  double CircumscribedRadius() const;
protected:
  static double epsilon; //!< Internal epsilon constant.
};

/*!
\brief Return the i-th vertex.
\param i Index.
*/
inline Vector Triangle::operator[] (int i) const
{
  return p[i];
}

//! Compute the barycenter of the triangle.
inline Vector Triangle::Center() const
{
  return (p[0] + p[1] + p[2]) / 3.0;
}

//! Compute the area of the triangle.
inline double Triangle::Area() const
{
  return 0.5 * Norm((p[0] - p[1]) / (p[2] - p[0]));
}

/*!
\brief Create a triangle.
\param a,b,c Vertices of the triangle.
*/
inline Triangle::Triangle(const Vector& a, const Vector& b, const Vector& c)
{
  p[0] = a;
  p[1] = b;
  p[2] = c;
}


class QString;

class Mesh
{
protected:
  std::vector<Vector> vertices; //!< Vertices.
  std::vector<Vector> normals;  //!< Normals.
  std::vector<int> varray;     //!< Vertex indexes.
  std::vector<int> narray;     //!< Normal indexes.
public:
  explicit Mesh();
  explicit Mesh(const std::vector<Vector>&, const std::vector<int>&);
  explicit Mesh(const std::vector<Vector>&, const std::vector<Vector>&, const std::vector<int>&, const std::vector<int>&);
  ~Mesh();

  void Reserve(int, int, int, int);

  Triangle GetTriangle(int) const;
  Vector Vertex(int) const;
  Vector Vertex(int, int) const;

  Vector Normal(int) const;

  int Triangles() const;
  int Vertexes() const;

  std::vector<int> VertexIndexes() const;
  std::vector<int> NormalIndexes() const;

  int VertexIndex(int, int) const;
  int NormalIndex(int, int) const;

  Vector operator[](int) const;

  Box GetBox() const;

  void Scale(double);

  void SmoothNormals();

  // Constructors from core classes
  //! Creat a mesh with a box
  explicit Mesh(const Box&);
  //! Creat a mesh with a disk
  explicit Mesh(const Disk&);
  //! Creat a mesh with a sphere
  explicit Mesh(const Sphere&);
  //! Creat a mesh with a cylinder
  explicit Mesh(const Cylinder&);
  //! Creat a mesh with a tore
  explicit Mesh(const Tore&);
  //! Creat a mesh with a capsule
  explicit Mesh(const Capsule&);
  //! Creat a mesh with a heightfield
  explicit Mesh(const HeightField&);


  //Transformations
  /*!
	\brief Merge with an other mesh
	\param a is the other mesh
  */
  void merge(const Mesh& a);
  /*!
	\brief apply a transformation to the mesh
	\param m is the matrix of transformation
  */
  void transformation(const Matrix& m);
  /*!
	\brief apply an homothety to the mesh
	\paramv is the vector of homothety
  */
  void homothetie(const Vector& v);

  /*!
	\brief apply a rotation to the mesh
	\param p is axis of rotation
	\param angle is the angle of rotation
  */
  void rotate(const Vector& p, double angle);
  /*!
	  \brief apply a rotation around the X axis
	  \param angle is the angle of rotation
  */
  void rotateX(double angle);
  /*!
	\brief apply a rotation around the Y axis
	\param angle is the angle of rotation
  */
  void rotateY(double angle);
  /*!
	  \brief apply a rotation around the Z axis
	  \param angle is the angle of rotation
  */
  void rotateZ(double angle);
  /*!
	\brief apply a translation
	\param v vector of translation
  */
  void translation(const Vector& v);

  //deformations

  /*!
	  \brief applies to the mesh a deformation= in a sphere in the direction of a vector 
	  \param s the sphere
	  \param v the vector 
  */
  void sphereWarp(const Sphere& s, const Vector& v);
  /*!
	 \brief applies to the mesh a homothety in a sphere
	 \param s the sphere
	 \param v the vector of homothety
  */
  void sphereZoom(const Sphere& s, const Vector& v);
  /*!
   \brief applies to the mesh a rotation in a sphere
   \param s the sphere
   \param v the axis of rotation
   \param angle the angle of rotation
  */
  void sphereRotate(const Sphere& s, const Vector& v,double angle);






  void Load(const QString&);
  void SaveObj(const QString&, const QString&) const;
protected:
  void AddTriangle(int, int, int, int);
  void AddSmoothTriangle(int, int, int, int, int, int);
  void AddSmoothQuadrangle(int, int, int, int, int, int, int, int);
  void AddQuadrangle(int, int, int, int);
};

/*!
\brief Return the set of vertex indexes.
*/
inline std::vector<int> Mesh::VertexIndexes() const
{
  return varray;
}

/*!
\brief Return the set of normal indexes.
*/
inline std::vector<int> Mesh::NormalIndexes() const
{
  return narray;
}

/*!
\brief Get the vertex index of a given triangle.
\param t Triangle index.
\param i Vertex index.
*/
inline int Mesh::VertexIndex(int t, int i) const
{
  return varray.at(t * 3 + i);
}

/*!
\brief Get the normal index of a given triangle.
\param t Triangle index.
\param i Normal index.
*/
inline int Mesh::NormalIndex(int t, int i) const
{
  return narray.at(t * 3 + i);
}

/*!
\brief Get a triangle.
\param i Index.
\return The triangle.
*/
inline Triangle Mesh::GetTriangle(int i) const
{
  return Triangle(vertices.at(varray.at(i * 3 + 0)), vertices.at(varray.at(i * 3 + 1)), vertices.at(varray.at(i * 3 + 2)));
}

/*!
\brief Get a vertex.
\param i The index of the wanted vertex.
\return The wanted vertex (as a 3D Vector).
*/
inline Vector Mesh::Vertex(int i) const
{
  return vertices[i];
}

/*!
\brief Get a vertex from a specific triangle.
\param t The number of the triangle wich contain the wanted vertex.
\param v The triangle vertex: 0, 1, or 2.
\return The wanted vertex (as a 3D Vector).
*/
inline Vector Mesh::Vertex(int t, int v) const
{
  return vertices[varray[t * 3 + v]];
}

/*!
\brief Get the number of vertices in the geometry.
\return The number of vertices in the geometry, in other words the size of vertices.
*/
inline int Mesh::Vertexes() const
{
  return int(vertices.size());
}

/*!
\brief Get a normal.
\param i Index of the wanted normal.
\return The normal.
*/
inline Vector Mesh::Normal(int i) const
{
  return normals[i];
}

/*!
\brief Get the number of triangles.
*/
inline int Mesh::Triangles() const
{
  return int(varray.size()) / 3;
}

/*!
\brief Get a vertex.
\param i The index of the wanted vertex.
\return The wanted vertex (as a 3D Vector).
\see vertex(int i) const
*/
inline Vector Mesh::operator[](int i) const
{
  return vertices[i];
}

