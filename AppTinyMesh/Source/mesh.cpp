#include "mesh.h"
/*!
\class Mesh mesh.h

\brief Core triangle mesh class.
*/



/*!
\brief Initialize the mesh to empty.
*/
Mesh::Mesh()
{
}

/*!
\brief Initialize the mesh from a list of vertices and a list of triangles.

Indices must have a size multiple of three (three for triangle vertices and three for triangle normals).

\param vertices List of geometry vertices.
\param indices List of indices wich represent the geometry triangles.
*/
Mesh::Mesh(const std::vector<Vector>& vertices, const std::vector<int>& indices) :vertices(vertices), varray(indices)
{
  normals.resize(vertices.size(), Vector::Z);
}

/*!
\brief Create the mesh.

\param vertices Array of vertices.
\param normals Array of normals.
\param va, na Array of vertex and normal indexes.
*/
Mesh::Mesh(const std::vector<Vector>& vertices, const std::vector<Vector>& normals, const std::vector<int>& va, const std::vector<int>& na) :vertices(vertices), normals(normals), varray(va), narray(na)
{
}

/*!
\brief Reserve memory for arrays.
\param nv,nn,nvi,nvn Number of vertices, normals, vertex indexes and vertex normals.
*/
void Mesh::Reserve(int nv, int nn, int nvi, int nvn)
{
  vertices.reserve(nv);
  normals.reserve(nn);
  varray.reserve(nvi);
  narray.reserve(nvn);
}

/*!
\brief Empty
*/
Mesh::~Mesh()
{
}

/*!
\brief Smooth the normals of the mesh.

This function weights the normals of the faces by their corresponding area.
\sa Triangle::AreaNormal()
*/
void Mesh::SmoothNormals()
{
  // Initialize 
  normals.resize(vertices.size(), Vector::Null);

  narray = varray;

  // Accumulate normals
  for (int i = 0; i < varray.size(); i += 3)
  {
    Vector tn = Triangle(vertices[varray.at(i)], vertices[varray.at(i + 1)], vertices[varray.at(i + 2)]).AreaNormal();
    normals[narray[i + 0]] += tn;
    normals[narray[i + 1]] += tn;
    normals[narray[i + 2]] += tn;
  }

  // Normalize 
  for (int i = 0; i < normals.size(); i++)
  {
    Normalize(normals[i]);
  }
}

/*!
\brief Add a smooth triangle to the geometry.
\param a, b, c Index of the vertices.
\param na, nb, nc Index of the normals.
*/
void Mesh::AddSmoothTriangle(int a, int na, int b, int nb, int c, int nc)
{
  varray.push_back(a);
  narray.push_back(na);
  varray.push_back(b);
  narray.push_back(nb);
  varray.push_back(c);
  narray.push_back(nc);
}

/*!
\brief Add a triangle to the geometry.
\param a, b, c Index of the vertices.
\param n Index of the normal.
*/
void Mesh::AddTriangle(int a, int b, int c, int n)
{
  varray.push_back(a);
  narray.push_back(n);
  varray.push_back(b);
  narray.push_back(n);
  varray.push_back(c);
  narray.push_back(n);
}



/*!
\brief Add a smmoth quadrangle to the geometry.

Creates two smooth triangles abc and acd.

\param a, b, c, d  Index of the vertices.
\param na, nb, nc, nd Index of the normal for all vertices.
*/
void Mesh::AddSmoothQuadrangle(int a, int na, int b, int nb, int c, int nc, int d, int nd)
{
  // First triangle
  AddSmoothTriangle(a, na, b, nb, c, nc);

  // Second triangle
  AddSmoothTriangle(a, na, c, nc, d, nd);
}

/*!
\brief Add a quadrangle to the geometry.

\param a, b, c, d  Index of the vertices and normals.
*/
void Mesh::AddQuadrangle(int a, int b, int c, int d)
{
  AddSmoothQuadrangle(a, a, b, b, c, c, d, d);
}

/*!
\brief Compute the bounding box of the object.
*/
Box Mesh::GetBox() const
{
  if (vertices.size() == 0)
  {
    return Box::Null;
  }
  return Box(vertices);
}

/*!
\brief Creates an axis aligned box.

The object has 8 vertices, 6 normals and 12 triangles.
\param box The box.
*/
Mesh::Mesh(const Box& box)
{
  // Vertices
  vertices.resize(8);

  for (int i = 0; i < 8; i++)
  {
    vertices[i] = box.Vertex(i);
  }

  // Normals
  normals.push_back(Vector(-1, 0, 0));
  normals.push_back(Vector(1, 0, 0));
  normals.push_back(Vector(0, -1, 0));
  normals.push_back(Vector(0, 1, 0));
  normals.push_back(Vector(0, 0, -1));
  normals.push_back(Vector(0, 0, 1));

  // Reserve space for the triangle array
  varray.reserve(12 * 3);
  narray.reserve(12 * 3);

  AddTriangle(0, 2, 1, 4);
  AddTriangle(1, 2, 3, 4);

  AddTriangle(4, 5, 6, 5);
  AddTriangle(5, 7, 6, 5);

  AddTriangle(0, 4, 2, 0);
  AddTriangle(4, 6, 2, 0);

  AddTriangle(1, 3, 5, 1);
  AddTriangle(3, 7, 5, 1);

  AddTriangle(0, 1, 5, 2);
  AddTriangle(0, 5, 4, 2);

  AddTriangle(3, 2, 7, 3);
  AddTriangle(6, 7, 2, 3);
}

Mesh::Mesh(const Disk& disk) {
    const int div = 100;
    double step = 2 * PI / div;
    double alpha;

    vertices.resize(1);
    vertices[0] = disk.getCenter();


    normals.push_back(Vector(0, 1, 0));
    //generation des points du disque
    for (int i = 0; i <= div; ++i) {
        alpha = i * step;
        Vector v(cos(alpha) * disk.getRadius(), 0, sin(alpha) * disk.getRadius());
        v += disk.getCenter();
        vertices.push_back(v);
    }

    //création des triangles
    varray.reserve(div * 3);
    narray.reserve(div * 3);

    for (int i = 1; i < div; i++) {
        AddTriangle(0,i,i+1,0);
    }
    AddTriangle(0, div, 1, 0);
}


Mesh::Mesh(const Cylinder& c) {

    const int div =50;
    double step = 2 * PI / div;
    double alpha;

    // creation du disque du bas


    vertices.push_back(c.getCenter() );


    normals.push_back(Vector(cos(alpha), 0, sin(alpha)));
    //generation des points du disque
    for (int i = 0; i <= div; ++i) {
        alpha = i * step;
        //on considere le centre de ce disque  comme étant celui du cylindre
        vertices.push_back(Vector(cos(alpha) * c.getRadius(), 0, sin(alpha) * c.getRadius()) + c.getCenter());
    }

    //création des triangles
    varray.reserve(div * 3 * 2);
    narray.reserve(div * 3 * 2);

        for (int i = 1; i < div; i++) {
            AddTriangle(0, i, i + 1, 0);
        }
        AddTriangle(0, div, 1, 0);
    
        
   
        Vector center(0, c.getHeight(), 0); // la coordonnee x du centre du disque du bas = celle du centre du disque du bas + la hauteur 
        center += c.getCenter();
        vertices.push_back(center);

        //creation du disque du haut
        normals.push_back(Vector(cos(alpha),1, sin(alpha)));
        
        //generation des points du disque
        for (int i = 0; i <= div; ++i) {
            alpha = i * step;
            
            vertices.push_back(Vector(cos(alpha) * c.getRadius(), 0, sin(alpha) * c.getRadius()) + center);
        }

        //création des triangle

        for (int i = 1; i < div; i++) {
            AddTriangle(div + 2, div + 2 + i, div + 2 + i + 1, 1);
        }
        AddTriangle(div + 2,div + 2 + div, div + 2 + 1, 1);

        // a ce stade on a créé deux disques séparé avec la hauteur du cylindre

        for (int i = 1; i <div ; i++)
        {
            AddTriangle(i, div + 3 + i, i + 1, 1);
            AddTriangle(i + 1, div + 3 + i, div+3+i+1, 1);

    }
    AddTriangle(div, div + div + 3, 1, 1);
    AddTriangle(1, div + div + 3, div + 3 + 1, 1);

      
}



Mesh::Mesh(const Sphere& sphere) {
    const int divBeta = 64;
    const int divAlpha = divBeta;
    double beta, alpha;


    varray.reserve(divBeta * divBeta);
    narray.reserve((divBeta * divBeta));
    int counter = 0;
    int normaleCounter = 0;
    int previousCircle[divBeta+1];//permet de stocker les indices des sommets du cercle précent
    int buffer[divBeta+1];

    for (int i = 0; i <= divAlpha; ++i) {
        alpha = -0.5 * PI + double(i) * PI / divAlpha;
        for (int j = 0; j <= divBeta; ++j)
        {
            beta = double(j) * 2.0 * PI / (divBeta);
            Vector v(cos(alpha) * cos(beta) * sphere.getRadius(), sin(alpha) * sphere.getRadius(), cos(alpha) * sin(beta) * sphere.getRadius());
            v += sphere.getCenter();
            vertices.push_back(v);
            buffer[j] = counter;
            counter++;
            if (i > 0) { // ce n'est pas le premier cercle
                if (j != 0) {//pas premier point de ce cercle
                    Vector normale = vertices[counter-1] - vertices[counter-2];
                    normale *= (vertices[counter - 2] - vertices[previousCircle[j-1]]);
                    normals.push_back(normale);
                    normaleCounter++;
                    AddTriangle(previousCircle[j - 1], counter - 2,counter -1, normaleCounter - 1);//triangle 1
                    normale = vertices[counter - 1] - vertices[previousCircle[j-1]];
                    normale *= (vertices[previousCircle[j]] - vertices[counter-1]);
                    normals.push_back(normale);
                    normaleCounter++;
                    AddTriangle(previousCircle[j], counter - 1, previousCircle[j - 1], normaleCounter - 1);//triangle 2
                    if (j == divBeta) {//relier le dernier point avec le premier 
                        Vector normale = vertices[counter - 1] - vertices[counter - divBeta - 1];
                        normale *= (vertices[counter - 1] - vertices[previousCircle[divBeta]]);
                        normals.push_back(normale);
                        normaleCounter++;
                        AddTriangle(counter - 1, previousCircle[divBeta], counter - divBeta - 1, normaleCounter - 1);
                        Vector normale2 = vertices[previousCircle[divBeta]] - vertices[counter - 1];
                        normale2 *= (vertices[counter - 1] - vertices[counter - divBeta - 1]);
                        normals.push_back(normale);
                        normaleCounter++;
                        AddTriangle(previousCircle[divBeta], counter - 1, counter - divBeta - 1, normaleCounter - 1);
                    }
                }
                
            }
        }
        for (int k = 0; k <= divBeta; k++){
            previousCircle[k] = buffer[k];
        }
    }
}

int wrapping(int x, int step) {
    return (x + 1) % step;
}
Mesh::Mesh(const Tore& tore) {
    int innerStep = 48;
    int outerStep = 48;
    double outerR;
    double innerR;
    double innerDistance;
    double x, y, z;

    //generation des sommets
    for (int i = 0; i < outerStep; i++) {
         outerR = double(i) / outerStep * PI *2.0;
        for (int j = 0; j < innerStep; j++) {
             innerR = double(j) / innerStep * PI * 2.0;
             innerDistance = cos(innerR) * tore.getInnerRadius();
             x = sin(outerR) * (innerDistance + tore.getOuterRadius());
             y = cos(outerR) * (innerDistance + tore.getOuterRadius());
             z = sin(innerR) * tore.getInnerRadius();
             vertices.push_back(Vector(x, y, z) + tore.getCenter());
        }
    }
    //generation des triangles
    varray.reserve(innerStep * outerStep * 6);
    narray.reserve(innerStep * outerStep * 6);
    int p1, p2, p3;
    int normaleCounter = 0;
    for (int i = 0; i < outerStep; i++) {
        for (int j = 0; j < innerStep; j++) {
            //triangle 1
            p1 = (i * innerStep) + j;
            p2 = (i * innerStep) + wrapping(j,innerStep);
            p3 = wrapping(i, outerStep) * innerStep + wrapping(j, innerStep);
            Vector normale = vertices[p1] - vertices[p2];
            normale *= Vector(vertices[p1] - vertices[p3]);
            normals.push_back(normale);
            AddTriangle(p1, p2, p3, normaleCounter);
            normaleCounter++;

            //triangle 2
            p1 = wrapping(i,outerStep)*innerStep + wrapping(j,innerStep);
            p2 = wrapping(i,outerStep) * innerStep + j;
            p3 = i * innerStep + j;
            normale = vertices[p1] - vertices[p2];
            normale *= Vector(vertices[p1] - vertices[p3]);
            normals.push_back(normale);
            AddTriangle(p1, p2, p3, normaleCounter);
            normaleCounter++;


        }
    }


}

Mesh::Mesh(const Capsule& capsule) {
   
    
    // on considere le centre de la capsule comme celui du cylindre ainsi que son rayon et sa hauteur/2
    Mesh cylindre(Cylinder(capsule.getCenter(), capsule.getRadius(), capsule.getHeight()- 2 * capsule.getRadius()));
    
    Vector temp(0, capsule.getHeight()- 2* capsule.getRadius(), 0);

    Mesh sphere1 = Mesh(Sphere(capsule.getCenter() , capsule.getRadius()));
    Mesh sphere2 = Mesh(Sphere(capsule.getCenter() + temp, capsule.getRadius()));
  
    cylindre.merge(sphere1);
    cylindre.merge(sphere2);
    this->merge(cylindre);
    
    

}

Mesh::Mesh(const HeightField& hf) {
    for (int i = 0; i < hf.getNx(); i++) {
        for (int j = 0; j < hf.getNy(); j++) {
            vertices.push_back(hf.getPoint(i, j));
        }
    }
    //varray.reserve(3);
    //narray.reserve(hf.getNx() * hf.getNy());
    int v1, v2, v3;
    Vector normale;
    for (int i = 1; i < hf.getNx(); i++) {
        for (int j = 1; j < hf.getNy(); j++) {
            v1 = (i - 1) * hf.getNy() + (j - 1);
            v2 = (i - 1) * hf.getNy() + j;
            v3 = i * hf.getNy() + j;
            normale = vertices[v1] - vertices[v2];
            normale *= Vector(vertices[v1] - vertices[v3]);
            normals.push_back(normale);
            AddTriangle(v1, v2, v3, normals.size() - 1);
            v2 = i * hf.getNy() + (j - 1);
            normale = vertices[v1] - vertices[v2];
            normale *= Vector(vertices[v1] - vertices[v3]);
            normals.push_back(normale);
            AddTriangle(v1, v2, v3, normals.size() - 1);
        }
    }
}



/*!
\brief Scale the mesh.
\param s Scaling factor.
*/
void Mesh::Scale(double s)
{
    // Vertexes
    for (int i = 0; i < vertices.size(); i++)
    {
        vertices[i] *= s;
    }

    if (s < 0.0)
    {
        // Normals
        for (int i = 0; i < normals.size(); i++)
        {
            normals[i] = -normals[i];
        }
    }
}



#include <QtCore/QFile>
#include <QtCore/QTextStream>
#include <QtCore/QRegularExpression>
#include <QtCore/qstring.h>

/*!
\brief Import a mesh from an .obj file.
\param filename File name.
*/
void Mesh::Load(const QString& filename)
{
  vertices.clear();
  normals.clear();
  varray.clear();
  narray.clear();

  QFile data(filename);

  if (!data.open(QFile::ReadOnly))
    return;
  QTextStream in(&data);

  // Set of regular expressions : Vertex, Normal, Triangle
  QRegularExpression rexv("v\\s*([-|+|\\s]\\d*\\.\\d+)\\s*([-|+|\\s]\\d*\\.\\d+)\\s*([-|+|\\s]\\d*\\.\\d+)");
  QRegularExpression rexn("vn\\s*([-|+|\\s]\\d*\\.\\d+)\\s*([-|+|\\s]\\d*\\.\\d+)\\s*([-|+|\\s]\\d*\\.\\d+)");
  QRegularExpression rext("f\\s*(\\d*)/\\d*/(\\d*)\\s*(\\d*)/\\d*/(\\d*)\\s*(\\d*)/\\d*/(\\d*)");
  while (!in.atEnd())
  {
    QString line = in.readLine();
    QRegularExpressionMatch match = rexv.match(line);
    QRegularExpressionMatch matchN = rexn.match(line);
    QRegularExpressionMatch matchT = rext.match(line);
    if (match.hasMatch())//rexv.indexIn(line, 0) > -1)
    {
      Vector q = Vector(match.captured(1).toDouble(), match.captured(2).toDouble(), match.captured(3).toDouble()); vertices.push_back(q);
    }
    else if (matchN.hasMatch())//rexn.indexIn(line, 0) > -1)
    {
      Vector q = Vector(matchN.captured(1).toDouble(), matchN.captured(2).toDouble(), matchN.captured(3).toDouble());  normals.push_back(q);
    }
    else if (matchT.hasMatch())//rext.indexIn(line, 0) > -1)
    {
      varray.push_back(matchT.captured(1).toInt() - 1);
      varray.push_back(matchT.captured(3).toInt() - 1);
      varray.push_back(matchT.captured(5).toInt() - 1);
      narray.push_back(matchT.captured(2).toInt() - 1);
      narray.push_back(matchT.captured(4).toInt() - 1);
      narray.push_back(matchT.captured(6).toInt() - 1);
    }
  }
  data.close();
}

/*!
\brief Save the mesh in .obj format, with vertices and normals.
\param url Filename.
\param meshName %Mesh name in .obj file.
*/
void Mesh::SaveObj(const QString& url, const QString& meshName) const
{
  QFile data(url);
  if (!data.open(QFile::WriteOnly))
    return;
  QTextStream out(&data);
  out << "g " << meshName << Qt::endl;
  for (int i = 0; i < vertices.size(); i++)
    out << "v " << vertices.at(i)[0] << " " << vertices.at(i)[1] << " " << vertices.at(i)[2] << QString('\n');
  for (int i = 0; i < normals.size(); i++)
    out << "vn " << normals.at(i)[0] << " " << normals.at(i)[1] << " " << normals.at(i)[2] << QString('\n');
  for (int i = 0; i < varray.size(); i += 3)
  {
    out << "f " << varray.at(i) + 1 << "//" << narray.at(i) + 1 << " "
      << varray.at(i + 1) + 1 << "//" << narray.at(i + 1) + 1 << " "
      << varray.at(i + 2) + 1 << "//" << narray.at(i + 2) + 1 << " "
      << "\n";
  }
  out.flush();
  data.close();
}



void Mesh::merge(const Mesh& a)
{
    int size = this->vertices.size();

    for (int i = 0; i < a.vertices.size(); i++)
    {
        this->vertices.push_back(a.vertices[i]);

    }

    for (int j = 0; j < a.varray.size(); j++)
    {
        this->varray.push_back(a.varray[j] + size);
    }
    size = this->normals.size();
    for (int k = 0; k < a.normals.size(); k++)
    {
        this->normals.push_back(a.normals[k]);
    }
    
    for (int m = 0; m < a.narray.size(); m++)
    {
        this->narray.push_back(a.narray[m]+size);
    }
}
void Mesh::transformation(const Matrix& m) {
    for (int i = 0; i < vertices.size(); i++) {
        vertices[i] = m * vertices[i];
    }
}
void Mesh::homothetie(const Vector& v) {
    Matrix m;
    m.homothetieMatrix(v);
    transformation(m);

}
void Mesh::rotate(const Vector& v, double angle) {
    Matrix m;
    m.rotateMatrix(v, angle);
    std::cout << m << std::endl;
    transformation(m);
}

void Mesh::rotateX(double angle) {
    rotate(Vector(1, 0, 0), angle);
}
void Mesh::rotateY(double angle) {
    rotate(Vector(0, 1, 0), angle);
}
void Mesh::rotateZ(double angle) {
    rotate(Vector(0, 0, 1), angle);
}

void Mesh::translation(const Vector& v) {
    for (int i = 0; i < vertices.size(); i++) {
        vertices[i] = vertices[i] + v;
    }
}

void Mesh::sphereWarp(const Sphere& s,const Vector& v) {
    double d;
    for (int i = 0; i < vertices.size(); i++) {
        d = Distance(s.getCenter(), vertices[i]);
        if (d < s.getRadius()) {
            double dr = 1 - (d / s.getRadius());
            if (dr != 0) { //smooth stepping
                dr = 3 * dr * dr - 2 * dr * dr * dr;
            }
            Vector deplacement(v * dr);
            vertices[i] += deplacement;
        }
    }
}
void Mesh::sphereZoom(const Sphere& s, const Vector& v) {
    double d;
    for (int i = 0; i < vertices.size(); i++) {
        d = Distance(s.getCenter(), vertices[i]);
        if (d < s.getRadius()) {
            double dr = 1 - (d / s.getRadius());
            if (dr != 0) { //smooth stepping
                dr = 3 * dr * dr - 2 * dr * dr * dr;
            }
            Vector zoom(v * dr);
            Matrix m;
            m.homothetieMatrix(zoom);
            vertices[i] = m * vertices[i];
        }
    }
}
void Mesh::sphereRotate(const Sphere& s, const Vector& v, double angle) {
    double d;
    for (int i = 0; i < vertices.size(); i++) {
        d = Distance(s.getCenter(), vertices[i]);
        if (d < s.getRadius()) {
            double dr = 1 - (d / s.getRadius());
            if (dr != 0) { //smooth stepping
                dr = 3 * dr * dr - 2 * dr * dr * dr;
            }
            Matrix m;
            m.rotateMatrix(v,angle * dr);
            vertices[i] = m * vertices[i];
        }
    }
}

