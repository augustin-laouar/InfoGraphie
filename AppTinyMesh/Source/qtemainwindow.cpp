#include "qte.h"
#include "implicits.h"
#include "ui_interface.h"
#include <chrono>

MainWindow::MainWindow() : QMainWindow(), uiw(new Ui::Assets)
{
	// Chargement de l'interface
    uiw->setupUi(this);

	// Chargement du GLWidget
	meshWidget = new MeshWidget;
	QGridLayout* GLlayout = new QGridLayout;
	GLlayout->addWidget(meshWidget, 0, 0);
	GLlayout->setContentsMargins(0, 0, 0, 0);
    uiw->widget_GL->setLayout(GLlayout);

	// Creation des connect
	CreateActions();

	meshWidget->SetCamera(Camera(Vector(10, 0, 0), Vector(0.0, 0.0, 0.0)));
}

MainWindow::~MainWindow()
{
	delete meshWidget;
}

void MainWindow::CreateActions()
{
	// Buttons
    connect(uiw->boxMesh, SIGNAL(clicked()), this, SLOT(BoxMeshExample()));
    connect(uiw->sphereImplicit, SIGNAL(clicked()), this, SLOT(SphereImplicitExample()));
	connect(uiw->diskMesh, SIGNAL(clicked()), this, SLOT(DiskMeshExample()));
	connect(uiw->sphereMesh, SIGNAL(clicked()), this, SLOT(SphereMeshExample()));
	connect(uiw->merge, SIGNAL(clicked()), this, SLOT(MergeExample()));
	connect(uiw->cylinderMesh, SIGNAL(clicked()), this, SLOT(CylinderExample()));
	connect(uiw->capsuleMesh, SIGNAL(clicked()), this, SLOT(CapsuleMeshExample()));
	connect(uiw->complexShap, SIGNAL(clicked()), this, SLOT(ComplexShapeExample()));
	connect(uiw->toreMesh, SIGNAL(clicked()), this, SLOT(ToreMeshExample()));
	connect(uiw->fieldMesh, SIGNAL(clicked()), this, SLOT(randomMountains()));
	connect(uiw->abyss_mesh, SIGNAL(clicked()), this, SLOT(randomAbyss()));
	connect(uiw->loadFromPath, SIGNAL(clicked()), this, SLOT(loadFielFromPath()));
	connect(uiw->ew_mesh, SIGNAL(clicked()), this, SLOT(earthWork()));
	connect(uiw->road_mesh, SIGNAL(clicked()), this, SLOT(roadAcrossMountains()));
	connect(uiw->mountain_mesh, SIGNAL(clicked()), this, SLOT(mountainRange()));
	connect(uiw->fieldFromParam, SIGNAL(clicked()), this, SLOT(fieldFromParameters()));
    connect(uiw->resetcameraButton, SIGNAL(clicked()), this, SLOT(ResetCamera()));
    connect(uiw->wireframe, SIGNAL(clicked()), this, SLOT(UpdateMaterial()));
    connect(uiw->radioShadingButton_1, SIGNAL(clicked()), this, SLOT(UpdateMaterial()));
    connect(uiw->radioShadingButton_2, SIGNAL(clicked()), this, SLOT(UpdateMaterial()));

	// Widget edition
	connect(meshWidget, SIGNAL(_signalEditSceneLeft(const Ray&)), this, SLOT(editingSceneLeft(const Ray&)));
	connect(meshWidget, SIGNAL(_signalEditSceneRight(const Ray&)), this, SLOT(editingSceneRight(const Ray&)));
}

void MainWindow::editingSceneLeft(const Ray&)
{
}

void MainWindow::editingSceneRight(const Ray&)
{
}


void MainWindow::DiskMeshExample()
{
	Mesh diskMesh = Mesh(Disk(Vector(0, 0, 0), 5));
	std::vector<Color> cols; //vecteurs des couleurs
	cols.resize(diskMesh.Vertexes()); //resize avec le nombre de sommet
	for (size_t i = 0; i < cols.size(); i++)
		cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 0.0);//calcul invraisemblable
	meshColor = MeshColor(diskMesh, cols, diskMesh.VertexIndexes());
	UpdateGeometry();
}

void MainWindow::SphereMeshExample() {
	Mesh sphereMesh = Mesh(Sphere(Vector(0,0,0), 3));
	std::vector<Color> cols; //vecteurs des couleurs
	cols.resize(sphereMesh.Vertexes()); //resize avec le nombre de sommet
	for (size_t i = 0; i < cols.size(); i++)
		cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 0.0);//calcul invraisemblable
	meshColor = MeshColor(sphereMesh, cols, sphereMesh.VertexIndexes());
	UpdateGeometry();
	
}


void MainWindow::CylinderExample() {
	Mesh cylindreMesh = Mesh(Cylinder(Vector(0, 0, 0), 3, 20));
	std::vector<Color> cols; //vecteurs des couleurs
	cols.resize(cylindreMesh.Vertexes()); //resize avec le nombre de sommet
	for (size_t i = 0; i < cols.size(); i++)
		cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 0.0);//calcul invraisemblable
	meshColor = MeshColor(cylindreMesh, cols, cylindreMesh.VertexIndexes());
	UpdateGeometry();


}
void MainWindow::ToreMeshExample() {

	Mesh toreMesh = Mesh(Tore(Vector(0,0,0),2,3));
	std::vector<Color> cols; //vecteurs des couleurs
	cols.resize(toreMesh.Vertexes()); //resize avec le nombre de sommet
	for (size_t i = 0; i < cols.size(); i++)
		cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378,1.0), 1.0, 0.0);//calcul invraisemblable
	meshColor = MeshColor(toreMesh, cols, toreMesh.VertexIndexes());
	UpdateGeometry();

}

void MainWindow::CapsuleMeshExample() {

	Mesh capsuleMesh = Mesh(Capsule(Vector(1,2,3), 4, 15));

	std::vector<Color> cols; //vecteurs des couleurs
	cols.resize(capsuleMesh.Vertexes()); //resize avec le nombre de sommet
	for (size_t i = 0; i < cols.size(); i++)
		cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 1.0, 0.0);//calcul invraisemblable
	meshColor = MeshColor(capsuleMesh, cols, capsuleMesh.VertexIndexes());
	UpdateGeometry();

}


void MainWindow::MergeExample() {
	Mesh sphereMesh = Mesh(Sphere(Vector(0, 0, 0), 3));
	Mesh toreMesh = Mesh(Tore(Vector(0, 0, 0), 3, 7));
	Mesh cylindreMesh = Mesh(Cylinder(Vector(15, -15, 0),3, 30));
	Mesh cylindreMesh2 = Mesh(Cylinder(Vector(-15, -15, 0), 3, 30));
	Mesh capsuleMesh = Mesh(Capsule(Vector(-15, -3, 0), 3, 10));
	Mesh capsuleMesh2 = Mesh(Capsule(Vector(15, -3, 0), 3, 10));

	capsuleMesh.rotateZ(PI / 2);
	capsuleMesh2.rotateZ(PI / 2);

	sphereMesh.merge(toreMesh);
	sphereMesh.merge(cylindreMesh);
	sphereMesh.merge(cylindreMesh2);
	sphereMesh.merge(capsuleMesh);
	sphereMesh.merge(capsuleMesh2);

	std::vector<Color> cols; //vecteurs des couleurs
	cols.resize(sphereMesh.Vertexes()); //resize avec le nombre de sommet
	for (size_t i = 0; i < cols.size(); i++)
		cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 1.0, 0.0);//calcul invraisemblable
	meshColor = MeshColor(sphereMesh, cols, sphereMesh.VertexIndexes());
	UpdateGeometry();


}

void MainWindow::ComplexShapeExample() {
	//(avant-arriere, gauche-droite, hauteur)
	Mesh complexShape = Mesh(Sphere(Vector(0, 0, 0), 3));
	complexShape.sphereWarp(Sphere(Vector(3, 0, 0), 1),Vector(2.5,0,0)); // creation du bec
	complexShape.sphereWarp(Sphere(Vector(1.5, 2, 1.5), 0.6), Vector(0, -0.8, -0.2)); // trou pour l'oeil 1
	complexShape.sphereWarp(Sphere(Vector(1.5, -2, 1.5), 0.6), Vector(0, 0.8, -0.2)); // trou pour l'oeil 2
	Mesh eye1(Sphere(Vector(1.5, 2, 1.5), 0.4));
	Mesh eye2(Sphere(Vector(1.5, -2, 1.5), 0.4));
	complexShape.merge(eye1);
	complexShape.merge(eye2);
	Mesh disk(Disk(Vector(0,0,2),3)); 
	disk.sphereWarp(Sphere(Vector(0, 0, 5), 1), Vector(0, 0, -1.5));
	disk.sphereWarp(Sphere(Vector(-2.18691, 0, 4.05364), 1), Vector(0.8, 0, -1.5));
	//disk.sphereWarp(Sphere(Vector(-2.62892, 0, 3.44526), 1), Vector(0, 0, -1.5));

	disk.sphereWarp(Sphere(Vector(2.18691, 0, 4.05364), 1), Vector(-0.8, 0, -1.5));
	//disk.sphereWarp(Sphere(Vector(2.62892, 0, 3.44526), 1), Vector(0, 0, -1.5));


	complexShape.merge(disk);

	std::vector<Color> cols; //vecteurs des couleurs
	cols.resize(complexShape.Vertexes()); //resize avec le nombre de sommet
	for (size_t i = 0; i < cols.size(); i++)
		cols[i] = Color(fmod(double(i) * 39.478378,3.0), fmod(double(i) * 39.478378, 3.0), fmod(double(i) * 39.478378,3.0), 0.0);//calcul invraisemblable
	meshColor = MeshColor(complexShape, cols, complexShape.VertexIndexes());
	UpdateGeometry();


}

void MainWindow::randomMountains() {
	int sizeX = 300;
	int sizeY = 300;
	QImage q(sizeX, sizeY,QImage::Format_Grayscale8);
	int color = 0;

	for (int i = 0; i < sizeX; i++) {
		for (int j = 0; j < sizeY; j++) { 
			q.setPixel(i,j, QColor(color,color,color ).rgb());
		}
	}
	HeightField hf(q, Vector(-sizeX/2, -sizeY/2, 0),Vector(sizeX/2 + 1,sizeY/2 +1,0), 1);
	hf.randomFieldInSquare(200, 1, 15, 35, 40, 50, 250, 50, 250); //genere 200 deformation aleatoirement
	Mesh hfMesh(hf);

	std::vector<Color> cols; //vecteurs des couleurs
	cols.resize(hfMesh.Vertexes()); //resize avec le nombre de sommet
	for (size_t i = 0; i < cols.size(); i++) {
		cols[i] = Color(fmod(double(i) * 39.478378, 3.0), fmod(double(i) * 39.478378, 3.0),1.0, 0.0);
	}
	meshColor = MeshColor(hfMesh, cols, hfMesh.VertexIndexes());
	UpdateGeometry();
	auto end = std::chrono::steady_clock::now();


}

void MainWindow::randomAbyss() {
	QImage q(300, 300, QImage::Format_Grayscale8);
	int color = 0;

	for (int i = 0; i < 300; i++) {
		for (int j = 0; j < 300; j++) {
			q.setPixel(i, j, QColor(color, color, color).rgb());
		}
	}
	HeightField hf(q, Vector(-150, -150, 0), Vector(149, 149, 0), 1);
	hf.randomFieldInSquare(200, -30, 15, 0, 40, 50, 250, 50, 250); //genere 200 deformation aleatoirement
	Mesh hfMesh(hf);

	std::vector<Color> cols; //vecteurs des couleurs
	cols.resize(hfMesh.Vertexes()); //resize avec le nombre de sommet
	for (size_t i = 0; i < cols.size(); i++) {
		cols[i] = Color(double(i * 2), fmod(double(i) * 39.478378, 2.0), 1.0, 0.0);
	}
	meshColor = MeshColor(hfMesh, cols, hfMesh.VertexIndexes());
	UpdateGeometry();


}



void MainWindow::earthWork() {
	QImage q(300, 300, QImage::Format_Grayscale8);
	int color = 0;

	for (int i = 0; i < 300; i++) {
		for (int j = 0; j < 300; j++) {
			q.setPixel(i, j, QColor(color, color, color).rgb());
		}
	}
	HeightField hf(q, Vector(-150, -150, 0), Vector(149, 149, 0), 1);
	hf.randomFieldInSquare(200, 1, 15, 35, 40, 50, 250, 50, 250); //genere 200 deformation aleatoirement
	hf.brush(Vector(0,0,0),70);
	Mesh hfMesh(hf);

	std::vector<Color> cols; //vecteurs des couleurs
	cols.resize(hfMesh.Vertexes()); //resize avec le nombre de sommet
	for (size_t i = 0; i < cols.size(); i++) {
		cols[i] = Color(double(i * 2), fmod(double(i) * 39.478378, 2.0), 1.0, 0.0);
	}
	meshColor = MeshColor(hfMesh, cols, hfMesh.VertexIndexes());
	UpdateGeometry();
}


void MainWindow::roadAcrossMountains() {
	QImage q(300, 300, QImage::Format_Grayscale8);
	int color = 0;

	for (int i = 0; i < 300; i++) {
		for (int j = 0; j < 300; j++) {
			q.setPixel(i, j, QColor(color, color, color).rgb());
		}
	}
	HeightField hf(q, Vector(-150, -150, 0), Vector(149, 149, 0), 1);
	hf.randomFieldInSquare(200, 1, 15, 35, 40, 50, 250, 50, 250); //genere 200 deformation aleatoirement
	for (int i = -150; i < 150; i += 10) {
		hf.brush(Vector(i, 0, 0), 10);
	}
	Mesh hfMesh(hf);

	std::vector<Color> cols; //vecteurs des couleurs
	cols.resize(hfMesh.Vertexes()); //resize avec le nombre de sommet
	for (size_t i = 0; i < cols.size(); i++) {
		cols[i] = Color(double(i * 2), fmod(double(i) * 39.478378, 2.0), 1.0, 0.0);
	}
	meshColor = MeshColor(hfMesh, cols, hfMesh.VertexIndexes());
	UpdateGeometry();
}


void MainWindow::mountainRange() {
	QImage q(300, 300, QImage::Format_Grayscale8);
	int color = 0;

	for (int i = 0; i < 300; i++) {
		for (int j = 0; j < 300; j++) {
			q.setPixel(i, j, QColor(color, color, color).rgb());
		}
	}
	HeightField hf(q, Vector(-150, -150, 0), Vector(149, 149, 0), 1);
	hf.randomFieldInSquare(200, 1, 15, 25, 40, 0, 300, 100, 180); //genere 200 deformation aleatoirement
	Mesh hfMesh(hf);

	std::vector<Color> cols; //vecteurs des couleurs
	cols.resize(hfMesh.Vertexes()); //resize avec le nombre de sommet
	for (size_t i = 0; i < cols.size(); i++) {
		cols[i] = Color(double(i * 2), fmod(double(i) * 39.478378, 2.0), 1.0, 0.0);
	}
	meshColor = MeshColor(hfMesh, cols, hfMesh.VertexIndexes());
	UpdateGeometry();
}




void MainWindow::loadFielFromPath() {
	QImage img = QImage(uiw->imagePath->text());
	if (img.isNull()) {
		uiw->infos->setText("Loading image failed. Verify the path is correct.");
		return;
	}
	uiw->infos->setText(uiw->imagePath->text() + " loaded");
	int sizeX = img.width();
	int sizeY = img.height();
	HeightField hf(img, Vector(-sizeX / 2, -sizeY / 2, 0), Vector(sizeX / 2 + 1, sizeY / 2 + 1, 0), 1);
	Mesh hfMesh(hf);

	std::vector<Color> cols; //vecteurs des couleurs
	cols.resize(hfMesh.Vertexes()); //resize avec le nombre de sommet
	for (size_t i = 0; i < cols.size(); i++) {
		cols[i] = Color(double(i * 2), fmod(double(i) * 39.478378, 2.0), 1.0, 0.0);
	}
	meshColor = MeshColor(hfMesh, cols, hfMesh.VertexIndexes());
	UpdateGeometry();

}

void MainWindow::fieldFromParameters() {
	int imageW = uiw->width->text().toInt();
	int imageH = uiw->width->text().toInt();
	Vector a(uiw->xA->text().toDouble(), uiw->yA->text().toDouble(), uiw->zA->text().toDouble());
	Vector b(uiw->xB->text().toDouble(), uiw->yB->text().toDouble(), uiw->zB->text().toDouble());
	double zScaling = uiw->zScaling->text().toDouble();
	int nbDeform = uiw->nbDeformations->text().toInt();
	double dhMin = uiw->minDH->text().toDouble();
	double dhMax = uiw->maxDH->text().toDouble();
	double drMin = uiw->minDR->text().toDouble();
	double drMax = uiw->maxDR->text().toDouble();
	int minX = uiw->dfXmin->text().toDouble();
	int maxX = uiw->dfXmax->text().toDouble();
	int minY = uiw->dfYmin->text().toDouble();
	int maxY = uiw->dfYmax->text().toDouble();

	QImage q(imageW, imageH, QImage::Format_Grayscale8);
	int color = 0;

	for (int i = 0; i < imageW; i++) {
		for (int j = 0; j < imageH; j++) {
			q.setPixel(i, j, QColor(color, color, color).rgb());
		}
	}
	HeightField hf(q, a, b, zScaling);
	if (nbDeform != 0 && ((dhMin == 0 && drMax == 0) || (drMin == 0 && drMax == 0))) {
		uiw->infos->setText("Generation failed. Parameters missing for deformations.");
		return;
	}
	if (nbDeform != 0 && ((minX == 0 && maxX == 0) || (minY == 0 && maxY == 0))) {
		uiw->infos->setText("Error : Range for deformations is 0 on axis X or Y.");
		return;
	}
	if (nbDeform != 0 && ((minX > maxX) || (minY>maxY))) {
		uiw->infos->setText("Error : range axis X or range axis Y.");
		return;
	}
	if (nbDeform != 0 && dhMin > dhMax) {
		uiw->infos->setText("Error : minimum height > maximum height.");
		return;
	}
	if (nbDeform != 0 && drMin > drMax) {
		uiw->infos->setText("Error : minimum radius > maximum radius.");
		return;
	}

	hf.randomFieldInSquare(nbDeform, dhMin, drMin, dhMax, drMax, minX, maxX, minY, maxY); //genere 200 deformation aleatoirement
	Mesh hfMesh(hf);
	uiw->infos->setText("Generation successfull.");
	std::vector<Color> cols; //vecteurs des couleurs
	cols.resize(hfMesh.Vertexes()); //resize avec le nombre de sommet
	for (size_t i = 0; i < cols.size(); i++) {
		cols[i] = Color(double(i * 2), fmod(double(i) * 39.478378, 2.0), 1.0, 0.0);
	}
	meshColor = MeshColor(hfMesh, cols, hfMesh.VertexIndexes());
	UpdateGeometry();

}


void MainWindow::BoxMeshExample()
{
	Mesh boxMesh = Mesh(Box(1.0));
	Sphere s(Vector(0.5, 0, 0), 1);
	boxMesh.sphereWarp(s, Vector(2, 2, 2));
	std::vector<Color> cols; //vecteurs des couleurs
	cols.resize(boxMesh.Vertexes()); //resize avec le nombre de sommet
	for (size_t i = 0; i < cols.size(); i++)
		cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 0.0);//calcul invraisemblable
	meshColor = MeshColor(boxMesh, cols, boxMesh.VertexIndexes());
	UpdateGeometry();
}


void MainWindow::SphereImplicitExample()
{
  AnalyticScalarField implicit;

  Mesh implicitMesh;
  implicit.Polygonize(31, implicitMesh, Box(2.0));
  std::vector<Color> cols;
  cols.resize(implicitMesh.Vertexes());
  for (size_t i = 0; i < cols.size(); i++)
    cols[i] = Color(0.8, 0.8, 0.8);

  meshColor = MeshColor(implicitMesh, cols, implicitMesh.VertexIndexes());
  UpdateGeometry();
}


void MainWindow::UpdateGeometry()
{
	meshWidget->ClearAll();
	meshWidget->AddMesh("BoxMesh", meshColor);

    uiw->lineEdit->setText(QString::number(meshColor.Vertexes()));
    uiw->lineEdit_2->setText(QString::number(meshColor.Triangles()));

	UpdateMaterial();
}

void MainWindow::UpdateMaterial()
{
    meshWidget->UseWireframeGlobal(uiw->wireframe->isChecked());

    if (uiw->radioShadingButton_1->isChecked())
		meshWidget->SetMaterialGlobal(MeshMaterial::Normal);
	else
		meshWidget->SetMaterialGlobal(MeshMaterial::Color);
}

void MainWindow::ResetCamera()
{
	meshWidget->SetCamera(Camera(Vector(-10.0), Vector(0.0)));
}
