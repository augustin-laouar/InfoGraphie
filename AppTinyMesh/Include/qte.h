#ifndef __Qte__
#define __Qte__

#include <QtWidgets/qmainwindow.h>
#include "realtime.h"
#include "meshcolor.h"

QT_BEGIN_NAMESPACE
	namespace Ui { class Assets; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
  Q_OBJECT
private:
  Ui::Assets* uiw;           //!< Interface

  MeshWidget* meshWidget;   //!< Viewer
  MeshColor meshColor;		//!< Mesh.

public:
  MainWindow();
  ~MainWindow();
  void CreateActions();
  void UpdateGeometry();

public slots:
  void editingSceneLeft(const Ray&);
  void editingSceneRight(const Ray&);
  //! Box mesh example
  void BoxMeshExample();
  //! Disk mesh example
  void DiskMeshExample();
  //! Cylinder mesh example
  void CylinderExample();
  //! Tore mesh example
  void ToreMeshExample();
  //! Sphere mesh example
  void SphereMeshExample();
  //! Sphere implicit mesh example
  void SphereImplicitExample();
  //! Capsule mesh example
  void CapsuleMeshExample();
  void ResetCamera();
  void UpdateMaterial();
  //! Example of merging 
  void MergeExample();
  //! Complex shape example (trying to creat a bird head)
  void ComplexShapeExample();
  //! Field mesh example. Create mountains randomly.
  void randomMountains();
  //! Field mesh example. Create holes randomly.
  void randomAbyss();
  //! Example of a field mesh after using a earthwork.
  void earthWork();
  //! Field mesh example. Create mountains randomly with a road crossing it.
  void roadAcrossMountains();
  //! Field mesh example. Create a moutain range randomly.
  void mountainRange();
  //! Load a field from an image given.
  void loadFielFromPath();
  //! Generate a field from parameters given.
  void fieldFromParameters();

};

#endif
