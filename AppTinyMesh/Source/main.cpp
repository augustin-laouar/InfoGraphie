#include "qte.h"
#include "mathematics.h"
#include <QtWidgets/qapplication.h>

int main(int argc, char *argv[])
{
	QApplication app(argc, argv);	
	MainWindow mainWin;
	mainWin.showMaximized();


	return app.exec();

}
