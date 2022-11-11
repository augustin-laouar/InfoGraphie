#include "basic_shapes.h"
HeightField::HeightField(const QImage& qi, const Vector& a, const Vector& b, double scaleZ) {
	this->a = a;
	this->b = b;
	nx = qi.width();
	ny = qi.height();
	scaleX = ((b[0] - a[0]) / (nx - 1));
	scaleY = ((b[1] - a[1]) / (ny - 1));
	z.reserve(nx * ny); // autant de point que de pixels dans l'image
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			QRgb pix = qi.pixel(i, j);
			z.push_back(qGray(pix) * scaleZ);
		}
	}
}

void HeightField::createRelief(const Vector& c, double r, double h) {
	Vector elevator(0, 0, h);
	Vector curr;
	double dist;
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			curr = getPoint(i, j);
			dist = sqrt((c[0] - curr[0]) * (c[0] - curr[0]) + (c[1] - curr[1]) * (c[1] - curr[1]) + (c[2] - curr[2]) * (c[2] - curr[2]));
			if (dist <= r) {
				double dr = 1 - (dist / r);
				if (dr != 0) { //smooth stepping
					dr = 3 * dr * dr - 2 * dr * dr * dr;
				}
				z[getIndex(i, j)] += h * dr;
			}
		}
	}
}

void HeightField::bump(const Vector& c, double d, double h) {
	double d2D;
	Vector curr;
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			curr = getPoint(i, j);
			d2D = sqrt((c[0] - curr[0]) * (c[0] - curr[0]) + (c[1] - curr[1]) * (c[1] - curr[1]));
			if (d2D <= d) {
				z[getIndex(i, j)] += h;
			}
		}
	}
}

void HeightField::terrasser(const Vector& c, double r, double r2, double h) {
	Vector curr;
	double dist2D;
	double ratio;
	double diff;
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			curr = getPoint(i, j);
			dist2D = sqrt((c[0] - curr[0]) * (c[0] - curr[0]) + (c[1] - curr[1]) * (c[1] - curr[1]));
			if (dist2D <= r) { // dans le rayon d'activation
				diff = curr[2] - h;
				if (diff > 0) {
					z[getIndex(i, j)] = curr[2] - diff;
				}
			}

			if (dist2D > r && dist2D <= r2) {
				ratio = 1 - (dist2D - r) / (r2 - r);
				if (ratio != 0) { //smooth stepping
					ratio = 3 * ratio * ratio - 2 * ratio * ratio * ratio;
				}
				diff = curr[2] - h;
				if (diff > 0) {
					z[getIndex(i, j)] = curr[2] - diff * ratio;
				}
			}

		}
	}
}


void HeightField::randomField(int nb, double minH, double minR, double maxH, double maxR) {
	int i, j;
	double d, h;
	std::random_device rd;
	std::default_random_engine eng(rd());
	std::uniform_real_distribution<double> distrH(minH, maxH);
	std::uniform_real_distribution<double> distrD(minR, maxR);
	for (int i2 = 0; i2 < nb; i2++) {
		i = rand() % nx;
		j = rand() % ny;
		h = distrH(eng);
		d = distrD(eng);
		createRelief(getPoint(i, j), d, h);
	}
}

void HeightField::randomFieldInSquare(int nb, double minH, double minD, double maxH, double maxD, int minX, int maxX, int minY, int maxY) {
	int i, j;
	double d, h;
	std::random_device rd;
	std::default_random_engine eng(rd());
	std::uniform_real_distribution<double> distrH(minH, maxH);
	std::uniform_real_distribution<double> distrD(minD, maxD);
	for (int i2 = 0; i2 < nb; i2++) {
		i = rand() % (maxX - minX) + minX;
		j = rand() % (maxY - minY) + minY;
		h = distrH(eng);
		d = distrD(eng);
		createRelief(getPoint(i, j), d, h);
	}
}