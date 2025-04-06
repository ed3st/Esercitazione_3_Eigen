#include <iostream>
#include <Eigen/Eigen>
#include <Eigen/LU>
#include <iomanip>

using namespace std;
using namespace Eigen;

/*Calcolo dell'errore relativo tr due vettori in norma 2.*/
double relative_error(VectorXd y, VectorXd z) { 
	double err = (z-y).norm()/y.norm();
	return err;
};

/*Risoluzione di  un sistema lineare(Ax = b) tramite la decomposizione LU con pivoting parziale + calcolo dell'errere relativo.
La funzione solve_system_PALU prende in input la matrice A, il vettore termine noto b e la soluzione esatta del sistema x.*/
int solve_system_PALU(MatrixXd A, VectorXd b, VectorXd x) { 
	if ( abs(A.determinant()) < 1e-15 ) { //Controllo che la matrice A sia invertibile.
		cerr << "A is a non-invertible matrix, PALU decomposition is not executable." << endl;
		return 1;
	}
	else {
		Vector2d xLu = A.partialPivLu().solve(b); //Risolvo il sistema.
		cout << "Sol with PALU decomposition:\n" << setprecision(15) << scientific << xLu << endl; //Stampo il risultato.
		cout << "Relative error with PALU decomposition:\t" << relative_error(x, xLu) << endl; //Calcolo e stampo  l'errore relativo.
	}
	return 0;
};

/*Risoluzione di  un sistema lineare(Ax = b) tramite la decomposizione LU con pivoting parziale + calcolo dell'errere relativo.
La funzione solve_system_PALU prende in input la matrice A, il vettore termine noto b e la soluzione esatta del sistema x.*/
int solve_system_QR(MatrixXd A, VectorXd b, VectorXd x) {
	Vector2d xQr = A.householderQr().solve(b); //Risolvo il sistema.
	cout << "Sol with QR decomposition:\n"<< xQr << endl; //Stampo la soluzione ottenuta.
	cout << "Relative error with QR decomposition:\t" << relative_error(x, xQr) << endl; //Calcolo e stampo l'errore relativo.
	return 0;
};


int main()
{
	Matrix2d A1; //1. A
	A1 << 5.547001962252291e-01, -3.770900990025203e-02,
		8.320502943378437e-01, -9.992887623566787e-01;
	Vector2d b1; //1. b
	b1 << -5.169911863249772e-01, 1.672384680188350e-01;
	Matrix2d A2; //2. A1
	A2 << 5.547001962252291e-01,-5.540607316466765e-01,
		8.320502943378437e-01,-8.324762492991313e-01;
	Vector2d b2; //2. b
	b2 << -6.394645785530173e-04, 4.259549612877223e-04;
	Matrix2d A3; //3. A
	A3 << 5.547001962252291e-01,-5.547001955851905e-01,
		8.320502943378437e-01,-8.320502947645361e-01;
	Vector2d b3; //3. b
	b3 <<  -6.400391328043042e-10, 4.266924591433963e-10;
	Vector2d x; //exact solution x
	x << -1.0e+0, -1.0e+00;
	
	/*1.*/
	cout << "1." << endl;
	solve_system_PALU(A1, b1, x);
	solve_system_QR(A1, b1, x);
	cout << "--------------------" << endl;
	/*2.*/
	cout << "2." << endl;
	solve_system_PALU(A2, b2, x);
	solve_system_QR(A2, b2, x);
	cout << "--------------------" << endl;
	/*3.*/
	cout << "3." << endl;
	solve_system_PALU(A3, b3, x);
	solve_system_QR(A3, b3, x);
	
    return 0;
}
