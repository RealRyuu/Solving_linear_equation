#define N 946
#define I 1
#include<cmath>
#include<chrono>
#include<fstream>
#include<cstdlib>

using namespace std::chrono;

void zero_vector(double* A) {
	for (int i = 0; i < N; i++) {
		A[i] = 0;
	}
}

void zero_matrix(double** A) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++)
			A[i][j] = 0;
	}
}

void vector_copy(double* TO, double* FROM) {
	for (int i = 0; i < N; i++)
		TO[i] = FROM[i]; //przepisywanie tu
}

void jacobi(double** A, double* b, double* X_prew, double* X_next) {
	for (int i = 0; i < N; i++) {
		X_next[i] = 0;
		for (int j = 0; j < N; j++) {
			if (j != i)
				X_next[i] -= A[i][j] * X_prew[j];
		}
		X_next[i] += b[i];
		X_next[i] /= A[i][i];
	}
}

void gauss_seidel(double** A, double* b, double* X_prew, double* X_next) {
	for (int i = 0; i < N; i++) {
		X_next[i] = 0;
		for (int j = 0; j < i; j++) {
			X_next[i] -= A[i][j] * X_next[j];
		}
		for (int j = i + 1; j < N; j++) {
			X_next[i] -= A[i][j] * X_prew[j];
		}
		X_next[i] += b[i];
		X_next[i] /= A[i][i];
	}
}

void gauss(double** A, double* b, double* X) {
	// reduce
	double m;
	for (int k = 0; k < N; k++) {
		for (int i = k + 1; i < N; i++) {
			m = A[i][k] / A[k][k];
			for (int j = 0; j < N; j++) {
				A[i][j] -= A[k][j] * m;
			}
			b[i] -= b[k] * m;
		}
	}
	// solve
	for (int i = N - 1; i >= 0; i--) {
		X[i] = b[i];
		for (int j = i + 1; j < N; j++) {
			X[i] -= A[i][j] * X[j];
		}
		X[i] /= A[i][i];
	}
}

void compute_residuum(double** A, double* b, double* X, double* res) {
	for (int i = 0; i < N; i++) {
		res[i] = 0;
		for (int j = 0; j < N; j++) {
			res[i] += A[i][j] * X[j];
		}
		res[i] -= b[i];
	}
}

double compute_norm2(double* res) {
	double norm = 0.0;
	for (int i = 0; i < N; i++)
		norm += res[i] * res[i];
	return sqrt(norm);
}

int main() {
	double a1, a2, a3;
	a1 = 5.0 + 6.0;
	a2 = a3 = -1.0;
	a1 = 3.0;
	
	/*double A[N][N] = {};
	double b[N][1] = {};
	double X_prew[N][1] = {};
	double X_next[N][1] = {};
	double res[N][1] = {};
	double norm;*/

	double** A = new double*[N];
	double* b = new double[N];
	double* X_prew = new double[N];
	double* X_next = new double[N];
	double* res = new double[N];
	double norm, norm2, norm3;

	high_resolution_clock::time_point start;
	high_resolution_clock::time_point end;
	double tmp;

	std::fstream file;
	file.open("rozbiezny2.csv", std::ios::out);

	for (int i = 0; i < N; i++) {
		A[i] = new double[N];
	}

	zero_matrix(A);
	zero_vector(res);

	for (int n = 0; n < N; n++) {
		A[n][n] = a1;
		if (n >= N - 1)
			continue;
		A[n][n + 1] = a2;
		A[n + 1][n] = a2;
		if (n >= N - 2)
			continue;
		A[n][n + 2] = a3;
		A[n + 2][n] = a3;
	}

	for (int n = 0; n < N; n++) {
		b[n] = sin((n + 1)*(0 + 1)/50);
	}

	int k = 0;
	double jacobi_time[I] = {0.0};
	double gauss_seidel_time[I] = { 0.0 };

	//goto skip;
	file << "jacobi time:" << ';' << ';' << "gauss-seidel time:" << '\n';
	for (int i = 0; i < I; i++) {
		goto skip;
		zero_vector(X_prew);
		zero_vector(X_next);
		do {
			k++;
			start = high_resolution_clock::now();
			jacobi(A, b, X_prew, X_next);
			end = high_resolution_clock::now();
			tmp = duration_cast<milliseconds>(end - start).count();
			jacobi_time[i] += tmp;

			//}

			//{obliczanie residium
			compute_residuum(A, b, X_next, res);

			norm = compute_norm2(res);
			file << norm << '\n';

			vector_copy(X_prew, X_next);
		} while (norm >= 10e-9);
		file << jacobi_time[i] << '\n';
	skip:
		zero_vector(X_prew);
		zero_vector(X_next);
		int g = 0;
		//start = high_resolution_clock::now();
		do {
			// gauss-siedl{
			g++;
			start = high_resolution_clock::now();
			gauss_seidel(A, b, X_prew, X_next); // <- gauss_seidel
			end = high_resolution_clock::now();
			tmp = duration_cast<milliseconds>(end - start).count();
			gauss_seidel_time[i] += tmp;;
			//}
			// przepisywanie{

			//}
			//{obliczanie residium
			compute_residuum(A, b, X_next, res);
			norm2 = compute_norm2(res);
			file << norm2 << '\n';

			vector_copy(X_prew, X_next);
		} while (norm2 >= 10e-9);
		file << gauss_seidel_time[i] << '\n';
	}

//skip:

	zero_vector(X_prew);
	zero_vector(X_next);

	double gauss_time;
	start = high_resolution_clock::now();
	gauss(A, b, X_next); //				<- gauss
	end = high_resolution_clock::now();
	gauss_time = duration_cast<milliseconds>(end - start).count();

	compute_residuum(A, b, X_next, res);
	norm3 = compute_norm2(res);
	file << '\n' << norm3 << '\n' << gauss_time;

	for (int i = 0; i < I; i++) {
		printf("%d\njacobi time: %lf\ngauss-seidl time: %lf\n\n",i, double(jacobi_time[i]), double(gauss_seidel_time[i]));
	}
	printf("gauss time: %lf\n", gauss_time);


	for (int i = 0; i < N; i++) {
		delete A[i];
	}

	delete b;
	delete X_prew;
	delete X_next;
	delete res;

	return 0;
}