#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <cmath>

double min(const double& a, const double& b) {
	if (a > b) {
		return b;
	}
	else {
		return a;
	}
}

double min(const double& a, const double& b, const double& c) {
	if (a < b && a < c) {
		return a;
	}
	else {
		if (b < a && b < c) {
			return b;
		}
		else {
			return c;
		}
	}
}

double max(const double& a, const double& b) {
	if (a > b) {
		return a;
	}
	else {
		return b;
	}
}

double sgn(const double& a) {
	if (a < 0) {
		return -1;
	}
	else {
		return 1;
	}
}

class TWD {

public:
	TWD() {
		T1 = 25;
		T2 = 4;
		L = 1;
		Nl = 100;
		h = L / Nl;
		n0 = 10000000000;
		k = 1.38 * 10000;
		m = 1.66;
		Ecut = 4.8 * sqrt(k * T1 / m);
		Ne = 100;
		t = 0.7 * h / Ecut;
		Nt = 100;
		T = Nt * t;
		E.resize(Ne);
		dE = 2 * Ecut / Ne;
		E[0] = -Ecut;
		for (double i = 1; i < Ne; ++i) {
			E[i] = E[i - 1] + dE;
		}
		f.resize(Ne);
		f_previous.resize(Ne);
		for (double j = 0; j < Ne; ++j) {
			f[j].resize(Nl);
			f_previous[j].resize(Nl);
		}
		for (int j = 0; j < Ne; ++j) {
			for (int i = 0; i < Nl; ++i) {
				f_previous[j][i] = n0 * 1 / (dE * speed_sum0(0, Ne)) * exp((-1) * m * E[j] * E[j] / (2 * k * T2));
			}
		}
		n.resize(Nl);
		Temp.resize(Nl);
		f12A.resize(Ne);
		f21A.resize(Ne);
		f12B.resize(Ne);
		f21B.resize(Ne);
		limiter_type = "MC limiter";
	}

	void print() {
		for (int j = 0; j < Ne; ++j) {
			for (int i = 0; i < Ne; ++i) {
				std::cout << f_previous[j][i] << " ";
			}
			std::cout << "\n";
		}
	}

	//основной алгоритм
	void Run() {
		for (double q = 0; q < Nt; ++q) {
			//1)вычисляем вспомогательные значения у левой границы для скоростей меньше 0
			for (double j = 0; j < Ne / 2; ++j) {
				double F0;
				F0 = max(0, 2 * f_previous[j][0] - f_previous[j][1]);
				double g = E[j] * t / h;
				double df = calculate_limiter(F0, f_previous[j][0], f_previous[j][1]);
				f12A[j] = f_previous[j][0] - (1 - g) / 2 * df;
			}
			//теперь на основе полученных значений можем вычислить вспомогательные значения 
			//2)у левой границы для скоростей больше 0
			for (double j = Ne / 2; j < Ne; ++j) {
				f12A[j] = EFsumA(0, Ne / 2) / speed_sum(Ne / 2, Ne) * exp(-m * E[j] * E[j] / (2 * k * T1));
			}
			//3)теперь делаем то же самое для значений f32, только теперь сначала вычисляем для положительных скоростей
			for (double j = Ne / 2; j < Ne; ++j) {
				double F0;
				F0 = max(0, 2 * f_previous[j][0] - f_previous[j][1]);
				double g = E[j] * t / h;
				double df = calculate_limiter(F0, f_previous[j][0], f_previous[j][1]);
				f21A[j] = f_previous[j][0] + (1 - g) / 2 * df;
			}
			for (double j = 0; j < Ne / 2; ++j) {
				f21A[j] = EFsum1(Ne / 2, Ne) / speed_sum(0, Ne / 2) * exp(-m * E[j] * E[j] / (2 * k * T1));
			}
			//теперь мы готовы вычислить граничные значения f у левой границы
			for (double j = 0; j < Ne; ++j) {
				f[j][0] = f_previous[j][0] - E[j] * t / h * (f21A[j] - f12A[j]);
			}
			//вычисляем значения внутри сосуда
			for (double j = 0; j < Ne; ++j) {
				double f12, f21;
				f12 = f12A[j];
				for (double i = 1; i < Nl - 2; ++i) {
					f21 = cfh(j, i);
					f[j][i] = f_previous[j][i] - E[j] * t / h * (f21 - f12);
					f12 = f21;
				}
			}
			//1) теперь вычислим вспомогательные значения у правой границы
			//сначала для положительных скоростей
			for (double j = Ne / 2; j < Ne; ++j) {
				double F0 = max(0, 2 * f_previous[j][Nl - 1] - f_previous[j][Nl - 2]);
				double g = E[j] * t / h;
				double df = calculate_limiter(f_previous[j][Nl - 2], f_previous[j][Nl - 1], F0);
				f21B[j] = f_previous[j][Nl - 1] + (1 - g) / 2 * df;
			}
			//2) теперь для отрицательных скоростей
			for (double j = 0; j < Ne / 2; ++j) {
				f21B[j] = EFsumB(Ne / 2, Ne) / speed_sum(0, Ne / 2) * exp(-m * E[j] * E[j] / (2 * k * T2));
			}
			//3) f32
			for (double j = 0; j < Ne / 2; ++j) {
				double F0 = max(0, 2 * f_previous[j][Nl - 1] - f_previous[j][Nl - 2]);
				double g = E[j] * t / h;
				double df = calculate_limiter(f_previous[j][Nl - 2], f_previous[j][Nl - 1], F0);
				f12B[j] = f_previous[j][Nl - 1] - (1 - g) / 2 * df;
			}
			for (double j = Ne / 2; j < Ne; ++j) {
				f12B[j] = EFsumB1(0, Ne / 2) / speed_sum(Ne / 2, Ne) * exp(-m * E[j] * E[j] / (2 * k * T2));
			}
			//значения при Nl - 2
			for (double j = 0; j < Ne; ++j) {
				double fe = cfh(j, Nl - 3);
				f[j][Nl - 2] = f_previous[j][Nl - 2] + E[j] * t / h * (f12B[j] - fe);
			}
			//и теперь значения у правой границы
			for (double j = 0; j < Ne; ++j) {
				f[j][Nl - 1] = f_previous[j][Nl - 1] - E[j] * t / h * (f21B[j] - f12B[j]);
			}
			f_previous = f;
			//каждый 10 шаг выводим значения
		/*	if (q == 10) {
				std::cout << "Time passed " << q * t << "\n" << "------------------------------------------------\n";
			/*	for (double j = 0; j < Ne; ++j) {
					std::cout << E[j] << " :  ";
					for (double i = 0; i < Nl; ++i) {
						std::cout << f[j][i] << " ";
					}
					std::cout << "\n";
				}  
			} */
		}
		//теперь, когда режим установился(прошло некоторое время), вычислим установившиеся значения концентрации
		std::cout << "Concentration\n" << "----------------------------------------------------------------------\n";
		for (double i = 0; i < Nl; ++i) {
			n[i] = n0 * f_sum0(i) * 1 / speed_sum0(0, Ne);
			std::cout << i * h << " :  " << n[i] << "\n";
		}
		//и температуры
		std::cout << "Temperature\n" << "----------------------------------------------------------------\n";
		for (double i = 0; i < Nl; ++i) {
			Temp[i] = T1 * f_sum(i) * 1 / f_sum0(i);
			std::cout << i * h << " :  " << Temp[i] << "\n";
		}
	}

private:

	double Ecut; //скорость обрезания
	double Ne;	//кол-во шагов разбиения скоростей
	double T; //время
	double L; //длина исследуемого отрезка
	double t; //шаг по времени
	double h; //шаг по координате
	double Nl; //кол-во шагов разбиения отрезка
	double Nt; //кол-во шагов разбиения времени
	std::string limiter_type;	//тип используемого в расчетах ограничителя
	double T1; //температура левой стенки
	double T2; //температура правой стенки
	double m; //масса молекулы
	double k; //постоянная больцмана
	double n0; //концентрация в начальный момент
	double dE; //шаг разбиения скоростей
	std::vector<std::vector<double>> f;	//координатная сетка в данный момент
	std::vector<std::vector<double>> f_previous; //координатная сетка в предыдущий момент
	std::vector<double> E;	//скоростная сетка
	std::vector<double> n; //концентрация
	std::vector<double> Temp; //температура
	std::vector<double> f12A; //промежуточные значения у границы А
	std::vector<double> f21A;
	std::vector<double> f12B;
	std::vector<double>f21B; //промежуточные значения у границы B

	void calculate_speed_field() {

	}

	double speed_sum0(const double& a, const double& b) {
		double result = 0;
		for (double i = a; i < b; ++i) {
			result += exp((-1) * m * E[i] * E[i] / (2 * k * T2));
		}
		return result;
	}

	double speed_sum(const double& a, const double& b) {
		double result = 0;
		for (double i = a; i < b; ++i) {
			result += abs(E[i]) * exp((-1) * m * E[i] * E[i] / (2 * k * T2));
		}
		return result;
	}

	//ограничители
	double minmod_limiter(const double& f0, const double& f1, const double& f2) {
		return min(abs(f2 - f1), abs(f1 - f0)) * sgn(f2 - f1);
	}

	double MC_limiter(const double& f0, const double& f1, const double& f2) {
		return min(abs(f2 - f0) / 2, abs(f2 - f1) * 2, abs(f1 - f0) * 2) * sgn(f2 - f1);
	}

	double Superbee_limiter(const double& f0, const double& f1, const double& f2) {
		return max(min(2 * abs(f2 - f1), abs(f1 - f0)), min(abs(f2 - f1), 2 * abs(f1 - f2)));
	}

	double van_Leer_limiter(const double& f0, const double& f1, const double& f2) {
		return 2 * (f2 - f1) * (f1 - f0) / (f2 - f0);
	}

	double calculate_limiter(const double& f0, const double& f1, const double& f2) {
		if (limiter_type == "minmod limiter") {
			return minmod_limiter(f0, f1, f2);
		}
		if (limiter_type == "MC limiter") {
			return MC_limiter(f0, f1, f2);
		}
		if (limiter_type == "Superbee_limiter") {
			return Superbee_limiter(f0, f1, f2);
		}
		if (limiter_type == "van Leer limiter") {
			return van_Leer_limiter(f0, f1, f2);
		}
	}

	//вычисляем промежуточные значения
	double cfh(const double& j, const double& i) {
		double g = E[j] * t / h;
		if (E[j] >= 0) {
			double df = calculate_limiter(f_previous[j][i - 1], f_previous[j][i], f_previous[j][i + 1]);
			return f_previous[j][i] + (1 - g) / 2 * df;
		}
		else {
			double df = calculate_limiter(f_previous[j][i], f_previous[j][i + 1], f_previous[j][i + 2]);
			return f_previous[j][i + 1] - (1 - g) / 2 * df;
		}
	}

	//начальные условия 
	void starter_conditions() {

	}
	//интеграл значений функции в точке по скоростям
	double f_sum0(const double& i) {
		double result = 0;
		for (double j = 0; j < Ne; ++j) {
			result += f[j][i];
		}
		return result;
	}
	double f_sum(const double& i) {
		double result = 0;
		for (double j = 0; j < Ne; ++j) {
			result += E[j] * E[j] * f[j][i];
		}
		return result;
	}
	double EFsumA(double a, double b) {
		double result = 0;
		for (int i = a; i < b; ++i) {
			result += abs(E[i]) * f12A[i];
		}
		return result;
	}
	double EFsum1(double a, double b) {
		double result = 0;
		for (int i = a; i < b; ++i) {
			result += abs(E[i]) * f21A[i];
		}
		return result;
	}

	double EFsumB1(double a, double b) {
		double result = 0;
		for (int i = a; i < b; ++i) {
			result += abs(E[i]) * f12B[i];
		}
		return result;
	}
	double EFsumB(double a, double b) {
		double result = 0;
		for (int i = a; i < b; ++i) {
			result += abs(E[i]) * f21B[i];
		}
		return result;
	}

};

int main() {
	TWD test;
	test.Run();
	return 0;
}