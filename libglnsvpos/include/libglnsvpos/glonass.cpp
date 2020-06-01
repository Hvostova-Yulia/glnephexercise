#include <iostream>
#include <cmath>
using namespace std;
int main()
{// Координаты
	X = -11998338.38;
	Y = 2268666.50;
	Z = 22399278.81;
	// Скорости
	VX = -2196.28239;
	VY = -2144.76013;
	VZ = -961.57646;
	// Ускорения
	AX = -0.0000056;
	AY = 0.0000009;
	AZ = -0.0000019;
	// Константы
	ae = 6378136;
	omega_z = 7.2921151467e-5;
	J02 = -1082.63 * 10 ^ -6;
	GM = 398600441.8e6;
	// Номер текущего четырёхлетия
	N4 = floor((2020 - 1996) / 4) + 1;
	// Номер текущих суток
	Nt = 365 * (2020 - 1996 - 4 * (N4 - 1)) + 31 + 25 + 1;
	// Расчёт текущей юлианской даты
	JD0 = 1461 * (N4 - 1) + Nt + 2450008.5 - (Nt - 3) / 25;
	// Момент по шкале МДВ, к которому привязаны эфемериды ГЛОНАСС
	tb = 13 * 60 * 60 + 45 * 60 + 18 + 10800;
	// Начальное время
	Toe = (12 + 3) * 60 * 60;
	// Конечное время
	Tof = (24 + 3) * 60 * 60;
	// Среднее звёздное время по Гринвичу
	// Время от эпохи 2000 года 1 января до текущей эпохи
	T_delta = (JD0 - 2451545) / 36525;
	// Угол поворота Земли, рад
	ERA = 2 * pi * (0.7790572732640 + 1.00273781191135448 * (JD0 - 2451545));
	GMST = ERA + 0.0000000703270726 + 0.0223603658710194 * T_delta + ...
		+ 0.0000067465784654 * T_delta ^ 2 - 0.0000000000021332 * T_delta ^ 3 + ...
		- 0.0000000001452308 * T_delta ^ 4 - 0.0000000000001784 * T_delta ^ 5;
	// Перевод в п/у инерциальную геоцентрическую систему
	S = GMST + omega_z * (tb - 3 * 60 * 60);
	Xate = X * cos(S) - Y * sin(S);
	Yate = X * sin(S) + Y * cos(S);
	Zate = Z;

	Vxate = VX * cos(S) - VY * sin(S) - omega_z * Yate;
	Vyate = VX * sin(S) + VY * cos(S) + omega_z * Xate;
	Vzate = VZ;

	Axte = AX * cos(S) - AY * sin(S);
	Ayte = AX * sin(S) + AY * cos(S);
	Azte = AZ;
	// Метод Рунге-Кутты
	// Начальные условия
	Y0 = [Xate Yate Zate Vxate Vyate Vzate];
	// Ситема диффуров
	int diffs(double tb)
	{
		double Xate1 = Y.X;// по х
		double Yate1 = Y.Y;// по у
		double Zate1 = Y.Z;// по z
		double r = sqrt(Xate1 ^ 2 + Yate1 ^ 2 + Zate1 ^ 2);
		double GM1 = GM / r ^ 2;
		double x01 = Xate / r;
		double y01 = Yate / r;
		double z01 = Zate / r;
		double p = ae / r;
		dY.X = Y.VX;
		dY.Y = Y.VY;
		dY.Z = Y.VZ;
		dY.VX = -GM1 * x01 - (double)1.5 * J02 * GM1 * x01 * p ^ 2 * (1 - (double)5 * z01 ^ 2);
		dY.VY = -GM1 * y01 - (double)1.5 * J02 * GM1 * y01 * p ^ 2 * (1 - (double)5 * z01 ^ 2);
		dY.VZ = -GM1 * z01 - (double)1.5 * J02 * GM1 * z01 * p ^ 2 * (3 - (double)5 * z01 ^ 2);
		return dY;
	}
	int Runge_Kutta(N, double h)
	{
		if (N == 0 || h == 0) return 0;
		for (K = 1; K < N; K++)
			K1 = diffs(0, Y[K - 1])
			Y2.X = Y[K - 1].X + h * K1.X / (double)2;
		Y2.Y = Y[K - 1].Y + h * K1.Y / (double)2;
		Y2.Z = Y[K - 1].Z + h * K1.Z / (double)2;
		Y2.VX = Y[K - 1].VX + h * K1.VX / (double)2;
		Y2.VY = Y[K - 1].VY + h * K1.VY / (double)2;
		Y2.VZ = Y[K - 1].VZ + h * K1.VZ / (double)2;

		K2 = diffs(0 + h / 2, Y2);

		Y3.X = Y[K - 1].X + h * K2.X / (double)2;
		Y3.Y = Y[K - 1].Y + h * K2.Y / (double)2;
		Y3.Z = Y[K - 1].Z + h * K2.Z / (double)2;
		Y3.VX = Y[K - 1].VX + h * K2.VX / (double)2;
		Y3.VY = Y[K - 1].VY + h * K2.VY / (double)2;
		Y3.VZ = Y[K - 1].VZ + h * K2.VZ / (double)2;

		K3 = diffs(0 + h / 2, Y3);

		Y4.X = Y[K - 1].X + h * K3.X;
		Y4.Y = Y[K - 1].Y + h * K3.Y;
		Y4.Z = Y[K - 1].Z + h * K3.Z;
		Y4.VX = Y[K - 1].VX + h * K3.VX;
		Y4.VY = Y[K - 1].VY + h * K3.VY;
		Y4.VZ = Y[K - 1].VZ + h * K3.VZ;

		K4 = diffs(0 + h, Y4);

		Y[K].X = Y[K - 1].X + h / (double)6 * (K1.X + 2 * K2.X + 2 * K3.X + K4.X);
		Y[K].Y = Y[K - 1].X + h / (double)6 * (K1.Y + 2 * K2.Y + 2 * K3.Y + K4.Y);
		Y[K].Z = Y[K - 1].X + h / (double)6 * (K1.Z + 2 * K2.Z + 2 * K3.Z + K4.Z);
		Y[K].VX = Y[K - 1].X + h / (double)6 * (K1.VX + 2 * K2.VX + 2 * K3.VX + K4.VX);
		Y[K].VY = Y[K - 1].X + h / (double)6 * (K1.VY + 2 * K2.VY + 2 * K3.VY + K4.VY);
		Y[K].VZ = Y[K - 1].X + h / (double)6 * (K1.VZ + 2 * K2.VZ + 2 * K3.VZ + K4.VZ);
		cout << endl;
		for (int K = 0; K <= N; K++) {
			cout << "Y[" << K << "]=" << Y[K] << " ";
		}
		return 0;
}// Выходной массив
	{
		int Y[K];
		int Y_out[K];
		cin >> K;

		for (int i = 0; i < N; i++)
		{
			cin >> Y[i];
		}

		for (int i = 0; i < N; i++)
		{
			if (Y[k] != K)
			{
				Y_out[i] = Y[i];
			}
		}

		for (int i = 0; i < N; i++)
		{
			cout << Y_out[i];
		}
	}

	}

// Учёт ускорений
		}
		double tau1 = (double)Toe - (double)tb;
		for (i = 0; i < N; i++) {

			AXTE = AX * (tau1 ^ 2) /(double)2;
			AYTE = AY * (tau1 ^ 2) /(double)2;
			AZTE = AZ * (tau1 ^ 2) /(double)2;

			delta_VX = AX * tau1;
			delta_VY = AY * tau1;
			delta_VZ = AZ * tau1;

			delta_Y = [AXTE AYTE AZTE delta_VX delta_VY delta_VZ];

			Y_out1[i] = Y_out[i] + delta_Y;
		}
}