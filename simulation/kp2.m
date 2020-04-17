clc
clear
%% Точный алгоритм пересчёта координат и составляющих вектора скорости центра масс НКА на заданный момент времени шкалы МВД
%% Эфемериды НКА 21 системы ГЛОНАСС в системе ПЗ-90.11
%% Дата 2020/02/26 13:45:18
tb = 13*60*60 + 45*60 + 18; %% с
Toe = 12*60*60;
Tof = 24*60*60;
Ts = 1;
ti = Toe:Ts:Tof;
X = -11998338.38; %% м
Y = 2268666.50; %% м
Z = 22399278.81; %% м
VX = -2196.28239; %% м/с
VY = -2144.76013; %% м/с
VZ = -961.57646; %% м/с
AX = -0.0000056; %% м/с^2
AY = 0.0000009; %% м/с^2
AZ = -0.0000019; %% м/с^2
%% Константы
J02 = 1082625.75 * 10^-9;
ae = 6378136; 
omegaz = 7.2921151467;
GM = 398600441.8 * 10^6; %% const гравитационного поля Земли
Gm = 4902.799*10^9; %% const гравитационного поля Луны
Gs = 13271244*10^13; %% const гравитационного поля Солнца
am = 3.84385243*10^5;
em = 0.054900489;
im = 0.0898041080;
as = 1.49598*10^8;
es = 0.016719;
re = 6371e3;
N4 = 7; %% номер текущего четырёхлетия
Nt = 57; %% номер текущих суток
%% Начальные условия
x0 = X;
y0 = Y;
z0 = Z;
VX0 = VX;
VY0 = VY;
VZ0 = VZ;
r0 = sqrt(X^2 + Y^2 + Z^2);
GM1 = GM/r0^2;
x01 = x0/r0;
y01 = y0/r0;
z01 = z0/r0;
p = ae/r0;
res0 = [ x0 y0 z0 VX0 VY0 VZ0 ];
%% Расчёт текущей юлианской даты
JD0 = 1461*(N4 - 1) + Nt + 245008.5 - (Nt -3)/25;
%% Номер юлианского дня для текущей даты
JDN = JD0 + 0.5;
%% Среднее звёздное время по Гринвичу
T_delta = (JD0 - 2451545)/36525;
ERA = 2*pi*(0.7790572732640 + 1.00273781191135448*(JD0 - 2451545));
GMST = ERA + 0.0000000703270726 + 0.0223603658710194*T_delta + 0.0000067465784654*T_delta^2 - 0.0000000000021332*T_delta^3 - 0.0000000001452308*T_delta^4 - 0.0000000000001784*T_delta^5;
%% Перевод в п/у инерциальную геоцентрическую систему ускорений
S = GMST + omegaz*(tb -10800);
Jxoms = AX*cos(S) - AY*sin(S);
Jyoms = AX*sin(S) + AY*cos(S);
Jzoms = AZ;
%% Ускорения от лунных и солнечных гравитационных возмущений
%% Лунные ускорения
T = (JD0 + (tb - 10800)/86400 -2451545)/36525; %% время от эпохи 2000 года
qum = 2.3555557435 + 8328.6914257190*T + 0.0001545547*T^2; %% средняя аномалия Луны
omegam = 2.1824391966 - 33.7570459536*T + 0.0000362262*T^2; %% средняя долгота восходящего узла Луны
gstrih = 1.4547885346 + 71.0176852437*T - 0.0001801481*T^2; %% средняя долгота перигея Луны
eps = 0.4090926006 -0.0002270711*T; %% средний наклон эклиптики к экватору
%% Эксцентричная аномалия Луны
Em0 = qum;
Em = qum + em*sin(Em0);
while abs(Em - Em0) < 10^-8
    Em = qum + em*sin(Em);
end
sintetam = ((sqrt(1 - em^2))*sin(Em))/(1 - em*cos(Em));
costetam = (cos(Em) - em)/(1 - em*cos(Em));
ksi11 = sin(omegam)*cos(omegam*(1 - cos(im)));
ksi12 = 1 - sin(omegam*(1 - cos(im)))^2;
ksi = 1 - cos(omegam*(1 - cos(im)))^2;
dzeta = cos(omegam)*sin(im);
eta = sin(omegam)*sin(im);
eta11 = ksi*cos(eps) - dzeta*sin(eps);
eta12 = ksi11*cos(eps) + eta*sin(eps);
dzeta11 = ksi*sin(eps) + dzeta*cos(eps);
dzeta12 = ksi11*sin(eps) - eta*cos(eps);
ksim = (sintetam*cos(gstrih) + costetam*sin(gstrih))*ksi11 + (costetam*cos(gstrih) - sintetam*sin(gstrih))*ksi12;
etam = (sintetam*cos(gstrih) + costetam*sin(gstrih))*eta11 + (costetam*cos(gstrih) - sintetam*sin(gstrih)*eta12);
dzetam = (sintetam*cos(gstrih) + costetam*sin(gstrih))*dzeta11 + (costetam*cos(gstrih) - sintetam*sin(gstrih)*dzeta12);
rm = am*(1 - em*cos(Em));
x0m = x0/rm; 
y0m = y0/rm;
z0m = z0/rm;
Gm1 = Gm/rm^2;
delta_m = ((ksim - x0m)^2 + (etam - y0m)^2 + (dzetam - z0m)^2)^(3/2);
jxom = Gm1*((ksim - x0m)/delta_m - ksim);
jyom = Gm1*((etam - y0m)/delta_m - etam);
jzom = Gm1*((dzetam - z0m)/delta_m - dzetam);
%% Солнечные ускорения
qus = 6.2400601269 + 628.3019551714*T - 0.0000026820*T^2;
%% Эксцентричная аномалия Солнца
Es0 = qus;
Es = qus + es*sin(Es0);
while abs(Es - Es0) < 10^-8
    Es = qus + es*sin(Es);
end
sintetas = ((sqrt(1 - es^2))*sin(Es))/(1 - es*cos(Es));
costetas = (cos(Es) - es)/(1 - es*cos(Es));
omegas = -7.6281824375 + 0.0300101976*T + 0.0000079741*T^2;
ksis = costetas*cos(omegas) - sintetas*sin(omegas);
etas = (sintetas*cos(omegas) + costetas*sin(omegas))*cos(eps);
dzetas = (sintetas*cos(omegas) + costetas*sin(omegas))*sin(eps);
rs = as*(1 -es*cos(Es));
x0s = x0/rs;
y0s = y0/rs;
z0s = z0/rs;
Gs1 = Gs/rs^2;
delta_s = ((ksis - x0s)^2 + (etas - y0s)^2 + (dzetas - z0s)^2)^(3/2);
jxos = Gs1*((ksis - x0s)/delta_s - ksis);
jyos = Gs1*((etas - y0s)/delta_s - etas);
jzos = Gs1*((dzetas - z0s)/delta_s - dzetas);
%% Пересчёт координат в ПЗ-90
x = x0*cos(S) + y0*sin(S);
y = -x0*sin(S) + y0*cos(S);
z = z0;
VX = VX0*cos(S) + VY*sin(S) + omegaz*y;
VY = -VX0*sin(S) + VY*cos(S) - omegaz*x;
VZ = VZ0;
S = GMST + omegaz*(ti - 10800);
%% Переход к сферическим координатам
r = sqrt(x^2 + y^2 + z^2);
teta = acos(z/r);
phi = atan2(y, x);
polarplot(phi,teta*180/pi, 'r')
ax = polaraxes;
polarplot(ax,pi,teta*180/pi, 'r')
%% Система дифф уравнений
function dres = diffs(t, res)
dres = res(:);
dres(1) = res(1);
dres(2) = res(2);
dres(3) = res(3);

dres(4) = -GM1*x01 - 1.5*J02*GM1*x01*(p^2)*(1 - 5*z01^2) + jxos + jxom;
dres(5) = -GM1*y01 - 1.5*J02*GM1*y01*(p^2)*(1 - 5*z01^2) + jyos + jyom;
dres(6) = -GM1*z01 - 1.5*J02*GM1*z01*(p^2)*(3 - 5*z01^2) + jzos + jzom;
%% Метод Рунге-Кутты
interval = Toe:1:Tof;
[t, res] = ode45(interval, res0);
[t, res] = ode45('diffs', tb:-Ts:ti(1), res0);
res1 = res(end:-1:2,:);
[t, res] = ode45('diffs', tb:Ts:ti(end), res0);
res1 = [res1, res];
end