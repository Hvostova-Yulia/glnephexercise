%% онстанты
J02 = 1082625.75 * 10^-9;
GM = 398600441.8 * 10^6; %% const гравитационного пол€ «емли
Gm = 4902.799*10^9; %% const гравитационного пол€ Ћуны
Gs = 13271244*10^13; %% const гравитационного пол€ —олнца
am = 3.84385243*10^5;
em = 0.054900489;
im = 0.0898041080;
as = 1.49598*10^8;
es = 0.016719;
re = 6371e3;
%% ”скорени€ от лунных и солнечных гравитационных возмущений
%% Ћунные ускорени€
T = (JD0 + (tb - 10800)/86400 -2451545)/36525; %% врем€ от эпохи 2000 года
qum = 2.3555557435 + 8328.6914257190*T + 0.0001545547*T^2; %% средн€€ аномали€ Ћуны
omegam = 2.1824391966 - 33.7570459536*T + 0.0000362262*T^2; %% средн€€ долгота восход€щего узла Ћуны
gstrih = 1.4547885346 + 71.0176852437*T - 0.0001801481*T^2; %% средн€€ долгота периге€ Ћуны
eps = 0.4090926006 -0.0002270711*T; %% средний наклон эклиптики к экватору
%% Ёксцентрична€ аномали€ Ћуны
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
%% —олнечные ускорени€
qus = 6.2400601269 + 628.3019551714*T - 0.0000026820*T^2;
%% Ёксцентрична€ аномали€ —олнца
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
%% —истема дифф уравнений
function dres = diffs(t, res)
ae = 6378136; 
x0=res(1);
y0=res(2);
z0=res(3);
r=sqrt(Xate^2+Yate^2+Zate^2);
GM1 = GM/r0^2;
x01 = x0/r0;
y01 = y0/r0;
z01 = z0/r0;
r0=ae/r;
p = ae/r0;

dres = res(:);
dres(1) = res(1);
dres(2) = res(2);
dres(3) = res(3);

dres(4) = -GM1*x01 - 1.5*J02*GM1*x01*(p^2)*(1 - 5*z01^2) + jxos + jxom;
dres(5) = -GM1*y01 - 1.5*J02*GM1*y01*(p^2)*(1 - 5*z01^2) + jyos + jyom;
dres(6) = -GM1*z01 - 1.5*J02*GM1*z01*(p^2)*(3 - 5*z01^2) + jzos + jzom;
end
