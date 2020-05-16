clc
clear
%% Точный алгоритм пересчёта координат и составляющих вектора скорости центра масс НКА на заданный момент времени шкалы МВД
%% Эфемериды НКА 21 системы ГЛОНАСС в системе ПЗ-90.11
%% Дата 2020/02/26 13:45:18
%Расчёт времени формата ГЛОНАСС
N4=floor((2020-1996)/4)+1;%Номер текущего четырёхлетия
Nt=365*(2020-1996-4*(N4-1))+31+25+1;%Номер текущих суток
tb=13*60*60 + 45*60 + 18 + 10800;%момент по шкале МДВ, к которому привязаны эфемериды ГЛОНАСС, в сек
Toe = (12+3)*60*60;
Tof = (24+3)*60*60;
Ts = 1;
ti = Toe:Ts:Tof;
X = -11998338.38; %% м
res = 2268666.50; %% м
Z = 22399278.81; %% м
VX = -2196.28239; %% м/с
VY = -2144.76013; %% м/с
VZ = -961.57646; %% м/с
AX = -0.0000056; %% м/с^2
AY = 0.0000009; %% м/с^2
AZ = -0.0000019; %% м/с^2
%% Константы
ae = 6378136; 
omega_z = 7.2921151467e-5;;
%% Расчёт текущей юлианской даты
JD0 = 1461*(N4 - 1) + Nt + 2450008.5 - (Nt -3)/25; 
%% Среднее звёздное время по Гринвичу
%% Время от эпохи 2000 года 1 января до текущей эпохи
T_delta=(JD0-2451545)/36525;
%% Угол поворота Земли, рад
ERA=2*pi*(0.7790572732640 + 1.00273781191135448*(JD0-2451545));
GMST=ERA+0.0000000703270726+0.0223603658710194*T_delta+...
    +0.0000067465784654*T_delta^2-0.0000000000021332*T_delta^3+...
    - 0.0000000001452308*T_delta^4-0.0000000000001784*T_delta^5;
%% Перевод в п/у инерциальную геоцентрическую систему 
S = GMST + omega_z*(tb - 3*60*60);
Xate=X*cos(S)-res*sin(S);
Yate=X*sin(S)+res*cos(S);
Zate=Z;

Vxate=VX*cos(S)-VY*sin(S)-omega_z*Yate;
Vyate=VX*sin(S)+VY*cos(S)+omega_z*Xate;
Vzate=VZ;

Axte=AX*cos(S)-AY*sin(S);
Ayte=AX*sin(S)+AY*cos(S);
Azte=AZ;
%% Метод Рунге-Кутты
%% Начальные условия 
res0 = [Xate Yate Zate Vxate Vyate Vzate];
[t, res] = ode45('diffs', tb:-Ts:ti(1), res0);
res1 = res(end:-1:2,:);
t1 = t(end:-1:2,:);
[t, res] = ode45('diffs', tb:Ts:ti(end), res0);
res1 = [res1; res];
t1 = [t1;t];
%% Учёт ускорений
tau1 = t1 - tb;
AXTE = AX*(tau1.^2)/2;
AYTE = AY*(tau1.^2)/2;
AZTE = AZ*(tau1.^2)/2;

delta_VX = AX*tau1;
delta_VY = AY*tau1;
delta_VZ = AZ*tau1;

delta_A = [AXTE AYTE AZTE delta_VX delta_VY delta_VZ];

res1 = res1 + delta_A;
%% Пересчёт координат в ПЗ-90
S = GMST + omega_z*(t1 - 3*60*60);
pz90(:,1) = res1(:,1).*cos(S) + res1(:,2).*sin(S);
pz90(:,2) = -res1(:,1).*sin(S) + res1(:,2).*cos(S);
pz90(:,3) = res1(:,3);
%% Координаты корпуса Е
%широта 55° 45' 24.0765",переводя получим 55.756687916667
N = 55.756687916667*pi/180;% широта [рад]
%долгота 37° 42' 11.0779" переводя в десятичные доли градуса получаем 37.703077194444
E=37.703077194444*pi/180;% долгота [рад]
H = 500; % высота [м]
cord_E = [N E H];
%Skyplot 
for i = 1:length(pz90(:,1))
    
    [X(i) Y(i) Z(i)] = ecef2enu(pz90(i,1),pz90(i,2),pz90(i,3),N,E,H,wgs84Ellipsoid,'radians');
    if Z(i) > 0
        r(i) = sqrt(X(i)^2 + Y(i)^2 + Z(i)^2);
        teta(i) = acos(X(i)/r(i));% 

        if X(i) > 0
            phi(i) = -atan(Y(i)/X(i))+pi/2;
        elseif (X(i)<0)&&(Y(i)>0)
            phi(i) = -atan(Y(i)/X(i))+3*pi/2;
        elseif (X(i)<0)&&(Y(i)<0)
            phi(i) = -atan(Y(i)/X(i))-pi/2;
        end
    else teta(i) = NaN;
        r(i) = NaN;
        phi(i) = NaN;
    end
end
%Построение графиков
figure(1)
[X,Y,Z]=sphere(50);
Rz=6371000;%радиус Земли
surf(Rz*X,Rz*Y,Rz*Z)
hold on
grid on
plot3(res1(:,1), res1(:,2), res1(:,3), 'b')
title('Траектория движения спутника ГЛОНАСС №21')
xlabel('Ось Х, м')
ylabel('Ось Y, м')
zlabel('Ось Z, м')
hold off
legend('Земля','ПЗ-90', 'Инерциальная СК');

figure(2)
surf(Rz*X,Rz*Y,Rz*Z)
hold on
grid on
plot3(pz90(:,1),pz90(:,2),pz90(:,3),'r')
title({'Траектория движения КА №21 ГЛОНАСС,' ; 'в системе координат ПЗ-90'})
xlabel('Ось Х, м')
ylabel('Ось Y, м')
zlabel('Ось Z, м')
hold off



%SkyPlot
figure (3)
pax = polaraxes;
polarplot(pax,phi,teta*180/pi,'r')
pax.ThetaDir = 'clockwise';
pax.ThetaZeroLocation = 'top';
title('SkyView спутника ГЛОНАСС №21')
th = hours(t1./3600-3);

figure(4);
grid on
hold on
plot(th,(-teta*180/pi+90),'DurationTickFormat','hh:mm:ss')
title('Угол места')
xlabel('Время в МДВ')
ylabel('Угол места спутника ГЛОНАСС №21, град')
