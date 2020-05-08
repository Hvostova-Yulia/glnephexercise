
%% Система дифф уравнений
function dres = diffs(t, res)
ae = 6378136; 
x0=res(1);
y0=res(2);
z0=res(3);
r=sqrt(x0^y0^2+z0^2);
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

dres(4) = -GM1*x01 - 1.5*J02*GM1*x01*(p^2)*(1 - 5*z01^2);
dres(5) = -GM1*y01 - 1.5*J02*GM1*y01*(p^2)*(1 - 5*z01^2);
dres(6) = -GM1*z01 - 1.5*J02*GM1*z01*(p^2)*(3 - 5*z01^2);
end
