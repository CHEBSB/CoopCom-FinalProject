R = 500;        % cell radius
% reference hexagon
X0 = linspace(-500, 500, 1001);
Y0 = zeros(1, 1001);
for i = 1:250
    Y0(i) = sqrt(3) * (R + X0(i));
end
for i = 251:750
    Y0(i) = sqrt(3) * R / 2;
end
for i = 751:1001
    Y0(i) = sqrt(3) * (R - X0(i));
end
Ym = -1 * Y0;
XL = linspace(-250, 250, 501);
% border of sectors
L1 = zeros(1, 1001);
L2 = zeros(1, 501);        
L3 = zeros(1, 501);
for i = 1:501
    L2(i) = sqrt(3) * XL(i);
    L3(i) = -sqrt(3) * XL(i);
end
% plot
plot(X0, Y0);
hold on
plot(X0, Ym);
plot(X0, L1);
plot(XL, L2);
plot(XL, L3);

