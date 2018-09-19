function Ex = jhang( b, y, z, s, flow, u )
%% "Jhang 2015" Longitudinal dispersion coefficient
% This function Zeng and Huai 2015  LDC
% use meter and qubic meter as your dimentions
a   = y .* ( b + z .* y);
T   =b + 2 .* z .* y;
p   =b + 2 .* y .* sqrt(1+z .^ 2);
r   =a  ./  p;
D   =a  ./  T; %#ok
switch nargin   
    case 5
    u   = flow ./ a;
end
g   =9.81;
w   =a  ./  y;
ustar   =sqrt(g .* r .* s);
Ex      = y .* ustar .* 3.563 .* ((u ./  (sqrt(g .* y))) .^ -0.4117) .* ((w  ./  y) .^ 0.6776) .* ((u  ./  ustar) .^ 1.0132);
end