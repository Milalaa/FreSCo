function imP = ImToPolar (imR, rMin, rMax, M, N)
% IMTOPOLAR converts rectangular image to polar form. The output image is 
% an MxN image with M points along the r axis and N points along the theta
% axis. The origin of the image is assumed to be at the center of the given
% image. The image is assumed to be grayscale.
% Bilinear interpolation is used to interpolate between points not exactly
% in the image.
%
% rMin and rMax should be between 0 and 1 and rMin < rMax. r = 0 is the
% center of the image and r = 1 is half the width or height of the image.
%
% V0.1 7 Dec 2007 (Created), Prakash Manandhar pmanandhar@umassd.edu

% function y=myfunction(a,b)其中a,b是输入函数的参数，y是函数返回的值。
% 当需要返回多个值时，可以将y看作一个数组，或者直接将函数的开头写成如function [x,y]=myfunction(x,y)的形式。


[Mr Nr] = size(imR); % size of rectangular image(input)
Om = (Mr+1)/2; % co-ordinates of the center of the image
On = (Nr+1)/2;
sx = (Mr-1)/2; % scale factors
sy = (Nr-1)/2;

imP  = zeros(M,  N);    % The output image is an MxN image with M points along the r axis and N points along the theta axis

delR = (rMax - rMin)/(M-1);
delT = 2*pi/N;

% loop in radius and 
for ri = 1:M
for ti = 1:N
    r = rMin + (ri - 1)*delR;
    t = (ti - 1)*delT;
    x = r*cos(t);   % In polar coordinates, the spacial-domain x
    y = r*sin(t);   % In polar coordinates, the spacial-domain y
    xR = x*sx + Om;  
    yR = y*sy + On; 
    imP (ri, ti) = interpolate (imR, xR, yR);   % Bilinear interpolation is used to interpolate between points not exactly in the image.
end
end

function v = interpolate (imR, xR, yR) % Bilinear interpolation is used to interpolate between points not exactly in the image.
    xf = floor(xR); % floor朝负无穷大方向取整；将x中元素取整，值y为不大于本身的最小整数。对于复数，分别对实部和虚部取整
    xc = ceil(xR);  % ceil表示向上取整的意思
    yf = floor(yR);
    yc = ceil(yR);
    if xf == xc & yc == yf
        v = imR (xc, yc);
    elseif xf == xc
        v = imR (xf, yf) + (yR - yf)*(imR (xf, yc) - imR (xf, yf));
    elseif yf == yc
        v = imR (xf, yf) + (xR - xf)*(imR (xc, yf) - imR (xf, yf));
    else
       A = [ xf yf xf*yf 1
             xf yc xf*yc 1
             xc yf xc*yf 1
             xc yc xc*yc 1 ];
       r = [ imR(xf, yf)
             imR(xf, yc)
             imR(xc, yf)
             imR(xc, yc) ];
       a = A\double(r);
       w = [xR yR xR*yR 1];
       v = w*a;
    end
