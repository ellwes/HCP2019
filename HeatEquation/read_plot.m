%fileID = fopen('o.ans','w');
%fwrite(fileID,magic(3),'double');
%fclose(fileID);
clear all, close all, clc
filename = 'o.ans'
fileID = fopen(filename);
A = fread(fileID,'double')
T = A(1);
n = sqrt(length(A(2:end)));
A = reshape(A(2:end), [n, n]);
x = linspace(0, 1, n);
y = linspace(0, 1, n);
[X, Y] = meshgrid(x, y);
contourf(X, Y, A, 'EdgeColor', 'None')
axis equal;
xlabel('x');
ylabel('y');
title(['temperature at time: ' num2str(T)])
fclose(fileID);