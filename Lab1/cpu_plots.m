clear all, close all

load('cputime.mat') % data we got
% time transformed into microseconds
T = T*1e6;

err = std(T); %error as standard deviation
m = mean(T);  %mean of each value

% errorplot
figure()
%errorbar(N, m, err)
plot(N, m, N, m, '*')
hold on
plot(N, m, '*r')
ylabel("microseconds", 'FontSize', 14)
xlabel("Number of multiplications", 'FontSize', 14)
title("CPU time")

% boxplot
boxplot_cpu(N, T)

% table
myTable = table(m', err');
myTable.Properties.RowNames = string(N');
myTable.Properties.VariableNames = {'Mean','Deviation'}