clear all, close all

N = [1000, 5000, 10000, 50000, 100000];

% time in (s)
T = [0.00000224 , 0.00000246 ,  0.00000230 , 0.00000232
    0.00001199 , 0.00001214 ,  0.00001240 , 0.00001186
    0.00002359 , 0.00002426 ,  0.00002478 , 0.00002657
    0.00012287 , 0.00011941 ,  0.00012020 , 0.00011826
    0.00023717 , 0.00023820 ,  0.00023655 , 0.00023649];

% time transformed into microseconds
T = T*1e6;

err = std(T'); %error as standard deviation
m = mean(T');  %mean of each value
errorbar(N, m, err)
hold on
plot(N, m, '*r')
ylabel("microseconds", 'FontSize', 14)
xlabel("Number of multiplications", 'FontSize', 14)
title("CPU time with cold-start average over multiple runs.")


T = table(T, m', err');
T.Properties.RowNames = string(N');
T.Properties.VariableNames = {'CPUTime','Mean','Deviation'}