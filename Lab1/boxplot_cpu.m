function [outputArg1,outputArg2] = boxplot_cpu(N, T)
%boxplot_cpu(N, T) creates a figure with multiple plots for different runs 
%           of CPU benchmarking.
%       N - number of performance measurements in a single run with length
%           K.
%       T - K x M matrix, where rows correspond to different number of
%           performance measurements and M is the number of tests ran.
%   Detailed explanation goes here 

    figure()
    %boxplot(T'./[1, 5, 10, 50, 100], N) <- test to see how different N behave
    %at N = 1000
    for i = 1:numel(N)
        subplot(2,3,i)
        boxplot(T(:,i), N(i))
        ylabel("Microseconds", 'FontSize', 14)
    end

    sgtitle("Time for computations")
end

