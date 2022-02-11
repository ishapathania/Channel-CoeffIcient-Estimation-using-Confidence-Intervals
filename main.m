% ECE 514 Mini-project Fall 2021 
clear;
close all;
clc;

seed = 576; %Name: Isha Pathania, Vowels: a=16,i=256 (256+16+16+16+256+16=576)  
rng(seed,'twister');

C=10;  % Given Channel Coefficient
m =1000;  % No. of trials
n1 = 10; % n=10 recieved samples
n2 = 100; % n=100 recieved samples
sigma_N = 20;  % Given standard deviation of noise N
sigma_a = sigma_N; % Case a: Variance of N is the same as X 
conf = 0.9; % Given 90% confidence interval

alpha = 1-conf;
case_a_n1 = 0;  % Number of cases where C lies in the confidence interval for case a, n=10
case_b_n1 = 0;  % Number of cases where C lies in the confidence interval for case b, n=10
case_c_n1 = 0;  % Number of cases where C lies in the confidence interval for case c, n=10
sample_mean_array = zeros(1,10);   % Sample means from first 10 trials
lower_conf_a = zeros(1,10);  % Confidence interval lower limit for first 10 trials for case a, n=10
upper_conf_a = zeros(1,10);  % Confidence interval upper limit for first 10 trials for case a, n=10
lower_conf_b = zeros(1,10);  % Confidence interval lower limit for first 10 trials for case b, n=10
upper_conf_b = zeros(1,10);  % Confidence interval upper limit for first 10 trials for case b, n=10
lower_conf_c = zeros(1,10);  % Confidence interval lower limit for first 10 trials for case c, n=10
upper_conf_c = zeros(1,10);  % Confidence interval upper limit for first 10 trials for case c, n=10
case_a_n2 = 0;  % Number of cases where C lies in the confidence interval for case a, n=100
case_b_n2 = 0;  % Number of cases where C lies in the confidence interval for case b, n=100
case_c_n2 = 0;  % Number of cases where C lies in the confidence interval for case c, n=100

%running for loop for n=10 trials
for i = 1:m
    N = sigma_N*randn(1,n1); % Generating noise N and scaling by sigma_N
    X = C+N;  
    sample_mean = sum(X)/n1; % Sample Mean for the generated data points when n=10

    %Case a : Estimating the confidence interval using the variance of X
    delta_a = norminv((1-alpha/2),0,sigma_a/sqrt(n1));
    y_alpha_a = delta_a * sqrt(n1)/sigma_a;  
    if ((sample_mean-delta_a <= C) && (C <= sample_mean+delta_a))
        case_a_n1 = case_a_n1+1;  %count of number of times C lies within the confidence interval
    end
    
    %Case b : Estimating the confidence interval using sample variance
    %without assuming measurements are gaussian
    sample_var = (sum((X-sample_mean).^2))/(n1-1);
    sigma_b = sqrt(sample_var);
    y_alpha_b = y_alpha_a;
    delta_b = y_alpha_b*sigma_b/sqrt(n1);  
    if ((sample_mean-delta_b <= C) && (C <= sample_mean+delta_b))
        case_b_n1 = case_b_n1+1;  %count of number of times C lies within the confidence interval
    end
   
    %Case c : Estimating the confidence interval using Sample variance while assuming that the measurements are gaussian
    y_alpha_c = tinv(1-alpha/2,n1-1);
    sigma_c = sigma_b;
    delta_c = y_alpha_c*sigma_c/sqrt(n1);  
    if ((sample_mean-delta_c <= C) && (C <= sample_mean+delta_c))
        case_c_n1 = case_c_n1+1;  %count of number of times C lies within the confidence interval
    end
    
    if (i<=10)
        sample_mean_array(i) = sample_mean;
        lower_conf_a(i) = sample_mean - delta_a;
        upper_conf_a(i) = sample_mean + delta_a;
        lower_conf_b(i) = sample_mean - delta_b;
        upper_conf_b(i) = sample_mean + delta_b;
        lower_conf_c(i) = sample_mean - delta_c;
        upper_conf_c(i) = sample_mean + delta_c;
    end

end

disp("For n=10 samples");
disp("The number of times C lies within the 90% confidence interval for Case a is " + case_a_n1);
disp("The number of times C lies within the 90% confidence interval for Case b is " + case_b_n1);
disp("The number of times C lies within the 90% confidence interval for Case c is " + case_c_n1);

figure;
x_data = linspace(1,10,10);
x_line = linspace(0,12,12);
y_line = ones(1, 12)*C;
plot(x_line, y_line, 'g'); hold on;
error_a = upper_conf_a - sample_mean_array;
errorbar(x_data, sample_mean_array, error_a, 'k.'); % Plots error bars at the specific data points given by x_data and y_data.
plot(x_data, sample_mean_array, 'b.', 'MarkerSize', 15); hold on;% Add data points to the plot.

xlim([0 11]);
ylim([-40 40]);
title('Case a (n=10 recieved samples)');
xlabel('Number of Trials');
legend('Channel Coefficient', 'Confidence Interval', 'Sample Mean');


figure;
x_data = linspace(1,10,10);
x_line = linspace(0,12,12);
y_line = ones(1, 12)*C;
plot(x_line, y_line, 'g'); hold on;
error_b = upper_conf_b - sample_mean_array;
errorbar(x_data, sample_mean_array, error_b, 'k.'); % Plots error bars at the specific data points given by x_data and y_data.
plot(x_data, sample_mean_array, 'b.', 'MarkerSize', 15); hold on;% Add data points to the plot.

xlim([0 11]);
ylim([-40 40]);
title('Case b (n=10 recieved samples)');
xlabel('Number of Trials');
legend('Channel Coefficient', 'Confidence Interval', 'Sample Mean');

figure;
x_data = linspace(1,10,10);
x_line = linspace(0,12,12);
y_line = ones(1, 12)*C;
plot(x_line, y_line, 'g'); hold on;
error_c = upper_conf_c - sample_mean_array;
errorbar(x_data, sample_mean_array, error_c, 'k.'); % Plots error bars at the specific data points given by x_data and y_data.
plot(x_data, sample_mean_array, 'b.', 'MarkerSize', 15); hold on;% Add data points to the plot.

xlim([0 11]);
ylim([-40 40]);
title('Case c (n=10 recieved samples)');
xlabel('Number of Trials');
legend('Channel Coefficient', 'Confidence Interval', 'Sample Mean');


%running for loop for n=100 trials
for i = 1:m
    N = sigma_N*randn(1,n2); % Generating noise N and scaling by sigma_N
    X = C+N;  
    sample_mean = sum(X)/n2; % Sample Mean for the generated data points when n=100

    %Case a : Estimating the confidence interval using the variance of X
    delta_a = norminv((1-alpha/2),0,sigma_a/sqrt(n2));
    y_alpha_a = delta_a * sqrt(n2)/sigma_a;  
    if ((sample_mean-delta_a <= C) && (C <= sample_mean+delta_a))
        case_a_n2 = case_a_n2+1;  %count of number of times C lies within the confidence interval
    end
    
    %Case b : Estimating the confidence interval using sample variance
    %without assuming measurements are gaussian
    sample_var = (sum((X-sample_mean).^2))/(n2-1);
    sigma_b = sqrt(sample_var);
    y_alpha_b = y_alpha_a;
    delta_b = y_alpha_b*sigma_b/sqrt(n2);  
    if ((sample_mean-delta_b <= C) && (C <= sample_mean+delta_b))
        case_b_n2 = case_b_n2+1;  %count of number of times C lies within the confidence interval
    end
   
    %Case c : Estimating the confidence interval using Sample variance while assuming that the measurements are gaussian
    y_alpha_c = tinv(1-alpha/2,n2-1);
    sigma_c = sigma_b;
    delta_c = y_alpha_c*sigma_c/sqrt(n2);  
    if ((sample_mean-delta_c <= C) && (C <= sample_mean+delta_c))
        case_c_n2 = case_c_n2+1;  %count of number of times C lies within the confidence interval
    end
end

fprintf("\n")
disp("For n=100 samples");
disp("The number of times C lies within the 90% confidence interval for Case a is " + case_a_n2);
disp("The number of times C lies within the 90% confidence interval for Case b is " + case_b_n2);
disp("The number of times C lies within the 90% confidence interval for Case c is " + case_c_n2);

alpha=0.91;
delta=norminv((1-alpha/2),0,sqrt(295.84/121));
disp(delta);