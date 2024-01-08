clc;
clear;
close all;

l = 1;
k = 6;

% Define the range of N values directly, avoiding the loop to generate it
N = 2.^25;

s = 1;

for n = N
    dx = l / n;
    x = linspace(0, l, n + 1);  % Use linspace for x vector
    
    phi_exact = sin(k * pi * x);
    
    % Calculate phi_bar directly using vectorized operations
    phi_bar = (-1 / (k * pi * 2 * dx)) * (cos(k * pi * (x + dx)) - cos(k * pi * (x - dx)));
    
    % Calculate phi_del_fd and phi_del_fv directly using vectorized operations
    phi_del_fd = (-1/8) * circshift(phi_exact, -1) + (3/4) * phi_exact + (3/8) * circshift(phi_exact, 2);
    phi_del_fv = (-1/6) * circshift(phi_bar, -1) + (5/6) * phi_bar + (1/3) * circshift(phi_bar, 2);
    
    % Calculate errors and orders of accuracy directly
    serr1_fd = abs(phi_del_fd - circshift(phi_exact, 1));
    serr1_fv = abs(phi_del_fv - circshift(phi_exact, 1));
    
    lerr1fd(s) = mean(serr1_fd);
    lerr1fv(s) = mean(serr1_fv);
    
    lerr2fd(s) = nthroot(mean(serr1_fd.^2), 2);
    lerr2fv(s) = nthroot(mean(serr1_fv.^2), 2);
    
    lerr3fd(s) = nthroot(mean(serr1_fd.^3), 3);
    lerr3fv(s) = nthroot(mean(serr1_fv.^3), 3);
    
    log_err1(s) = log10(lerr1fd(s));
    log_err2(s) = log10(lerr1fv(s));
    log_h(s) = log10(dx);
    
    fn(s) = n;
    s = s + 1;
end

% Plot errors vs. grid spacing
figure(1);
plot(log_h, log_err1, 'b', 'LineWidth', 2);
hold on;
plot(log_h, log_err2, 'r', 'LineWidth', 2);
title('NORM OF ERRORS VS GRID SPACING');
xlabel('log h');
ylabel('log E1');
legend('fd', 'fv');
hold off;

% Calculate and plot actual orders of accuracy
for i = 1:length(N)
    pfd(i) = log10((lerr3fd(i) - lerr2fd(i)) / (lerr2fd(i) - lerr1fd(i))) / log10(1/2);
    pfv(i) = log10((lerr3fv(i) - lerr2fv(i)) / (lerr2fv(i) - lerr1fv(i))) / log10(1/2);
end

figure(2);
plot(log_h, pfd, 'b', 'LineWidth', 2);
hold on;
plot(log_h, pfv, 'r', 'LineWidth', 2);
title('ACTUAL ORDER OF ACCURACY VS GRID SPACING');
xlabel('log h');
ylabel('p');
legend('fd', 'fv');
hold off;