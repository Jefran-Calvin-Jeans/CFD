clc
clear all
close all

l = 1;
k = 6;
N = 2.^(0:24); % Powers of 2

s = 1;
for n = N
    dx = l / n;
    x = (0:dx:l)';
    phi_exact = sin(k * pi * x);
    phi_bar = (-1 / (k * pi * 2 * dx)) * (cos(k * pi * (x + dx)) - cos(k * pi * (x - dx)));

    i = 2:2:(n + 1);
    phi_del_fd = (-1 / 6) * phi_exact(i - 2) + (5 / 6) * phi_exact(i) + (1 / 3) * phi_exact(i + 2);
    phi_del_fv = (-1 / 8) * phi_bar(i - 2) + (3 / 4) * phi_bar(i) + (3 / 8) * phi_bar(i + 2);

    serr1_fd = abs(phi_del_fd - phi_exact(i + 1));
    serr1_fv = abs(phi_del_fv - phi_exact(i + 1));
    serr11_fd = abs(phi_del_fd - phi_exact(i + 1) / phi_exact(i + 1));
    serr11_fv = abs(phi_del_fv - phi_exact(i + 1) / phi_exact(i + 1));

    lesum1fd = sum(serr1_fd);
    lesum1fv = sum(serr1_fv);
    lesum11fd = sum(serr11_fd);
    lesum11fv = sum(serr11_fv);
    lesum2fd = sum(serr11_fd.^2);
    lesum2fv = sum(serr11_fv.^2);
    lesum3fd = sum(serr11_fd.^3);
    lesum3fv = sum(serr11_fv.^3);

    lerr1fd(s) = (1 / (n / 2)) * lesum1fd;
    lerr1fv(s) = (1 / (n / 2)) * lesum1fv;
    lerr11fd(s) = (1 / (n / 2)) * lesum11fd;
    lerr11fv(s) = (1 / (n / 2)) * lesum11fv;
    lerr2fd(s) = nthroot((1 / (n / 2)) * lesum2fd, 2);
    lerr2fv(s) = nthroot((1 / (n / 2)) * lesum2fv, 2);
    lerr3fd(s) = nthroot((1 / (n / 2)) * lesum3fd, 3);
    lerr3fv(s) = nthroot((1 / (n / 2)) * lesum3fv, 3);

    log_err1(s) = log10(lerr1fd(s));
    log_err2(s) = log10(lerr1fv(s));
    log_h(s) = log10(dx);
    fn(s) = n;

    s = s + 1;
end

figure(1)
plot(log_h, log_err1, 'b', 'DisplayName', 'fd');
hold on
plot(log_h, log_err2, 'r', 'DisplayName', 'fv');
title('NORM OF ERRORS VS GRID SPACING');
xlabel('log h');
ylabel('log E1');
legend('Location', 'Best');
hold off

for i = 1:25
    pfd(i) = log10((lerr3fd(i) - lerr2fd(i)) / (lerr2fd(i) - lerr11fd(i))) / log10(1 / 2) * (2 / 3);
    pfv(i) = log10((lerr3fv(i) - lerr2fv(i)) / (lerr2fv(i) - lerr11fv(i))) / log10(1 / 2) * (2 / 3);
end

figure(2)
plot(log_h, pfd, 'b', 'DisplayName', 'fd');
hold on
plot(log_h, pfv, 'r', 'DisplayName', 'fv');
title('ACTUAL ORDER OF ACCURACY VS GRID SPACING');
xlabel('log h');
ylabel('p');
legend('Location', 'Best');
hold off