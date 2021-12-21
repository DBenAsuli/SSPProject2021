%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SSP Project 2021
% MATLAB Exercise
% Part 1
% Dvir Ben Asuli, Assaf Gadish
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Question 1
x_a_mone = [1 -0.4 -0.69 -0.964 -0.714];
x_b_mehane = [1 0.6 0.34 -0.034 -0.03807 -0.4212];
% Section a
x_a_mone_roots = roots(x_a_mone);
s_xx_a_zeros = add_conj_inverse_roots(x_a_mone_roots);
fprintf('S_xx_a zeros:\n');
disp(s_xx_a_zeros)
zplane(s_xx_a_zeros, []);

figure
x_b_mehane_roots = roots(x_b_mehane);
s_xx_b_poles = add_conj_inverse_roots(x_b_mehane_roots);

fprintf('s_xx_b poles:\n');
disp(s_xx_b_poles)
zplane([], s_xx_b_poles);

% Section b
N = 1024;
%w = linspace(0, pi, N);
s_xx_a_mone = poly(s_xx_a_zeros);
s_xx_a_mehane = [1];
s_xx_a = freqz(s_xx_a_mone, s_xx_a_mehane);
w = linspace(0, pi, length(s_xx_a));

figure;
plot(w, abs(s_xx_a));
title('Spectrum of the MA(4) process, x_{a}[n]')
xlabel('\omega [rad]');
ylabel('|S_{xx_a}(e^{j\omega})|');

s_xx_b_mone = [1];
s_xx_b_mehane = poly(s_xx_b_poles);
s_xx_b = freqz(s_xx_b_mone, s_xx_b_mehane);
w = linspace(0, pi, length(s_xx_b));

figure;
plot(w, abs(s_xx_b));
title('Spectrum of the AR(5) process, x_{b}[n]')
xlabel('\omega [rad]');
ylabel('|S_{xx_b}(e^{j\omega})|');



function [all_roots] = add_conj_inverse_roots(roots)
    conj_inverse_roots = 1 ./ conj(roots);
    all_roots = transpose([roots' conj_inverse_roots']);
end