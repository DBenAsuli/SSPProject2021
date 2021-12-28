%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SSP Project 2021
% MATLAB Exercise
% Part 1
% Dvir Ben Asuli, Assaf Gadish
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_a_mone = [1 -0.4 -0.69 -0.964 -0.714];
x_b_mehane = [1 0.6 0.34 -0.034 -0.03807 -0.4212];

%%% Question 1
[Sxx_a, Sxx_b] = question1(x_a_mone, x_b_mehane);

%%% Question 2
% question2(x_a_mone, x_b_mehane);

%%% Question 3
question3(x_a_mone, x_b_mehane, Sxx_a, Sxx_b);

function [] = question3(x_a_mone, x_b_mehane, Sxx_a, Sxx_b)
    % Section a: biased correlogram estimator
    [xa, xb] = generate_xa_xb(x_a_mone, x_b_mehane);
    question3_parts_abc(xa, Sxx_a, 'x_a');
    question3_parts_abc(xb, Sxx_b, 'x_b');
end

function [] = question3_parts_abc(x, Sxx, varname)
    N = 400;
    K = 1024;
    Rxx = estimate_Rxx(x, N);
    Sxx_pos_correlogram = get_positive_fft(Rxx, K);
    figure;
    plot_positive_spectrum(Sxx_pos_correlogram, 'S_{xx}^B')

    % Section b: periodogram estimator
    X_k = get_positive_fft(x, K);
    Sxx_pos_periodogram = 1 / N * (X_k .* conj(X_k));
    figure;
    plot_positive_spectrum(Sxx_pos_periodogram, 'S_{xx}^P')

    % Section c: a,b with original
    figure;
    % 'S_{xx}^B'
    title(strcat('S_{xx} for', {' '}, varname, '[n]'))
    hold on
    plot_positive_spectrum_graph_only(Sxx_pos_correlogram);
    % 'S_{xx}^P'
    plot_positive_spectrum_graph_only(Sxx_pos_periodogram);
    plot_positive_spectrum_graph_only(Sxx);
    hold off
    legend('biased correlogram', 'periodogram', 'true spectrum')
end

function [] = plot_positive_spectrum(Sxx, s_name)
    plot_positive_spectrum_graph_only(Sxx);
    title(strcat('|', s_name, '(e^{j\omega})|(\omega)'))
    xlabel('\omega [rad]');
    ylabel(strcat('|', s_name, '(e^{j\omega})| [dB]'));
end

function [] = plot_positive_spectrum_graph_only(Sxx)
    w = linspace(0, pi, length(Sxx));
    plot(w, abs(Sxx));
end


function ffted_x_positive = get_positive_fft(x, K)
    x_padded = [x' zeros(1, 2 * K - size(x, 1))];
    ffted_x = fftshift(fft(x_padded));
    ffted_x_positive = ffted_x(K + 1 : 2 * K);
end

function Rxx = estimate_Rxx(x, N)
    Rxx = zeros((2 * N) - 1, 1);
    for l = (0: N - 1)
        curr = estimate_Rxx_l(x, l, N);
        Rxx(N - l) = curr;
        Rxx(N + l) = curr;
    end
end

function Rxx_l = estimate_Rxx_l(x, l, N)
    abs_l = abs(l);
    Rxx_l = 0;
    for n = (1 : N - abs_l)
        Rxx_l = Rxx_l + x(n) * x(n + abs_l);
    end
    Rxx_l = Rxx_l / N;
end    

function [s_xx_a, s_xx_b] = question1(x_a_mone, x_b_mehane)
    % Question 1 section a
    x_a_mone_roots = roots(x_a_mone);
    s_xx_a_zeros = add_conj_inverse_roots(x_a_mone_roots);
    fprintf('S_xx_a zeros:\n');
    figure;
    disp(s_xx_a_zeros)
    zplane(s_xx_a_zeros, []);
    
    x_b_mehane_roots = roots(x_b_mehane);
    s_xx_b_poles = add_conj_inverse_roots(x_b_mehane_roots);
    fprintf('s_xx_b poles:\n');
    figure;
    disp(s_xx_b_poles)
    zplane([], s_xx_b_poles);

    % Question 1 section b
    N = 1024;
    %w = linspace(0, pi, N);
    s_xx_a_mone = poly(s_xx_a_zeros);
    s_xx_a_mehane = [1];
    s_xx_a = freqz(s_xx_a_mone, s_xx_a_mehane, N);
    w = linspace(0, pi, length(s_xx_a));
    
    figure;
    plot(w, abs(s_xx_a));
    title('Spectrum of the MA(4) process, x_{a}[n]')
    xlabel('\omega [rad]');
    ylabel('|S_{xx_a}(e^{j\omega})|');
    
    s_xx_b_mone = [1];
    s_xx_b_mehane = poly(s_xx_b_poles);
    s_xx_b = freqz(s_xx_b_mone, s_xx_b_mehane, N);
    w = linspace(0, pi, length(s_xx_b));
    
    figure;
    plot(w, abs(s_xx_b));
    title('Spectrum of the AR(5) process, x_{b}[n]')
    xlabel('\omega [rad]');
    ylabel('|S_{xx_b}(e^{j\omega})|');
end

function [] = question2(x_a_mone, x_b_mehane)
    [xa, xb] = generate_xa_xb(x_a_mone, x_b_mehane);
        
    figure;
    plot(xa);
    title('x_a[n]');
    xlabel('n');
    
    figure;
    plot(xb);
    title('x_b[n]');
    xlabel('n');
end

function [xa, xb] = generate_xa_xb(x_a_mone, x_b_mehane)
    wa = randn(404, 1);
    wb = randn(2000, 1);

    xa_padded = filter(x_a_mone, 1, wa);
    xa = xa_padded(5:404);
    xb_padded = filter(1, x_b_mehane, wb);
    xb = xb_padded(1601:2000);
end



function [all_roots] = add_conj_inverse_roots(roots)
    conj_inverse_roots = 1 ./ conj(roots);
    all_roots = transpose([roots' conj_inverse_roots']);
end




