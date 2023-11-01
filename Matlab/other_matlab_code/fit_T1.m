t = [0.5, 1, 2, 3, 4, 5, 6]';
y_data = [4.5496, 7.1583, 11.2189, 13.9505, 18.9563, 19.7784, 19.6695]';
fun = @(params, t) params(1) * (1 - exp(-t / params(2)));

% Initial guess for parameters (M and T1)
initial_guess = [1, 2]; % You can adjust these initial values

% Fit the data to the model using lsqcurvefit
fit_params = lsqcurvefit(fun, initial_guess, t, y_data);

% Extract the fitted parameters
fitted_M = fit_params(1);
fitted_T1 = fit_params(2);

% Plot the original data, true model, and fitted model
figure;
plot(t, y_data, 'LineWidth', 2);
hold on;
plot(t, fun(fit_params, t), 'r--', 'LineWidth', 2);
xlabel('TR (s)');
ylabel('Amplitude');
legend('Data', 'Fitted Model');
title('T1 fitted curve');
grid on;

% Display the fitted parameters
fprintf('Fitted M: %.4f\n', fitted_M);
fprintf('Fitted T1: %.4f\n', fitted_T1);