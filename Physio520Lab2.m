%% Before: Basic Random Variable Analysis
% Generate 1000 random variables from a standard normal distribution
x = randn(1,1000);

% Calculate and display the mean and standard deviation
disp(['Mean of x: ', num2str(mean(x))]);
disp(['Standard deviation of x: ', num2str(std(x))]);

% Plot histogram of the random variables
figure;
histogram(x, 10);
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Frequency', 'Interpreter', 'latex', 'FontSize', 18);
title('Histogram of Random Variables');







%% Exercise 1: Statistical Properties of Random Variables
N = 1000; % Length of the series
x = randn(1, N);

% Calculate and display expected value of x_i
disp(['Expected value of x_i: ', num2str(mean(x))]);

% Calculate and display expected value of x_i^2
disp(['Expected value of x_i^2: ', num2str(mean(x.^2))]);

% Calculate and display expected value of x_i * x_{i+1}
disp(['Expected value of x_i * x_{i+1}: ', num2str(mean(x(1:N-1) .* x(2:N)))]);

%% Exercise 2: Plotting Random Walks
% Plot a single random walk
N = 99; % Length of random walk
x = randn(1, N);
X = zeros(1, N+1); % Initialize position
X(1) = 0; % Initial position

for i = 2:N+1
    X(i) = X(i-1) + x(i-1);
end

figure;
plot(1:N+1, X, 'LineWidth', 1.5);
xlabel('$i$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$X_i$', 'Interpreter', 'latex', 'FontSize', 18);
title('Single Random Walk');

% Plot multiple random walks
M = 10; % Number of random walks
x = randn(M, N);
X = zeros(M, N+1); % Initialize positions
X(:,1) = 0; % Initial positions

for j = 1:M
    for i = 2:N+1
        X(j,i) = X(j,i-1) + x(j,i-1);
    end
end

figure;
plot(1:N+1, X', 'LineWidth', 1.5);
xlabel('$i$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$X_i$', 'Interpreter', 'latex', 'FontSize', 18);
title('Multiple Random Walks');

%% exercise 3
% This is the code to plot M random walks
M = 1000;
N = 99; % length of random walk
x = randn(M,N);
X = zeros(M, N+1); % initial position
for j = 1:M
    for i = 2:N+1
        X(j,i) = X(j,i-1) + x(j,i-1);
    end
end

% Construct histogram of final positions
bins = -40:2:40;
Xend = X(:,N+1);
[f, edges] = histcounts(Xend, bins);

% Plot first 100 trajectories
figure(1); clf;
subplot(1,2,1);
plot(1:N+1, X(1:100,:)','k');
xlabel('$i$','Interpreter','latex','fontsize',18);
ylabel('$X_i$','Interpreter','latex','fontsize',18);
ylim([-40 40]);

% Plot frequency distribution of final position
subplot(1,2,2);
barh(bins(1:end-1) + diff(bins)/2, f, 'histc');
set(gca, 'ytick', []);
xlabel('Frequency');
ylabel('Final Position');

% Mean and std of final position
mean_Xend = mean(Xend);
std_Xend = std(Xend);

disp(['Mean of final position: ', num2str(mean_Xend)]);
disp(['Standard deviation of final position: ', num2str(std_Xend)]);
