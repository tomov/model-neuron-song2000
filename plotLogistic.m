logistic = @(Vthresh, lambda, V) 1 ./ (1 + exp(- (V - Vthresh) / lambda));

V = -60:0.01:-48;
Vthresh = -54;
lambda = 0.01;

figure;
plot(V, logistic(Vthresh, lambda, V));
ylabel('P_{spike}');
xlabel('V (mV)');
text(-59, 0.6, '\lambda = 0.01');
xlim([-60 -48]);