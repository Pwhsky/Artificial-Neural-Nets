%% Kohonen
clear
N = 1000; % Input size
M = 100; % Output size

% Generate the randomized input pattern and weights
x = 0.5*rand(N,2) + 0.5*(randi(3,N,1)>=[3 2]);
w = 2*rand(M,2)-1;

% Plot the input data
figure(1)
clf
scatter(x(:,1),x(:,2))

% Initial conditions for ordering phase
sigma_0 =100;
tau_sigma = 300;
eta_0 = 0.1;
T_order = 1e3;

% Ordering functions and neighborhood function
sigmaOrder = @(t) sigma_0* exp(-t/tau_sigma);
etaOrder = @(t) eta_0*exp(-t/tau_sigma);
Lambda = @(i_0,sig) exp(-(((1:M)'-i_0).^2)./(2*sig.^2));

% Convergence constants
sigmaConv = 0.9;
etaConv = 0.01;
T_conv = 2e4;

% Run the ordering network
for t = 1:T_order
    distance = x(randi(N),:)-w;
    [~,i_0] = min(sum(distance.^2,2));
    
    w = w + etaOrder(t).*Lambda(i_0,sigmaOrder(t)).*distance;
end

% Plot the result after ordering
figure(1)
clf
subplot(1,2,1)
plot(w(:,1),w(:,2),'k','Linewidth',2)
hold on
plot([0 0 1 1 0.5 0.5 0],[0 1 1 0.5 0.5 0 0],'k--','linewidth',2)
ylabel('$\xi_2$','interpreter','latex')
set(get(gca,'YLabel'),'Rotation',0,'FontSize',20,...
    'VerticalAlignment','middle','HorizontalAlignment','right');
xlabel('$\xi_1$','interpreter','latex')
set(get(gca,'xlabel'),'fontsize',20)
%%
% Run the convergence network
for t= 1:T_conv
    distance = x(randi(N),:)-w;
    [~,i_0] = min(sum(distance.^2,2));
    
    w = w + etaConv.*Lambda(i_0,sigmaConv).*distance;   
end

% Plot the result after convergence
subplot(1,2,2)
plot(w(:,1),w(:,2),'k','Linewidth',2)
hold on
plot([0 0 1 1 0.5 0.5 0],[0 1 1 0.5 0.5 0 0],'k--','linewidth',2)
ylabel('$\xi_2$','interpreter','latex')
set(get(gca,'YLabel'),'Rotation',0,'FontSize',20,...
    'VerticalAlignment','middle','HorizontalAlignment','right');
xlabel('$\xi_1$','interpreter','latex')
set(get(gca,'xlabel'),'fontsize',20)