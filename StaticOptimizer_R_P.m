frmax = 1200*5; % Rectus femoris 
fimax = 5*1500; % Iliopsoas
fgmax = 3000*5; % gluteals
fhmax = 3000*5; % Hamstring Fhmax rh
ftmax = 2500; % Tibialis Anterior
fgamax = 3000; % Gastrocneius
rh = 0.081; %radius of the hip joint 
rk = 0.035; %radius of knee joint
rankle = 0.052; %radius of the ankle joint

%% Setting up the problem
% x1 = activation Rectus femoris
% x2 = activation Iliopsoas
% x3 = activation gluteals
% x4 = activation hamstring
% x5 = activation Tibialis Anterior
% x6 = activation gastrocneius
load('ForceTorque.mat')
objective = @(x) x(1)^2 + x(2)^2 + x(3)^2 + x(4)^2 + x(5)^2 + x(6)^2;
% initial guess
x0 = [0.2 0.2 0.2 0.2 0.2 0.2];
% variable bounds
lb = [0.1 0.1 0.1 0.1 0.1 0.11];
ub = [1 1 1 1 1 1];
% show initial objective
disp(['Initial Objective: ' num2str(objective(x0))])

%% linear constraints
% Mh3 = - Fhmax*x4*rh - Fgmax*x3*rh + Fimax*x2*rh + Frmax*x1*rh
% Mk3 = -Fhmax*x4*rk + Frmax*x1*rk
% Ma3 = -Fgamax*x6*ra + Ftmax*x5*ra
A = [];
b = [];
Aeq = [frmax*rh fimax*rh -fgmax*rh -fhmax*rh 0 0; frmax*rk 0 0 -fhmax*rk 0 0; 0 0 0 0 ftmax*rankle -fgmax*rankle];
% nonlinear constraints
nonlincon = @nlcon;
%options
options = optimoptions(@fmincon,'Algorithm','interior-point');
%moments
mh3 = totalTable.("M_h3");
mk3 = totalTable.("M_k3");
ma3 = totalTable.("M_a3");

%% Managing the output table
x1 = [];
x2 = [];
x3 = [];
x4 = [];
x5 = [];
x6 = [];
totalActivation = [];
eflag = [];

%% Running the optimizer

for i = 1:1:361
    beq = [mh3(i); mk3(i); ma3(i)];

    % optimize with fmincon
    %[X,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN]
    % = fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
    [X, FVAL, EXITFLAG, OUTPUT, LAMBDA, GRAD, HESSIAN] = fmincon(objective,x0,A,b,Aeq,beq,lb,ub,nonlincon, options);
    x1(end+1) = X(1);
    x2(end+1) = X(2);
    x3(end+1) = X(3);
    x4(end+1) = X(4);
    x5(end+1) = X(5);
    x6(end+1) = X(6);
    totalActivation(end+1) = FVAL;
    eflag(end+1) = EXITFLAG;
end

%% Displaying the results
% show final objective
X = [x1; x2; x3; x4; x5; x6];
figure();
tiledlayout(6,1)

for k = 1:6
    ax=nexttile;
    plot(ax,0:360,X(k,:))

end

figure();
tiledlayout(3,1)
ax=nexttile
plot(ax,0:360,totalTable.M_h3)
ax=nexttile
plot(ax,0:360,totalTable.M_k3)
ax=nexttile
plot(ax,0:360,totalTable.M_a3)


figure();
tiledlayout(2,1)
ax=nexttile
plot(ax,0:360,totalTable.f_x)
ax=nexttile
plot(ax,0:360,totalTable.f_y)

figure();
plot(0:360,totalActivation)

function [c,ceq] = nlcon(x)
  c = [];
  ceq = [];

end