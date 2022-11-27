clear; 
clc;
close all
%%
syms thetav theta1 theta2 theta3 theta4 lv l1 l2 l3 l4 x_e y_e 
data = table2array(readtable('cycling_data.txt'));
% load('finalThetas.mat')
theta = 90 - data(:,1);
beta = deg2rad(data(:,2));

xpos_init = -0.1;
ypos_init = -0.7;


xp = xpos_init + 0.17*cosd(theta);
yp = ypos_init + 0.17*sind(theta);
pos = [xp,yp];

%%
syms xpos ypos
HomogeneousTrans = [cos(thetav),-sin(thetav),0,lv;
                   sin(thetav),cos(thetav),0,0;
                   0, 0, 1, 0;
                   0, 0, 0, 1]; %

H1 = subs(HomogeneousTrans,[thetav lv],[theta1 l1]);
H2 = subs(HomogeneousTrans,[thetav lv],[theta2 l2]);
H3 = subs(HomogeneousTrans,[thetav lv],[theta3 l3]);
H4 = subs(HomogeneousTrans,[thetav lv],[theta4 l4]);


Htotal = H1*H2*H3;
Htotal_v = subs(Htotal,[l1 l2 l3], [0 0.46,0.44]);

FinPosition_v=Htotal_v*[0;0;0;1];
eqn = solveThetas(FinPosition_v(1:2), [xpos;ypos]);
x = solve(eqn, [theta1, theta2]);

%%

newTotalThetas = zeros([2,2,361]);
for i = 1:361

    newTotalThetas(:,:,i) = [double(subs(x.theta1, [xpos,ypos],[pos(i,1), pos(i,2)])), ...
        double(subs(x.theta2, [xpos,ypos],[pos(i,1), pos(i,2)]))];
end

%% find consecutive points
finalThetas = zeros([361,2]);
finalThetas(1,:) = newTotalThetas(1,:,1);

for i = 1:360
    diff = abs(finalThetas(i,:) - newTotalThetas(:,:,i+1));
    [~,I] = min(diff);
    finalThetas(i+1,1) = newTotalThetas(I(1),1,i+1);
    finalThetas(i+1,2) = newTotalThetas(I(2),2,i+1);
end

%% Plot
fig = figure;

for i =1:361
        
    xlim([-.6,.6]);
    ylim([-1,.2])
%     t1 = newTotalThetas(1,1,i);
%     t2 = newTotalThetas(1,2,i);
    t1 = finalThetas(i,1);
    t2 = finalThetas(i,2);
    posH=[0 0];
    
    posK=subs(H1*H2*[0;0;0;1],[theta1 l1 l2],[t1,0, 0.46]);
    posK=posK(1:2).';
    
    posA=subs(H1*H2*H3*[0;0;0;1],[theta1,theta2 l1 l2 l3],[t1 t2 0,0.46,0.44]);
    posA=posA(1:2).';
    
    posT=[double(posA(1)) + 0.07*cos(beta(i)), double(posA(2)) + 0.07*sin(beta(i))]; %subs(H1*H2*H3*H4*[0;0;0;1],[theta1,theta2,theta3,l1,l2,l3,l4],[t1 t2 deg2rad(data(361-i,2))-t1-t2,0,0.46,0.44,0.07]);
    posT=posT(1:2).';
    
    hold on
    plot(posH(1),posH(2),'bo','MarkerSize',4)
    plot([posH(1);posK(1)],[posH(2);posK(2)],'k-o','MarkerSize',4)
    plot([posK(1);posA(1)],[posK(2);posA(2)],'m-o','MarkerSize',4)
    plot([posA(1);posT(1)],[posA(2);posT(2)],'r-o','MarkerSize',4)
    plot(xpos_init,ypos_init,'S','MarkerSize',8);
    drawnow;
    pause(0.01);
    clf(fig);
end

%% convert FT and FN to Fx and Fy 
Ft = -data(:,3);
Fn = -data(:,4);
Fx = Ft.*cos(beta) + Fn.*cos(pi/2 - beta);
Fy = Ft.*sin(beta) - Fn.*sin(pi/2 - beta);

%% variables to pass

combined = [finalThetas, beta, Fx, Fy, pos];
%theta1, theta2, theta3 (just defined from horizontal), fx, fy, x/y ankle position
%save('output.mat', 'combined')

%%

function [c,ceq] = NonLFunc(x,FinPositionFunc,finalpos)
    c=[];
    xval = x(3:5);
    ceq = FinPositionFunc(xval)-finalpos;
end

function F = solveThetas(FinPosition_v, pos)
F(1) = FinPosition_v(1) - pos(1);
F(2) = FinPosition_v(2) - pos(2);
end
