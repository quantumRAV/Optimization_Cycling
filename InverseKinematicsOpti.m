close all
syms thetav theta1 theta2 theta3 theta4 lv l1 l2 l3 l4 x_e y_e 

xArr = sym('x',[5 1]);


HomogeneousTrans = [cos(thetav),-sin(thetav),0,lv;
                   sin(thetav),cos(thetav),0,0;
                   0, 0, 1, 0;
                   0, 0, 0, 1]; %

H1 = subs(HomogeneousTrans,[thetav lv],[theta1 l1]);
H2 = subs(HomogeneousTrans,[thetav lv],[theta2 l2]);
H3 = subs(HomogeneousTrans,[thetav lv],[theta3 l3]);
H4 = subs(HomogeneousTrans,[thetav lv],[theta4 l4]);

Htotal = H1*H2*H3*H4;

Htotal_v = subs(Htotal,[l1 l2 l3 l4 theta4], [0 0.46,0.44,0.07, 0]);

FinPosition = simplify(Htotal*[x_e;y_e;0;1]); %symbolic
FinPosition_v=Htotal_v*[0;0;0;1];
FinPosition_v=FinPosition_v(1:2);

costFuncSym = xArr(1)^2 + xArr(2)^2;

FinPositionFunc = matlabFunction(FinPosition_v,'Vars',{[theta1,theta2,theta3]});
costFunc = matlabFunction(costFuncSym,'Vars',{[xArr(1),xArr(2),xArr(3),xArr(4),xArr(5)]});

xpos =0.45;
ypos = -0.6; %x=0.5, y = -0.8312: straight leg

lb = [xpos-0.001 ypos-0.001 pi pi 0];
ub = [xpos+0.001 ypos+0.001 2*pi 2*pi pi];

options = optimoptions("fmincon",...
    "Algorithm","interior-point",...
    "EnableFeasibilityMode",true,...
    "SubproblemAlgorithm","cg",'MaxFunctionEvaluations',3e4);
gg= fmincon(costFunc,[xpos, ypos, 3*pi/2, 3*pi/2, 0],[],[],[],[],lb,ub,@(x)NonLFunc(x,FinPositionFunc,[xpos;ypos]),options); %sensitive to initial condition.  3pi/2 3pi/2 and 0 is better than pi,pi and 0 for reaching more extreme points

%% Plot
posH=[0 0];

posK=subs(H1*H2*[0;0;0;1],[theta1 l1 l2],[gg(3),0, 0.46]);
posK=posK(1:2).';

posA=subs(H1*H2*H3*[0;0;0;1],[theta1,theta2 l1 l2 l3],[gg(3) gg(4) 0,0.46,0.44]);
posA=posA(1:2).';

posT=subs(H1*H2*H3*H4*[0;0;0;1],[theta1,theta2,theta3,l1,l2,l3,l4],[gg(3) gg(4) gg(5),0,0.46,0.44,0.07]);
posT=posT(1:2).';

figure()
hold on
plot(posH(1),posH(2),'bo','MarkerSize',4)
plot([posH(1);posK(1)],[posH(2);posK(2)],'k-o','MarkerSize',4)
plot([posK(1);posA(1)],[posK(2);posA(2)],'m-o','MarkerSize',4)
plot([posA(1);posT(1)],[posA(2);posT(2)],'r-o','MarkerSize',4)
plot(xpos,ypos,'S','MarkerSize',8);

xlim([-0.1,1]);
ylim([-1,0.11])

function [c,ceq] = NonLFunc(x,FinPositionFunc,finalpos)
    c=[];
    xval = x(3:5);
    ceq = FinPositionFunc(xval)-finalpos;
%    gg = ceq;
%    posH=[0 0];

% posK=subs(H1*H2*[0;0;0;1],[theta1 l1 l2],[gg(3),0, 0.46]);
% posK=posK(1:2).';
% 
% posA=subs(H1*H2*H3*[0;0;0;1],[theta1,theta2 l1 l2 l3],[gg(3) gg(4) 0,0.46,0.44]);
% posA=posA(1:2).';
% 
% posT=subs(H1*H2*H3*H4*[0;0;0;1],[theta1,theta2,theta3,l1,l2,l3,l4],[gg(3) gg(4) gg(5),0,0.46,0.44,0.07]);
% posT=posT(1:2).';
% 
%     figure(1)
%     hold on
%     plot(posH(1),posH(2),'bo','MarkerSize',4)
%     plot([posH(1);posK(1)],[posH(2);posK(2)],'k-o','MarkerSize',4)
%     plot([posK(1);posA(1)],[posK(2);posA(2)],'m-o','MarkerSize',4)
%     plot([posA(1);posT(1)],[posA(2);posT(2)],'r-o','MarkerSize',4)
%     plot(0,0.9,'S','MarkerSize',8);
    

end


