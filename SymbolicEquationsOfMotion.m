clc 
clearvars
syms t thetaT 

lin_accel_O = sym('lin_accel_O',[3 1]); %acceleration of point O %change all these to symfunmatrix of variable t in order to perform time differentiation
alpha_OP = sym('alpha_OP',[3 1]); %angular acceleration of link OP
r_O_P = sym('r_O_P',[3 1]); %vector from O to point P on link OP
omega_OP = sym('omega_OP',[3 1]); %vector from O to point P on link OP
theta_OP = sym('theta_OP'); %Angle of link OP with frame B centered at point O, with respect to frame A.
l_OP = sym('l_OP'); %length of the link from O to P.

lin_accel = lin_accel_O + cross(alpha_OP,r_O_P) + cross(omega_OP,cross(omega_OP,r_O_P));

HomogenousTrans = [cos(thetaT),-sin(thetaT),0,l_OP;
                   sin(thetaT),cos(thetaT),0,0;
                   0, 0, 1, 0;
                   0, 0, 0, 1]; %homogenous transform transform vector from frame B to frame A when frame B rotates relative to Frame A with angle theta_OP

RotTrans =  [cos(thetaT),-sin(thetaT);
                   sin(thetaT),cos(thetaT)]; %rotation transform transform vector from frame B to frame A when frame B rotates relative to Frame A with angle theta_OP

frame_abbrev = {"h","k","a","t","p","c"}; %h=hip, k = knee, a = ankle, t = toe, p = pedal, c = center
link_length = {"h_k","k_a","a_t","t_p","p_c"}; % h_k = length from hip to knee, k_a = knee to ankle, a_t = ankle to toe, t_p = toe to pedal, p_c = pedal to crank center
numLinks = length(link_length);

thetaVec = sym(zeros(numLinks,3));
omegaVec = sym(zeros(numLinks,3));
alphaVec = sym(zeros(numLinks,3));
linAccelVec = sym(zeros(numLinks,3)); %need to initialize to zero because the first frame (in this case the hip) is assumed to not be moving.
linAccelCOMVec = sym(zeros(numLinks,3));
H_Mat = sym(repmat(eye(4),1,1,numLinks)); %homogenous transform to transform from frame{k} to the base frame at the origin
Rot_Mat = sym(repmat(eye(2),1,1,numLinks)); %rotation transform to transform from frame{k} to the base frame at the origin
rVec = sym(zeros(numLinks,3)); %vector from the base of the link to the end of the link
rCOMVec = sym(zeros(numLinks,3)); %vector from the base of the link to the COM of the link
I_mat = sym(repmat(zeros(3),1,1,numLinks)); %for inertia matrix

gVec = sym([0;-str2sym("g");0]);
sumFvec = sym(zeros(numLinks,3));
sumMoments = sym(zeros(numLinks,3));

symVec=sym(zeros(1,numLinks*2*3));
thetaSymVec_t=sym(zeros(1,numLinks)); %function of time
thetaSymVec=sym(zeros(1,numLinks)); %not a function of time

for k=1:numLinks
    thetaVec(k,3) = str2sym("theta"+"_"+frame_abbrev{k}+"(t)"); %this gives the angular position of frame{k} relative to the frame {k-1}
    thetaSymVec_t(k) = thetaVec(k,3);
    thetaSymVec(k) = str2sym("theta"+"_"+frame_abbrev{k});

    omegaVec(k,:) = diff(thetaVec(k,:),t,1);
    omegaVec(k,:) = sum(omegaVec(max(k-1,1):k,:),1); %omega for the frame is the sum of all the omegas of frame {k} and all previous frames
    alphaVec(k,:) = diff(omegaVec(k,:),t);   %simplyy differentiate because there are no additional acceleration terms resulting from movement of coordinate frame as the motion is assumed to happen in the plane
    I_mat(:,:,k) = sym("I_"+link_length{k},[3 3]);

    H_T = subs(HomogenousTrans,{thetaT, l_OP},[thetaVec(k,3),sym("L_"+link_length{k})]);
    H_Mat(:,:,k) = H_Mat(:,:,max(k-1,1))*H_T;
    rot_T = subs(RotTrans,{thetaT},[thetaVec(k,3)]);
    Rot_Mat(:,:,k) = Rot_Mat(:,:,max(k-1,1))*rot_T;
    

    rVec_H = Rot_Mat(:,:,k)*[H_T(1,4);0];
    rVec(k,:)=[rVec_H(:,1).',0];
    linAccelVec(k,:) = subs(lin_accel,[lin_accel_O,alpha_OP,omega_OP,r_O_P],[linAccelVec(max(k-1,1),:).',alphaVec(k,:).',omegaVec(k,:).',rVec(k,:).']);

    rVecCOM_H = Rot_Mat(:,:,k)*[H_T(1,4)/2;0];
    rCOMVec(k,:)=[rVecCOM_H(:,1).',0];
    linAccelCOMVec(k,:) = subs(lin_accel,[lin_accel_O,alpha_OP,omega_OP,r_O_P],[linAccelVec(max(k-1,1),:).',alphaVec(k,:).',omegaVec(k,:).',rCOMVec(k,:).']);

    

    %% Compute moment and sum of forces given that we've computed the linear and angular acceleration
    F_reaction_joint1 = sym("F_"+frame_abbrev{k},[3,1]) ;
    M_reaction_joint1 = sym("M_"+frame_abbrev{k},[3,1]) ;
    F_reaction_joint2 = -sym("F_"+frame_abbrev{k+1},[3,1]);
    M_reaction_joint2 = -sym("M_"+frame_abbrev{k+1},[3,1]) ;
    massSym =  sym("m_"+link_length{k});
    F_ext = sym("F_ext_"+frame_abbrev{k},[3,1]);
    M_ext = sym("M_ext_"+frame_abbrev{k},[3,1]);

    Fsum = F_reaction_joint1 + F_reaction_joint2 + massSym*gVec + F_ext ==linAccelVec(k,:).';
    Fsum(3)=0;%No forces in z direction because the motion is planar
    sumFvec(k,:) = Fsum.';

    Msum = I_mat(:,:,k)*alphaVec(k,:).' + cross(omegaVec(k,:).',I_mat(:,:,k)*omegaVec(k,:).') == cross(-rCOMVec(k,:).',F_reaction_joint1) + cross(rVec(k,:).'-rCOMVec(k,:).',F_reaction_joint2) + M_ext + M_reaction_joint1 + M_reaction_joint2;
    Msum(1)=0;%only moments around the z axis for planar motion
    Msum(2)=0;%only moments around the z axis for planar motion
    sumMoments(k,:) = Msum.';

    symVec((k-1)*6+1:(k-1)*6+3)=F_reaction_joint1.';
    symVec((k-1)*6+4:(k-1)*6+6)=M_reaction_joint1.';




end
totaleqs=[reshape(sumFvec,[15,1]);reshape(sumMoments,[15,1])];
%% get static equilibria
staticEqs = subs(totaleqs,thetaSymVec_t,thetaSymVec);
textStr="";
for jj=1:length(staticEqs)

    textStr = textStr+newline+"& "+latex(simplify(staticEqs(jj,1))) + " \\";

end


%writelines(textStr,"StaticEqEquations.txt")

%%print
%latex(simplify(staticEqs(30,1)))

%% equations for muscle forces and torques
createVarNameFcn = @(prefix,suffix) sprintf('%s_%s',prefix,suffix);

muscleSuffix ={'RF' 'IP' 'G' 'H' 'TA' 'GA'};
%muscleRadii = str2sym(cellfun(@(r)sym(createVarNameFcn('r',r),[6,1]),muscleSuffix,'UniformOutput',false));
muscleRadii = sym("r_",[3,6]);
muscleForce =  str2sym(cellfun(@(r)createVarNameFcn('F',r),muscleSuffix,'UniformOutput',false));
muscleActivation = str2sym(cellfun(@(r)createVarNameFcn('a',r),muscleSuffix,'UniformOutput',false));
MuscleTorques = sym({'M_h','M_k','M_t'});

maxMuscleForces = [1200*5, 1500*5, 3000*5, 3000*5, 2500, 3000];
MuscleTorqueEquations = MuscleTorques.'-[muscleRadii(1,1) muscleRadii(1,2) -muscleRadii(1,3) -muscleRadii(1,4) 0 0;muscleRadii(2,1) 0 0 -muscleRadii(2,4)  0 0;0 0 0 0 muscleRadii(3,5) -muscleRadii(3,6) ]*(muscleActivation.*muscleForce).';
MuscleTorqueEquations_P = subs(MuscleTorqueEquations,muscleRadii,repmat([0.081; 0.035;0.052],[1,6]));  %use same radii for now
MuscleTorqueEquations_P = subs(MuscleTorqueEquations_P,muscleForce, maxMuscleForces);

if true  %set to true if you want to regenerate the constraint functions for muscle optimization
    currdir = [pwd filesep]; % You might need to use currdir = pwd
    filename_cons = [currdir,'constraint_MuscleOptimization.m'];
    matlabFunction([],MuscleTorqueEquations_P,'file',filename_cons,'vars',{[MuscleTorques,muscleActivation]},'outputs',{'c','ceq'});
end

%% create table to input parameters
letterDictStruct.FirstLetter=containers.Map(["F","L","M","g","m","theta"],["Force","Length","Moment","gravity","mass","angle"]);
letterDictStruct.SecondLetter=containers.Map(["a","c","h","k","p","t","ext","1","2","3"],["ankle","crank center","hip","knee","pedal center","toe","external","x axis","y axis", "z axis"]);
svar_arr = symvar(staticEqs);
if false  %don't overwrite table
    createVariableTable("StaticVarTable.csv",svar_arr, letterDictStruct);
end

%% solve equations for the reaction forces and torques
parmTable = readtable("StaticVarTable.csv");

pTable = parmTable((parmTable.SolveSymbolically=="Y"),:);

StaticEqs_solved = solve(staticEqs,str2sym(pTable.VarName));

pTable2 = parmTable((parmTable.SolveSymbolically=="N" & ~isnan(parmTable.Value)),:); %for parameters
subEqs = subs(StaticEqs_solved,str2sym(pTable2.VarName).',pTable2.Value.');

%create matlab function from symbolic equations
paramOutputs = fieldnames(subEqs);  %names of the param outputs
totEqs = sym(zeros(length(paramOutputs),1));

for j= 1:length(paramOutputs)

    totEqs(j) = subEqs.(paramOutputs{j});


end

paramNames = parmTable.VarName(parmTable.SolveSymbolically=="P")
JointForceTorque_Fcn = matlabFunction(totEqs,'Vars',{'F_t1','F_t2','theta_h','theta_k','theta_a'})

if false %set to true if you want to regenerate a separate file for calculating joint forces and torques 
    JointForceTorque_Fcn = matlabFunction(totEqs,'Vars',{'F_t1','F_t2','theta_h','theta_k','theta_a'},'File','JointForceTorque_FileFcn')
end

%% read table of parameters and solve equations using the simplified structure subEqs


%convert matrix to table
KautzDataStruct = load("output.mat");
KautzData = KautzDataStruct.combined;
KautzData(:,3) = -KautzData(:,1) -KautzData(:,2)+KautzData(:,3);  %convert the ankle angles in the convention defined above
KautzData = array2table(KautzData,'VariableNames',{'theta_h','theta_k','theta_a','f_x','f_y','x','y'});
KautzData.f_x_negated=-KautzData.f_x; %need to negate because the negative sign is assumed in the EOMs.  Recall that these forces that are read in were converted to be in the global x and y frames
KautzData.f_y_negated=-KautzData.f_y; %need to negate because the negative sign is assumed in the EOMs.
parmTable = readtable("StaticVarTable.csv");

%get variables that are non-zero
pTable = parmTable((parmTable.SolveSymbolically=="P"),:); %for parameters

resultStruct=struct();
for vv=1:length(paramOutputs)

    resultStruct.(paramOutputs{vv}) = zeros(size(KautzData,1),1);

end

for vv=1:length(muscleSuffix)

    resultStruct.("Activation_"+muscleSuffix{vv}) = zeros(size(KautzData,1),1);

end


exitflags = zeros(size(KautzData,1),1)


for j=1:size(KautzData,1)
    
    %substitue F_t1 and F_t2 to the f_x and f_y_respectively.  Substitute
    %theta_h, theta_k,theta_a 
    resSolution = JointForceTorque_Fcn(KautzData.f_x_negated(j),KautzData.f_y_negated(j),KautzData.theta_h(j),KautzData.theta_k(j),KautzData.theta_a(j));

    for vv=1:length(paramOutputs)

        resultStruct.(paramOutputs{vv})(j) = resSolution(vv);

    end

    %%perform optimization to get the activations
    MuscleTorques = [resSolution(paramOutputs=="M_h3"),resSolution(paramOutputs=="M_k3"), resSolution(paramOutputs=="M_a3")] ;
    objFunc = @(x) sum(x.^2) + sum(x);
    x0=zeros(1,6);
    options = optimoptions('fmincon','Algorithm','sqp');
    lb=zeros(1,6)+0.1;
    ub=ones(1,6);
    [xval_activation,fval,exitflags(j),output,lambda,grad]=fmincon(objFunc,x0,[],[],[],[],lb,ub,@(x)constraint_MuscleOptimization([MuscleTorques,x]),options);
    for vv=1:length(muscleSuffix)

        resultStruct.("Activation_"+muscleSuffix{vv})(j) = xval_activation(vv);

    end

    %check exit flag


end

disp("Finished")
resultTable = struct2table(resultStruct);

totalTable = [KautzData,resultTable];
writetable(totalTable,"results_w_Activation.xlsx")



%% plots

thighL = parmTable.Value(parmTable.VarName=="L_h_k");
shankL = parmTable.Value(parmTable.VarName=="L_k_a");
footL =parmTable.Value(parmTable.VarName=="L_a_t");

Ho_func = matlabFunction(HomogenousTrans,'vars',[thetaT,l_OP]);



close all
numMuscles = 6;
fig=figure();
set(gcf, 'Units', 'normalized');
set(gcf, 'Position', [0 0.1 0.8 0.8]);
set(gcf,'color','w');
set(0, 'DefaultAxesFontName', 'Arial')
tiledlayout(numMuscles,2);

leg_ax = nexttile([numMuscles,1]);
set(leg_ax,'XColor', 'none','YColor','none');

axColl = {};
animatedLineCol = {}

MuscleAVid = VideoWriter('MuscleActivation'); %open video file
MuscleAVid.FrameRate = 30;  
open(MuscleAVid)


muscleNames = {'Rectus Femoris (RF)' 'Iliopsoas (IP)' 'Gluteals (G)' 'Hamstrings (H)' 'Tibialis Anterior (TA)' 'Gastrocnemius (GA)'};

for vv = 1:numMuscles
    axCol1{vv} = nexttile;
    set(axCol1{vv},'FontSize',12)
    animatedLineCol{vv} = animatedline(axCol1{vv});
    xlim(axCol1{vv},[0,360]);
    ylim(axCol1{vv},[0,max(totalTable.(sprintf("Activation_%s",muscleSuffix{vv})))*maxMuscleForces(vv)]);
    %ylim(axCol1{vv},[0,inf]);
    title(muscleNames{vv});
    ylabel('Force (N)','FontSize',14,'FontWeight','bold');


end


aLine_thigh = animatedline(leg_ax,'LineWidth',3,'Marker','o','MarkerFaceColor','r','MarkerSize',9);
aLine_shank = animatedline(leg_ax,'LineWidth',3,'Marker','o','MarkerFaceColor','b','MarkerSize',9);
aLine_foot = animatedline(leg_ax,'LineWidth',3,'Marker','o','MarkerFaceColor','g','MarkerSize',9);
xlim(leg_ax,[-.6,.6]);
ylim(leg_ax,[-1,.2]);

xlabel(axCol1{6},'Crank Angle of Bicycle (deg)', 'FontSize',14,'FontWeight','bold');


angles=[0:360];


for j=1:length(angles)
    for v=1:length(muscleSuffix)
        
        activationV = totalTable.(sprintf("Activation_%s",muscleSuffix{v}));
        %plot(axCol1{v},angles(1:j),activationV(1:j));
        addpoints(animatedLineCol{v},angles(j),activationV(j).*maxMuscleForces(v)   );
    end


    %plot cyclist
    posH    = [0;0];
    posK = Ho_func(totalTable.theta_h(j),0)*Ho_func(totalTable.theta_k(j),thighL)*[0;0;0;1];
    posA = Ho_func(totalTable.theta_h(j),0)*Ho_func(totalTable.theta_k(j),thighL)*Ho_func(totalTable.theta_a(j),shankL)*[0;0;0;1];
    posT = Ho_func(totalTable.theta_h(j),0)*Ho_func(totalTable.theta_k(j),thighL)*Ho_func(totalTable.theta_a(j),shankL)*Ho_func(0,footL)*[0;0;0;1];

    clearpoints(aLine_thigh);
    clearpoints(aLine_shank);
    clearpoints(aLine_foot);

    addpoints(aLine_thigh,[posH(1);posK(1)],[posH(2);posK(2)]);
    addpoints(aLine_shank,[posK(1);posA(1)],[posK(2);posA(2)]);
    addpoints(aLine_foot,[posA(1);posT(1)],[posA(2);posT(2)]);

    drawnow limitrate;
    pause(0.0027);
    frame = getframe(gcf); %get frame
    writeVideo(MuscleAVid, frame);

    %clf(fig);

end

close(MuscleAVid);



%% Plot torque

numMuscles = 6;
fig=figure();
set(gcf, 'Units', 'normalized');
set(gcf, 'Position', [0 0.1 0.8 0.8]);
set(gcf,'color','w');
set(0, 'DefaultAxesFontName', 'Arial')
tiledlayout(numMuscles,2);

ax=nexttile(1);
plot(ax,angles,totalTable.M_h3);
ylabel('(Nm)','FontSize',12,'FontWeight','bold');
title('Hip Torque','FontSize',12,'FontWeight','bold');

ax=nexttile(3);
plot(ax,angles,totalTable.M_k3);
ylabel('(Nm)','FontSize',12,'FontWeight','bold');
title('Knee Torque','FontSize',12,'FontWeight','bold');

ax=nexttile(5);
plot(ax,angles,totalTable.M_a3);
ylabel('(Nm)','FontSize',12,'FontWeight','bold');
title('Ankle Torque','FontSize',12,'FontWeight','bold');

ax=nexttile(7);
plot(ax,angles,totalTable.f_x);
ylabel('(N)','FontSize',12,'FontWeight','bold');
title('Foot x-force','FontSize',12,'FontWeight','bold');

ax=nexttile(9);
plot(ax,angles,totalTable.f_y);
ylabel('(N)','FontSize',12,'FontWeight','bold');
title('Foot y-force (N)','FontSize',12,'FontWeight','bold');
xlabel('Crank Angle of Bicycle (deg)', 'FontSize',12,'FontWeight','bold');


axColl = {};




muscleNames = {'Rectus Femoris (RF)' 'Iliopsoas (IP)' 'Gluteals (G)' 'Hamstrings (H)' 'Tibialis Anterior (TA)' 'Gastrocnemius (GA)'};

for vv = 1:numMuscles
    axCol1{vv} = nexttile((vv)*(2));
    
    
    
    activationV = totalTable.(sprintf("Activation_%s",muscleSuffix{vv}));
        %plot(axCol1{v},angles(1:j),activationV(1:j));
    plot(axCol1{vv},angles,activationV,'k-');
    ylabel(axCol1{vv} ,'Activation','FontSize',14,'FontWeight','bold');
    title(muscleNames{vv});

    set(axCol1{vv},'FontSize',12)
    xlim(axCol1{vv},[0,360]);
    ylim(axCol1{vv},[0,1]);


end
xlabel(axCol1{6},'Crank Angle of Bicycle (deg)', 'FontSize',14,'FontWeight','bold');

%% plot the activations side by side

PrilutskyData = readtable("PrilutskyData_AdjustedActivation.csv");

numMuscles = 6;
fig=figure();
set(gcf, 'Units', 'normalized');
set(gcf, 'Position', [0 0.1 0.8 0.8]);
set(gcf,'color','w');
set(0, 'DefaultAxesFontName', 'Arial')
tiledlayout(numMuscles,2);
muscleNames = {'Rectus Femoris (RF)' 'Iliopsoas (IP)' 'Gluteals (G)' 'Hamstrings (H)' 'Tibialis Anterior (TA)' 'Gastrocnemius (GA)'};

PrilutskyPos =[1,5,7,9,11];
PrilutskyMuscles={"RF","GLM","HA","TA","GA"};
numMuscles_Pri = 5;

axColl = {};

for vv = 1:numMuscles_Pri
    axCol1{vv} = nexttile(PrilutskyPos(vv));
    pT = PrilutskyData(strcmp(PrilutskyData.Muscle,PrilutskyMuscles{vv}),:);
    activationV_P =interp1(pT.angle,pT.AdjustedActivation,angles,'linear') ;
        %plot(axCol1{v},angles(1:j),activationV(1:j));
    plot(axCol1{vv},angles,activationV_P/100,'k-');
    ylabel(axCol1{vv} ,'Activation','FontSize',14,'FontWeight','bold');
    title(PrilutskyMuscles{vv});

    set(axCol1{vv},'FontSize',12)
    xlim(axCol1{vv},[0,360]);
    ylim(axCol1{vv},[0,max(activationV_P/100)*1.1]);


end
xlabel(axCol1{5},'Crank Angle of Bicycle (deg)', 'FontSize',14,'FontWeight','bold');




axColl = {};

for vv = 1:numMuscles
    axCol1{vv} = nexttile((vv)*(2));
    
    
    
    activationV = totalTable.(sprintf("Activation_%s",muscleSuffix{vv}));
        %plot(axCol1{v},angles(1:j),activationV(1:j));
    plot(axCol1{vv},angles,activationV,'k-');
    ylabel(axCol1{vv} ,'Activation','FontSize',14,'FontWeight','bold');
    title(muscleNames{vv});

    set(axCol1{vv},'FontSize',12)
    xlim(axCol1{vv},[0,360]);
    ylim(axCol1{vv},[0,max(activationV)]*1.1);


end
xlabel(axCol1{6},'Crank Angle of Bicycle (deg)', 'FontSize',14,'FontWeight','bold');



function[outTable] = createVariableTable(fnamestr,sym_vars, letterDictStruct)

    numVars = length(sym_vars);
    out_s.Var=sym(zeros(numVars,1));
    out_s.VarName=cell([numVars,1]);
    out_s.definition=cell([numVars,1]);
    out_s.Value=NaN([numVars,1]);
    for j=1:length(sym_vars)
        s_var =(sym_vars(j));
        out_s.Var(j)=s_var;
        out_s.VarName{j}=string(s_var);
        s_var=string(s_var); %convert to string
        spVar = split(s_var,"_");
        DefS="";
        for vv=1:length(spVar)
            if vv ==1
                DefS = DefS+letterDictStruct.FirstLetter(spVar(vv))+":";
            else
                %get number which represents component at end
                returnv=regexp(spVar(vv),".*(\d)","tokens");
                if ~isempty(returnv)
                    St1=strsplit(spVar(vv),returnv{1});
                    DefS = DefS + letterDictStruct.SecondLetter(St1{1})+"_";
                    DefS = DefS + letterDictStruct.SecondLetter(returnv{1});
                else
                    DefS = DefS+ letterDictStruct.SecondLetter(spVar(vv)) +"_";
                end
                
                

            end

        end

        out_s.definition{j} = DefS;

    end

    outTable = struct2table(out_s);
    writetable(outTable,fnamestr,'Delimiter',','); 

end
