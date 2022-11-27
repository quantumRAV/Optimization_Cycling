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


writelines(textStr,"StaticEqEquations.txt")

%%print
%latex(simplify(staticEqs(30,1)))


%% solve system of equations for the reaction forces:

%sumv=reshape(sumFvec,[15 1]);
%qg=solve(sumv(1:10).',str2sym({'F_a1','F_k2','F_k1','F_a2','F_h1','F_h2','F_t1','F_t2','F_p1','F_p2','F_c1','F_c2'}));  %solve system of equations for the unknown equations


totaleqs=[reshape(sumFvec,[15,1]);reshape(sumMoments,[15,1])];
qg2=solve(totaleqs.',symVec);
%sol=subs(totaleqs,qg2); %substitute out the reaction forces and moments.



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

%% read table of parameters
parmTable = readtable("StaticVarTable.csv");

%get variables that are non-zero
pTable = parmTable(~isnan(parmTable.Value),:);

newStaticEqns = simplify(subs(staticEqs,str2sym(pTable.VarName).',pTable.Value.'));
newStaticEqns = newStaticEqns(find(newStaticEqns~=symtrue)); %get rid of redundant equations

%% solve static equations for reaction forces and torques, give force at the feet
nSE = newStaticEqns(find(newStaticEqns~=0));  %only get the non zero equations after simplification
returnSolution =(solve(nSE));

%returnSolution = vpa(returnSolution,5);


%% read table of parameters and solve equations

%convert matrix to table
KautzDataStruct = load("output.mat");
KautzData = KautzDataStruct.combined;
KautzData(3,:) = -KautzData(1,:) -KautzData(2,:)+KautzData(3,:);  %convert the ankle angles in the convention defined above
KautzData = array2table(KautzData,'VariableNames',{'theta_h','theta_k','theta_a','f_x','f_y','x','y'});

parmTable = readtable("StaticVarTable.csv");

%get variables that are non-zero
pTable = parmTable(~isnan(parmTable.Value),:);

if false
    for j=1:size(KautzData,1)
    
        %substitue F_t1 and F_t2 to the f_x and f_y_respectively.  Substitute
        %theta_h, theta_k,theta_a to be 
        pTable(find(pTable.VarName=="F_t1"|pTable.VarName=="F_t2"|pTable.VarName=="theta_a"|pTable.VarName=="theta_h"|pTable.VarName=="theta_k"),'Value')={KautzData.f_x(j);KautzData.f_y(j);KautzData.theta_a(j);KautzData.theta_h(j);KautzData.theta_k(j)};
    
        newStaticEqns = simplify(subs(staticEqs,str2sym(pTable.VarName).',pTable.Value.'));
        newStaticEqns = newStaticEqns(find(newStaticEqns~=symtrue)); %get rid of redundant equations
        
        %% solve static equations for reaction forces and torques, give force at the feet
        nSE = newStaticEqns(find(newStaticEqns~=0));  %only get the non zero equations after simplification
        returnSolution =(solve(nSE));
    
    
    end
end

%% read table of parameters and solve equations using the simplified structure subEqs

KautzDataStruct = load("output.mat");
KautzData = KautzDataStruct.combined;
KautzData(:,3) = -KautzData(:,1) -KautzData(:,2)+KautzData(:,3);  %convert the ankle angles in the convention defined above
KautzData = array2table(KautzData,'VariableNames',{'theta_h','theta_k','theta_a','f_x','f_y','x','y'});
KautzData.f_x=-KautzData.f_x; %need to negate because the negative sign is assumed in the EOMs.  Recall that these forces that are read in were converted to be in the global x and y frames
KautzData.f_y=-KautzData.f_y; %need to negate because the negative sign is assumed in the EOMs.
parmTable = readtable("StaticVarTable.csv");

%get variables that are non-zero
pTable = parmTable((parmTable.SolveSymbolically=="P"),:); %for parameters

resultStruct=struct();
fstructn = fieldnames(subEqs);
for vv=1:length(fstructn)

    resultStruct.(fstructn{vv}) = zeros(size(KautzData,1),1);

end


for j=1:size(KautzData,1)

    %substitue F_t1 and F_t2 to the f_x and f_y_respectively.  Substitute
    %theta_h, theta_k,theta_a to be 
    pTable(find(pTable.VarName=="F_t1"|pTable.VarName=="F_t2"|pTable.VarName=="theta_a"|pTable.VarName=="theta_h"|pTable.VarName=="theta_k"),'Value')={KautzData.f_x(j);KautzData.f_y(j);KautzData.theta_a(j);KautzData.theta_h(j);KautzData.theta_k(j)};

    newStaticEqns = subs(subEqs,str2sym(pTable.VarName).',pTable.Value.');
    %newStaticEqns = newStaticEqns(find(newStaticEqns~=symtrue)); %get rid of redundant equations
    
    %% solve static equations for reaction forces and torques, give force at the feet
    %nSE = newStaticEqns(find(newStaticEqns~=0));  %only get the non zero equations after simplification
    returnSolution =newStaticEqns;

    for vv=1:length(fstructn)

        resultStruct.(fstructn{vv})(j) = returnSolution.(fstructn{vv});

    end



end

disp("Finished")
resultTable = struct2table(resultStruct);

totalTable = [KautzData,resultTable];

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
