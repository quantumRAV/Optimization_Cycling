function totEqs = JointForceTorque_FileFcn(F_t1,F_t2,theta_h,theta_k,theta_a)
%JointForceTorque_FileFcn
%    totEqs = JointForceTorque_FileFcn(F_t1,F_t2,THETA_H,THETA_K,THETA_A)

%    This function was generated by the Symbolic Math Toolbox version 9.1.
%    27-Nov-2022 18:37:15

t2 = cos(theta_a);
t3 = cos(theta_h);
t4 = cos(theta_k);
t5 = sin(theta_a);
t6 = sin(theta_h);
t7 = sin(theta_k);
t8 = F_t2.*t3.*t4.*(1.1e+1./2.5e+1);
t9 = F_t1.*t3.*t7.*(1.1e+1./2.5e+1);
t10 = F_t1.*t4.*t6.*(1.1e+1./2.5e+1);
t11 = F_t2.*t6.*t7.*(1.1e+1./2.5e+1);
t15 = F_t2.*t2.*t3.*t4.*(7.0./1.0e+2);
t16 = F_t1.*t2.*t3.*t7.*(7.0./1.0e+2);
t17 = F_t1.*t2.*t4.*t6.*(7.0./1.0e+2);
t18 = F_t1.*t3.*t4.*t5.*(7.0./1.0e+2);
t19 = F_t2.*t2.*t6.*t7.*(7.0./1.0e+2);
t20 = F_t2.*t3.*t5.*t7.*(7.0./1.0e+2);
t21 = F_t2.*t4.*t5.*t6.*(7.0./1.0e+2);
t22 = F_t1.*t5.*t6.*t7.*(7.0./1.0e+2);
t29 = t3.*t4.*1.2236994e+1;
t30 = t6.*t7.*1.2236994e+1;
t32 = t2.*t6.*t7.*3.742515e-1;
t33 = t3.*t5.*t7.*3.742515e-1;
t34 = t4.*t5.*t6.*3.742515e-1;
t35 = t2.*t3.*t4.*3.742515e-1;
t12 = -t9;
t13 = -t10;
t14 = -t11;
t23 = -t16;
t24 = -t17;
t25 = -t18;
t26 = -t19;
t27 = -t20;
t28 = -t21;
t31 = -t30;
t36 = -t32;
t37 = -t33;
t38 = -t34;
totEqs = [F_t1;F_t2+1.06929e+1;F_t1;F_t2;F_t1;F_t2+1.185048e+2;F_t1;F_t2+4.49298e+1;F_t1;F_t2;t15+t22+t23+t24+t25+t26+t27+t28+t35+t36+t37+t38;t3.*3.7589958e+1+t8+t12+t13+t14+t15+t22+t23+t24+t25+t26+t27+t28+t29+t31+t35+t36+t37+t38+F_t2.*t3.*(2.3e+1./5.0e+1)-F_t1.*t6.*(2.3e+1./5.0e+1);t8+t12+t13+t14+t15+t22+t23+t24+t25+t26+t27+t28+t29+t31+t35+t36+t37+t38;0.0;0.0];