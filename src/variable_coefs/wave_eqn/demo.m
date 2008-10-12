%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This MATLAB script demonstrates the use of 
% solveWaveEqn1dKPY() and solveWaveEqn1dKPY_OTS()
%  
% Kevin T. Chu
% 2008 August
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
format long
format compact

% simulation parameters
debug_on = 0;

% time integration parameters
% NOTE: t_init = 0.0
t_final = 10.0;
t_final = 1.0;

% set dx and dt
x_lo = -1.0;
x_hi =  1.0;
N = 200;
dx = (x_hi-x_lo)/N;
dt_KPY = dx/4;
dt_VarMesh = dx/4;

% solve 1d wave equation using KPY time integration on transformed domain
% with OTS
[u_KPY_transformed_OTS, u_exact_transformed, x_transformed] = ...
  solveWaveEqn1dKPY_Transformed_OTS(N, t_final, debug_on);

% solve 1d wave equation using KPY time integration on transformed domain
% without OTS
[u_KPY_transformed, u_exact_transformed, x_transformed] = ...
  solveWaveEqn1dKPY_Transformed(N, dt_KPY, t_final, debug_on);

% solve 1d wave equation using KPY time integration on variable mesh
% without OTS
[u_KPY_VarMesh, u_exact_VarMesh, x_VarMesh] = ...
  solveWaveEqn1dKPY_VarMesh(N, dt_VarMesh, t_final, debug_on);

% solve 1d wave equation using KPY time integration on variable mesh
% with OTS
[u_KPY_VarMesh_OTS, u_exact_VarMesh, x_VarMesh] = ...
  solveWaveEqn1dKPY_VarMesh_OTS(N, t_final, debug_on);

% solve 1d wave equation using standard KPY time integration 
[u_KPY, u_exact, x] = solveWaveEqn1dKPY(N, dt_KPY, t_final, debug_on);

% compute error
err_KPY_transformed_OTS = u_KPY_transformed_OTS-u_exact_transformed;
err_L_inf_KPY_transformed_OTS = norm(err_KPY_transformed_OTS,'inf')
err_KPY_transformed = u_KPY_transformed-u_exact_transformed;
err_L_inf_KPY_transformed = norm(err_KPY_transformed,'inf')
err_KPY_VarMesh = u_KPY_VarMesh-u_exact_VarMesh;
err_KPY_VarMesh_OTS = u_KPY_VarMesh_OTS-u_exact_VarMesh;
err_L_inf_VarMesh_OTS = norm(err_KPY_VarMesh_OTS,'inf')
err_L_inf_VarMesh = norm(err_KPY_VarMesh,'inf')
err_KPY = u_KPY-u_exact;
err_L_inf_KPY = norm(err_KPY,'inf')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% KPY-Transformed OTS
figure(1); clf;
plot(x_transformed,u_KPY_transformed_OTS,'bo')
hold on;
plot(x_transformed,u_exact_transformed,'r')
title('KPY-OTS Solution')

figure(2); clf;
plot(x_transformed,err_KPY_transformed_OTS);
title('Error in KPY-OTS Solution')

% KPY-Transformed no OTS
figure(3); clf;
plot(x_transformed,u_KPY_transformed,'bo')
hold on;
plot(x_transformed,u_exact_transformed,'r')
title('KPY-transformed Solution')

figure(4); clf;
plot(x_transformed,err_KPY_transformed);
title('Error in KPY-transformed Solution')

% KPY-VariableMesh with OTS
figure(5); clf;
plot(x_VarMesh,u_KPY_VarMesh_OTS,'bo')
hold on;
plot(x_VarMesh,u_exact_VarMesh,'r')
title('KPY-VarMesh Solution')

figure(6); clf;
plot(x_VarMesh,err_KPY_VarMesh_OTS);
title('Error in KPY-VarMesh Solution')

% KPY-VariableMesh no OTS
figure(7); clf;
plot(x_VarMesh,u_KPY_VarMesh,'bo')
hold on;
plot(x_VarMesh,u_exact_VarMesh,'r')
title('KPY-VarMesh Solution')

figure(8); clf;
plot(x_VarMesh,err_KPY_VarMesh);
title('Error in KPY-VarMesh Solution')

% Standard KPY 
figure(9); clf;
plot(x,u_KPY,'bo')
hold on;
plot(x,u_exact,'r')
title('KPY Solution')

figure(10); clf;
plot(x,err_KPY);
title('Error in KPY Solution')

