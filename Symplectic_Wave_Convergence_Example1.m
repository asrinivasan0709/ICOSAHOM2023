clc; clear; close all;
% Symplectic Integrator, Example
% ASrinivasan, 23Dec2022
% 1D Wave eqn
addpath('...\mole-master\mole_MATLAB')

NE = 100*[4, 8, 16, 32, 64, 128]';  

for i = 1:size(NE, 1)
    fprintf('Executing iteration %i of %i ... \n', i, size(NE, 1)); 

    NElem = NE(i, 1);
    xmin = -30; xmax = 30; % +/- 5
    dh = (xmax-xmin)/NElem; 
    tEnd = 24;   % 200
    CFL = 0.5; cV = 1; 
    dt =  CFL*dh/cV; 
    
    % init val
    xGrid = [xmin xmin+dh/2:dh:xmax-dh/2 xmax]';
    xNod = [xmin:dh:xmax]';
    muu = 0.; sg = 0.1;
%     u0 = 1/sqrt(sg*2*pi)*exp(-(xGrid-muu).^2/2/sg);
%     u0 = sech(15*(xGrid + 1)).^2; 
    u0 = exp(-100*(xGrid - 0.5).^2); % works for convergence
    v0 = zeros(size(xGrid,1),1);
    u0 = [u0;v0]; 
   

%     RRK, 4th order
    fprintf(' ... Executing RRK ... \n '); 
    uRRK = sparse(2*(NElem+2), tEnd/dt); 

    tic;
    tStart1 = 0; tEnd1 = tEnd/2;
    [uRRK1, tRRK1, gamRRK1, eRRK1] = RRK(1, NElem, dh, dt, tStart1, tEnd1, u0, cV); 

    tStart2 = tRRK1(end); tEnd2 = tEnd;
    u02 = uRRK1(:, end);     
    [uRRK2, tRRK2, gamRRK2, eRRK2] = RRK(1, NElem, dh, dt, tStart2, tEnd2, u02, cV); 
    
    uRRK = [uRRK2(:, end-5:end)]; 
    tRRK = [tRRK1; tRRK2];
    eRRK = [eRRK1; eRRK2]; 

    TimeRRK(i, 1) = toc;

    fEvalRRK(i, 1) = size(tRRK,1)*4*NElem;
    clear uRRK1 tRRK1 gamRRK1 eRRK1 uRRK2 tRRK2 gamRRK2 eRRK2; 

%     RK4, 4th order
    fprintf(' ... Executing RK4 ... \n '); 
    uRK4 = sparse(2*(NElem+2), tEnd/dt); 

    tic;
    tStart1 = 0; tEnd1 = tEnd/2;
    [uRK41, tRK41, gamRK41, eRK41] = RRK(0, NElem, dh, dt, tStart1, tEnd1, u0, cV);

    tStart2 = tRK41(end); tEnd2 = tEnd;
    u02 = uRK41(:, end); 
    [uRK42, tRK42, gamRK42, eRK42] = RRK(0, NElem, dh, dt, tStart2, tEnd2, u02, cV);
    
    uRK4 = [uRK42(:, end-5:end)]; 
    tRK4 = [tRK41; tRK42];
%     gamRK4 = [gamRK41; gamRK42]; 
    eRK4 = [eRK41; eRK42]; 

    TimeRK4(i, 1) = toc;
    
    fEvalRK4(i, 1) = size(tRK4,1)*4*NElem;
    clear uRK41 tRK41 gamRK41 eRK41 uRK42 tRK42 gamRK42 eRK42; 

% Forest Ruth
    fprintf(' ... Executing Forest Ruth ... \n '); 
    uFRuth = sparse(2*(NElem+2), tEnd/dt); 
    tic;
    tStart1 = 0; tEnd1 = tEnd/2;
    [uFRuth1, EFRuth1, tFRuth1] = FRuth(NElem,dh,dt,tStart1,tEnd1,u0, cV);

    tStart2 = tFRuth1(end); tEnd2 = tEnd;
    u02 = uFRuth1(:, end); 
    [uFRuth2, EFRuth2, tFRuth2] = FRuth(NElem,dh,dt,tStart2,tEnd2,u02, cV);

    uFRuth = [uFRuth2(:, end-5:end)]; 
    tFRuth = [tFRuth1; tFRuth2]; EFRuth = [EFRuth1; EFRuth2];

    TimeFRuth(i, 1) = toc; 
    
    fEvalFRuth(i, 1) = size(tFRuth,1)*3*NElem;
    clear uFRuth1 EFRuth1 tFRuth1 uFRuth2 EFRuth2 tFRuth2; 

%     PEFRL
    fprintf(' ... Executing PEFRL ... \n '); 
    uPEFRL = sparse(2*(NElem+2), tEnd/dt); 

    tic;
    tStart1 = 0; tEnd1 = tEnd/2;
    [uPEFRL1, EPEFRL1, tPEFRL1] = PEFRL(NElem,dh,dt,tStart1,tEnd1,u0, cV); 

    tStart2 = tPEFRL1(end); tEnd2 = tEnd;
    u02 = uPEFRL1(:, end); 
    [uPEFRL2, EPEFRL2, tPEFRL2] = PEFRL(NElem,dh,dt,tStart2,tEnd2,u02, cV); 

    uPEFRL = [uPEFRL2(:, end-5:end)];
    tPEFRL = [tPEFRL1; tPEFRL2]; EPEFRL = [EPEFRL1; EPEFRL2];

    TimePEFRL(i, 1) = toc; 

    fEvalPEFRL(i, 1) = size(tPEFRL,1)*4*NElem; 
    clear uPEFRL1 EPEFRL1 tPEFRL1 uPEFRL2 EPEFRL2 tPEFRL2; 

%    COMP4
    fprintf(' ... Executing COMP4 ... \n '); 
    uCOMP4 = sparse(2*(NElem+2), tEnd/dt); 

    tic;

    tStart1 = 0; tEnd1 = tEnd/2;  
    [uCOMP41, ECOMP41, tCOMP41] = COMP4(NElem,dh,dt,tStart1,tEnd1,u0, cV); 

    tStart2 = tCOMP41(end); tEnd2 = tEnd;
    u02 = uCOMP41(:, end); 
    [uCOMP42, ECOMP42, tCOMP42] = COMP4(NElem,dh,dt,tStart2,tEnd2,u02, cV); 

    uCOMP4 = [uCOMP42(:, end-5:end)];
    tCOMP4 = [tCOMP41; tCOMP42]; ECOMP4 = [ECOMP41; ECOMP42];

    TimeCOMP4(i, 1) = toc; 
    
    fEvalCOMP4(i, 1) = size(tCOMP4,1)*5*NElem;
    clear uCOMP41 ECOMP41 tCOMP41 uCOMP42 ECOMP42 tCOMP42; 


%     RRK root solver, 4th order
    fprintf(' ... Executing RRK-root solver ... \n '); 
    uRRKr = sparse(2*(NElem+2), tEnd/dt); 

    tic;
    tStart1 = 0; tEnd1 = tEnd/2;
    [uRRK1r, tRRK1r, gamRRK1r, eRRK1r, MaxIt1r] = RRK4Root(NElem, dh, dt, tStart1, tEnd1, u0, cV); 

    tStart2 = tRRK1r(end); tEnd2 = tEnd;
    u02 = uRRK1r(:, end);     
    [uRRK2r, tRRK2r, gamRRK2r, eRRK2r, MaxIt2r] = RRK4Root(NElem, dh, dt, tStart2, tEnd2, u02, cV);
    
    uRRKr = [uRRK2r(:, end-5:end)]; 
    tRRKr = [tRRK1r; tRRK2r];
    eRRKr = [eRRK1r; eRRK2r]; 

    TimeRRKr(i, 1) = toc;
    
    nPerRun = mean([MaxIt1r; MaxIt2r]); 
    fEvalRRKr(i, 1) = size(tRRKr,1)*nPerRun*8*NElem;
    clear uRRK1r tRRK1r gamRRK1r eRRK1r uRRK2r tRRK2r gamRRK2r eRRK2r; 



    % Calculate convergence
    if cV*tEnd < xmax
        uEx4P = exp(-100*(xGrid + cV*tEnd - 0.5).^2);
        uEx4M = exp(-100*(xGrid - cV*tEnd - 0.5).^2);

%         uEx4P = 1/sqrt(sg*2*pi)*exp(-(xGrid + cV*tEnd -muu).^2/2/sg);
%         uEx4M = 1/sqrt(sg*2*pi)*exp(-(xGrid - cV*tEnd -muu).^2/2/sg);

        uEx = 0.5*(uEx4P + uEx4M); 
    else
        uEx = u0(1:NElem+2, 1); 

    end



    [yRK4] = tInterp(tRK4, uRK4, tEnd, NElem);
    normRK4(i, 1) = norm(yRK4(1:NElem+2, 1) - uEx, 'inf'); 

    [yRRK4] = tInterp(tRRK, uRRK, tEnd, NElem);
    normRRK4(i, 1) = norm(yRRK4(1:NElem+2, 1) - uEx, 'inf'); 

    [yFRuth] = tInterp(tFRuth, uFRuth, tEnd, NElem);
    normFRuth(i, 1) = norm(yFRuth(1:NElem+2, 1) - uEx, 'inf'); 

    [yPEFRL] = tInterp(tPEFRL, uPEFRL, tEnd, NElem);
    normPEFRL(i, 1) = norm(yPEFRL(1:NElem+2, 1) - uEx, 'inf'); 

    [yCOMP4] = tInterp(tCOMP4, uCOMP4, tEnd, NElem);
    normCOMP4(i, 1) = norm(yCOMP4(1:NElem+2, 1) - uEx, 'inf'); 

    [yRRK4r] = tInterp(tRRKr, uRRKr, tEnd, NElem);
    normRRK4r(i, 1) = norm(yRRK4r(1:NElem+2, 1) - uEx, 'inf'); 


end


% Calculate order of convergence
for i = 1:size(NE, 1) - 1
    CvRK4(i, 1) = 1/log(2)*log(normRK4(i,1)/normRK4(i+1,1)); 
    CvRRK4(i, 1) = 1/log(2)*log(normRRK4(i,1)/normRRK4(i+1,1)); 
    CvFRuth(i, 1) = 1/log(2)*log(normFRuth(i,1)/normFRuth(i+1,1)); 
    CvPEFRL(i, 1) = 1/log(2)*log(normPEFRL(i,1)/normPEFRL(i+1,1)); 
    CvCOMP4(i, 1) = 1/log(2)*log(normCOMP4(i,1)/normCOMP4(i+1,1)); 
    CvRRK4r(i, 1) = 1/log(2)*log(normRRK4r(i,1)/normRRK4r(i+1,1)); 

end

CvOut = [CvRK4, CvRRK4, CvFRuth, CvPEFRL, CvCOMP4, CvRRK4r]; %
normOut = [normRK4, normRRK4, normFRuth, normPEFRL, normCOMP4, normRRK4r];  %
TimeOut = [TimeRK4, TimeRRK, TimeFRuth, TimePEFRL, TimeCOMP4, TimeRRKr]; %

% Plots
plot(xGrid, uRRK(1:NElem+2, end));
hold on
plot(xGrid, uRK4(1:NElem+2, end));
hold on
plot(xGrid, uFRuth(1:NElem+2, end));
hold on
plot(xGrid, uPEFRL(1:NElem+2, end));
hold on
plot(xGrid, uCOMP4(1:NElem+2, end), '-.', 'LineWidth', 2);
hold on
plot(xGrid, uRRKr(1:NElem+2, end));
hold on
plot(xGrid, uEx, '-k', 'LineWidth', 1); % u0(1:NElem+2, 1)
title('u vs x'); legend('RRK', 'RK4', 'FRuth', 'PEFRL', 'COMP4', 'RRK-root', 'Exact') %
xlabel('x')
ylabel('u(t)')
set(gca,'FontSize',10)
% 

pn1 = 1; pn2 = 120; 
figure
plot(tRRK(pn1:pn2:end), abs( (eRRK(pn1:pn2:end) - eRRK(1,1))/eRRK(1,1)), '-ok', 'MarkerSize', 3)
hold on
plot(tRK4(pn1:pn2:end), abs( (eRK4(pn1:pn2:end) - eRK4(1,1))/eRK4(1,1)), '-*k', 'MarkerSize', 3)
hold on
plot(tFRuth(pn1:pn2:end), abs( (EFRuth(pn1:pn2:end) - EFRuth(1,1))/EFRuth(1,1) ), '-^k', 'MarkerSize', 3)
hold on
plot(tPEFRL(pn1:pn2:end), abs( (EPEFRL(pn1:pn2:end) - EPEFRL(1,1))/EPEFRL(1,1)  ), '-.sk', 'MarkerSize', 3)
hold on
plot(tCOMP4(pn1:pn2:end), abs( (ECOMP4(pn1:pn2:end) - ECOMP4(1,1))/ECOMP4(1,1) ), '--xk', 'MarkerSize', 3)
hold on
plot(tRRKr(pn1:pn2:end), abs( (eRRKr(pn1:pn2:end) - eRRKr(1,1))/eRRKr(1,1) ), '-ob', 'MarkerSize', 2)

set(gca, 'YScale', 'log');  
legend('RRK', 'RK4', 'FRuth', 'PEFRL', 'COMP4', 'RRK-root', 'location', 'best') %
title('Energy En-E0 vs time')
xlabel('time(s)')
ylabel('|(E_n -  E_0)/E_0|')
set(gca,'FontSize',10); grid on;
% % 


figure
plot(NE, TimeRRK, '-ok'); 
hold on; plot(NE, TimeRK4, '-*k');
hold on; plot(NE, TimeFRuth, '-^k');
hold on; plot(NE, TimePEFRL, '-.sk');
hold on; plot(NE, TimeCOMP4, '--xk');
hold on; plot(NE, TimeRRKr, '-ob'); 

legend('RRK', 'RK4', 'FRuth', 'PEFRL', 'COMP4', 'RRK-root', 'location', 'best') %
title('Computation Time Comparison')
xlabel('Mesh Spatial Size, N'); ylabel('Computational Time [s]')
set(gca,'FontSize',10); grid on; 

% Function Evals

figure
loglog(fEvalRRK, normRRK4, '-ok'); 
hold on; loglog(fEvalRK4, normRK4, '-*k');
hold on; loglog(fEvalFRuth, normFRuth, '-^k');
hold on; loglog(fEvalPEFRL, normPEFRL, '-.sk');
hold on; loglog(fEvalCOMP4, normCOMP4, '--xk');
hold on; loglog(fEvalRRKr, normRRK4r, '-ob'); 
hold on; loglog([1e7, 4e7], [1e-2, 1/16*1e-2], '-.r', 'LineWidth', 1.5); 
legend('RRK', 'RK4', 'FRuth', 'PEFRL', 'COMP4', 'RRK-root', ...
    'order = 4', 'location', 'best') %
title('Computation Work Comparison')
xlabel('# of Functional Evaluations'); ylabel('|| U - u_{exact} ||_{\infty}')
set(gca,'FontSize',10); grid on; 





function [uRRK, tRRK, gamRRK, eRRK] = RRK(RKFlag,NElem,dh,dt,tStart,tEnd,u0, cV)
    
%     uRRK = zeros(2*(NElem+2), (tEnd)/dt); 


% Relaxation RK4
     
    Z1 = zeros(NElem+2,NElem+2);
    D = div(4, NElem, dh);
    G = grad(4, NElem, dh);
    L = cV^2*D*G; 
    
    AMat = [Z1, speye(NElem+2,NElem+2);
            L, Z1]; 
    
    AMat = sparse(AMat); 
               
    C(1,1) = 0;
    C(2,1) = 1/2;
    C(3,1) = 1/2;
    C(4,1) = 1;

    A(2,1) = 1/2;
    A(3,1) = 0;   A(3,2) = 1/2;
    A(4,1) = 0;   A(4,2) = 0;     A(4,3) = 1;     

    B(1,1) = 1/6; 
    B(2,1) = 1/3;
    B(3,1) = 1/3;
    B(4,1) = 1/6;

     
    uRRK(:,1) = u0; % initial value
    y = u0; 
    tRRK(1,1) = tStart; % initial time    
    gamRRK(1,1) = 1;
    ut = cV*G*y(1:NElem+2,1);     
    E0 = ut'*ut + ...
            y(NElem+3:end,1)'*y(NElem+3:end,1); 
    eRRK(1,1) = E0; 
    
    iCount = 2; 
    t = tStart + dt; 
   
    while t <= tEnd+dt
        z1 = y; 
        [k1] = AMat*z1; %
        
        % Step 2        
        z2 =  y + dt*A(2,1)*k1;       
        [k2] = AMat*z2;

        % Step 3    
        z3 = y + dt*(A(3,1)*k1 + A(3,2)*k2); 
        [k3] = AMat*z3 ;
        
        % Step 4    
        z4 = y + dt*(A(4,1)*k1 + A(4,2)*k2 + A(4,3)*k3);         
        [k4] = AMat*z4 ; 
        
        BKsum = B(1,1)*k1 + B(2,1)*k2 + B(3,1)*k3 + ...
            B(4,1)*k4; % + B(5,1)*k5;   
    
        switch RKFlag
            case 1
                U = cV*y(1:NElem+2,1); V = y(NElem+3:end,1);
                dU = cV*BKsum(1:NElem+2,1); dV = BKsum(NElem+3:end,1);
                
                AA = U'*G'*G*U + V'*V;
                BB = 2*(U'*G'*G*dU + V'*dV);
                CC = dU'*G'*G*dU + dV'*dV;
                
                gam = (E0 - AA - BB)/(dt* CC);

            case 0
                gam = 1;
        end
        uNew = y + gam*dt*BKsum;  
        uRRK(:,iCount) = uNew;   
        y = uRRK(:,iCount);
        tRRK(iCount,1) = t;
        gamRRK(iCount,1) = gam; 
        
        ut1 = cV*G*uNew(1:NElem+2,1); 
       eRRK(iCount,1) = ut1'*ut1 + ...
            uNew(NElem+3:end,1)'*uNew(NElem+3:end,1) ;   %- E0
        
        iCount = iCount+1; 
        t = t + gam*dt;           
        
    end   

end


function [uFRuth, EFRuth, tFRuth] = FRuth(NElem,dh,dt,tStart,tEnd,u0, cV)

% Forest Ruth
% uFRuth = zeros(2*(NElem+2), (tEnd - tStart)/dt); 

U = u0(1:NElem+2,1); V = u0(NElem+3:end,1); 
uFRuth(:,1) = [U; V]; % initial value

D = div(4, NElem, dh);
G = grad(4, NElem, dh);
L = cV^2*D*G; 

rt = -(1-2^(1/3)-2^(-1/3))/6; % root from Ruth's paper
c1 = rt + 1/2; c2 = -rt; 
c3 = c2; c4 = rt + 1/2;

d1 = 2*rt + 1; d2 = -4*rt - 1;
d3 = d1; d4 = 0; 

ut = cV*G*U;
E0 = ut'*ut + V'*V;
EFRuth(1,1) = E0; 
tFRuth(1,1) = tStart; % initial time    
iCount = 2; 
t = tStart + dt; 

    while t <= tEnd+dt 
        p1 = U + c1*dt*V; 
        q1 = V + d1*dt*L*p1;         
        
        % Step 2        
        p2 = p1 + c2*dt*q1; 
        q2 = q1 + d2*dt*L*p2; 
    
    %
        % Step 3    
        p3 = p2 + c3*dt*q2; 
        q3 = q2 + d3*dt*L*p3; 
    
    
        % Step 4    
        p4 = p3 + c4*dt*q3; 
        q4 = q3; % + d4*dt*L*p4; 
                
        uNew = [p4; q4];                  
        
        ut1 = cV*G*p4; 
        EFRuth(iCount,1) = ut1'*ut1 + q4'*q4 ; %- E0
        uFRuth(:,iCount) = uNew;   
        U = uFRuth(1:NElem+2, iCount);
        V = uFRuth(NElem+3:end, iCount);

        tFRuth(iCount,1) = t;
        
        iCount = iCount+1; 
        t = t + dt;              
        
    end

end


function [uPEFRL, EPEFRL, tPEFRL] = PEFRL(NElem,dh,dt,tStart,tEnd,u0, cV)

% Forest Ruth, Position Extended, Omelyan et al
% uPEFRL = zeros(2*(NElem+2), tEnd/dt); 

U = u0(1:NElem+2,1); V = u0(NElem+3:end,1); 
uPEFRL(:,1) = [U; V]; % initial value

D = div(4, NElem, dh);
G = grad(4, NElem, dh);
L = cV^2*D*G; 

Zeta = 0.1786178958448091E+00;
Lmb = -0.2123418310626054E+00;
Xi = -0.6626458266981849E-1; 

ut = cV*G*U;
E0 = ut'*ut + V'*V;
EPEFRL(1,1) = E0; 
tPEFRL(1,1) = tStart; % initial time    
iCount = 2; 

Zdt = Zeta*dt*V; 

t = tStart + dt; 

    while t <= tEnd+dt 
        p1 = U + Zdt; % Zeta*dt*V; 
        q1 = V + (1-2*Lmb)/2*dt*L*p1; 
        
        % Step 2        
        p2 = p1 + Xi*dt*q1; 
        q2 = q1 + Lmb*dt*L*p2; 
    
    %
        % Step 3    
        p3 = p2 + (1-2*(Xi + Zeta))*dt*q2; 
        q3 = q2 + Lmb*dt*L*p3; 
    
    
        % Step 4    
        p4 = p3 + Xi*dt*q3; 
        q4 = q3 + (1-2*Lmb)/2*dt*L*p4; 
        
        Zdt = Zeta*dt*q4;
        p4 = p4 + Zdt; % Zeta*dt*q4; 
        
        uNew = [p4; q4];                  
    
        ut1 = cV*G*p4;     
        EPEFRL(iCount,1) = ut1'*ut1 + q4'*q4; % - E0; 
        uPEFRL(:,iCount) = uNew;   
        U = uPEFRL(1:NElem+2, iCount);
        V = uPEFRL(NElem+3:end, iCount);
        tPEFRL(iCount,1) = t;
        
        iCount = iCount+1; 
        t = t + dt;              
        
    end

end


function [uCOMP4, ECOMP4, tCOMP4] = COMP4(NElem,dh,dt,tStart,tEnd,u0, cV)

% uCOMP4 = zeros(2*(NElem+2), tEnd/dt); 

U = u0(1:NElem+2,1); V = u0(NElem+3:end,1); 
uCOMP4(:,1) = [U; V]; % initial value

D = div(4, NElem, dh);
G = grad(4, NElem, dh);
L = cV^2*D*G; 

alp1 = (146 + 5*sqrt(19))/540;  bet5 = alp1;
alp2 = (-2 + 10*sqrt(19))/135;  bet4 = alp2; 
alp3 = 1/5;                     bet3 = alp3;
alp4 = (-23 - 20*sqrt(19))/270; bet2 = alp4;
alp5 = (14 - sqrt(19))/108;     bet1 = alp5; 

ut = cV*G*U;
E0 = ut'*ut + V'*V;
ECOMP4(1,1) = E0; 
tCOMP4(1,1) = tStart; % initial time    
iCount = 2; 
t = tStart + dt; 

    while t <= tEnd+dt 
        p1 = U + bet1*dt*V; 
        q1 = V + (bet1 + alp1)*dt*L*p1;         
        
        % Step 2        
        p2 = p1 + (bet2 + alp1)*dt*q1; 
        q2 = q1 + (bet2 + alp2)*dt*L*p2; 
    
    %
        % Step 3    
        p3 = p2 + (bet3 + alp2)*dt*q2; 
        q3 = q2 + (bet3 + alp3)*dt*L*p3; 
    
    
        % Step 4    
        p4 = p3 + (bet4 + alp3)*dt*q3; 
        q4 = q3 + (bet4 + alp4)*dt*L*p4; 

        % Step 5    
        p5 = p4 + (bet5 + alp4)*dt*q4; 
        q5 = q4 + (bet5 + alp5)*dt*L*p5;
        
        p5 = p5 + alp5*dt*q5; 

        uNew = [p5; q5];                  
        
        ut1 = cV*G*p5; 
        ECOMP4(iCount,1) = ut1'*ut1 + q5'*q5; % - E0; 
        uCOMP4(:,iCount) = uNew;   
        U = uCOMP4(1:NElem+2, iCount);
        V = uCOMP4(NElem+3:end, iCount);

        tCOMP4(iCount,1) = t;
        
        iCount = iCount+1; 
        t = t + dt;              
        
    end

end



function [uRRK, tRRK, gamRRK, eRRK, MaxItCount] = RRK4Root(NElem,dh,dt,tStart,tEnd,u0, cV)

% Relaxation RK, using Secant root finder for gamma

Z1 = zeros(NElem+2,NElem+2);
D = div(4, NElem, dh);
G = grad(4, NElem, dh);
robin = robinBC(4, NElem, dh, 1, 0);
L = cV^2*D*G + robin; 

AMat = [Z1, speye(NElem+2,NElem+2);
        L, Z1]; 

AMat = sparse(AMat); 

C(1,1) = 0;   C(2,1) = 1/2;
C(3,1) = 1/2; C(4,1) = 1;

A(2,1) = 1/2;
A(3,1) = 0;   A(3,2) = 1/2;
A(4,1) = 0;   A(4,2) = 0;     A(4,3) = 1;    

B(1,1) = 1/6; B(2,1) = 1/3;
B(3,1) = 1/3; B(4,1) = 1/6;


uRRK(:,1) = u0; % initial value
y = u0; 
tRRK(1,1) = tStart; % initial time    
gamRRK(1,1) = 1;
ut = cV*G*y(1:NElem+2,1); 
E0 = ut'*ut + ...
        y(NElem+3:end,1)'*y(NElem+3:end,1); 
eRRK(1,1) = E0; 

iCount = 2; 
t = tStart + dt; MaxIt = 30; MinErr = 1E-12; 

    while t <= tEnd+dt    
        x0 = 0.99999; x1 = 1.000001; 

        % Secant root finder for gamma
        for i = 1:MaxIt
            u_x0 = RKSolver(y, AMat, dt, x0);
            u_x1 = RKSolver(y, AMat, dt, x1);
            
            ut0 = cV*G*u_x0(1:NElem+2,1); 
            ut1 = cV*G*u_x1(1:NElem+2,1); 
            
            E_x0 = abs(ut0'*ut0 + u_x0(NElem+3:end,1)'*u_x0(NElem+3:end,1) - E0);
            E_x1 = abs(ut1'*ut1 + u_x1(NElem+3:end,1)'*u_x1(NElem+3:end,1) - E0);

            x2 = x1 - E_x1*(x1 - x0)/(E_x1 - E_x0);    
            x0 = x1; x1 = x2;         
            if abs(x1 - x0) <= MinErr, break, end
    
        end

        MaxItCount(iCount, 1) = i; 
        gam = x0; uNew = u_x0; 

        uRRK(:,iCount) = uNew;   
        y = uRRK(:,iCount);
        tRRK(iCount,1) = t;
        gamRRK(iCount,1) = gam; 
        
        ut2 = cV*G*uNew(1:NElem+2,1); 
        eRRK(iCount,1) = ut2'*ut2 + ...
            uNew(NElem+3:end,1)'*uNew(NElem+3:end,1);  
        
        iCount = iCount+1; 
        t = t + gam*dt;           
        
        
    end

end


function uNew = RKSolver(y, AMat, dt, gam)

C(1,1) = 0;   C(2,1) = 1/2;
C(3,1) = 1/2; C(4,1) = 1;

A(2,1) = 1/2;
A(3,1) = 0;   A(3,2) = 1/2;
A(4,1) = 0;   A(4,2) = 0;     A(4,3) = 1;    

B(1,1) = 1/6; B(2,1) = 1/3;
B(3,1) = 1/3; B(4,1) = 1/6;

    z1 = y;        
    [k1] = AMat*z1;    
    
    % Step 2        
    z2 =  y + dt*A(2,1)*k1;        
    [k2] = AMat*z2;   

    % Step 3    
    z3 = y + dt*(A(3,1)*k1 + A(3,2)*k2);             
    [k3] = AMat*z3;

    % Step 4    
    z4 = y + dt*(A(4,1)*k1 + A(4,2)*k2 + A(4,3)*k3);         
    [k4] = AMat*z4;
    
    uNew = y  + gam*dt*(B(1,1)*k1 + B(2,1)*k2 + B(3,1)*k3 + ...
        B(4,1)*k4);          

end



function [yR] = tInterp(tR, uR, tEnd, NElem)
% interpolate for t = tEnd, used for convergence

    t1 = tR(end-1, 1); t2 = tR(end, 1);
    u1 = uR(1:NElem+2, end-1); u2 = uR(1:NElem+2, end); 
    yR = u1 + (u2 - u1)*(tEnd - t1)/(t2 - t1); 

end

