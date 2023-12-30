clc; clear; close all;
% Symplectic Integrator, Example
% ASrinivasan, 23Jan2023
% 1D Shallow Water wave eqn
% https://arxiv.org/pdf/2009.09641.pdf
addpath('C:\Anand\Acer_Data\SDSU\MOLE\mole-master\mole_MATLAB')

NE = 100*[1, 2, 4, 8, 16]';  

for i = 1:size(NE, 1)
    fprintf('Executing iteration %i of %i ... \n', i, size(NE, 1)); 

    NElem = NE(i, 1);
    xmin = -30; xmax = 30;  
    dh = (xmax-xmin)/NElem; 
    tEnd = 15;   % 200
    CFL = 0.1; cV = 1; 
    dt =  CFL*dh/cV; 
    
    % init val
    xGrid = [xmin xmin+dh/2:dh:xmax-dh/2 xmax]';
    xNod = [xmin:dh:xmax]';
%     u0 = 15/4*(cosh(3*sqrt(2/5)*(xGrid)) - 2).*sech(3/sqrt(10)*(xGrid)).^4; 
%     v0 = 15/2*sech(3/sqrt(10)*(xGrid)).^2; 

    u0 = 1 + 0.5/5*exp(-1*xGrid.^2);    
    v0 = zeros(size(xGrid,1),1);

%     muu = 0.; sg = 0.25;
%     u0 = 1 + 1/sqrt(sg*2*pi)*exp(-(xGrid-muu).^2/2/sg);
%     v0 = zeros(size(xGrid,1),1);

    u0 = [u0;v0]; 
   

%     RRK, 4th order
    fprintf(' ... Executing RRK ... \n '); 
    uRRK = sparse(2*(NElem+2), floor(tEnd/dt) ); 

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
    clear uRRK1 tRRK1 gamRRK1 eRRK1 uRRK2 tRRK2 gamRRK2 eRRK2; %

%     RK4, 4th order
    fprintf(' ... Executing RK4 ... \n '); 
%     uRK4 = sparse(2*(NElem+2), floor(tEnd/dt) ); 

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
    uFRuth = sparse(2*(NElem+2), floor(tEnd/dt)); 
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
    uPEFRL = sparse(2*(NElem+2), floor(tEnd/dt)); 

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
    uCOMP4 = sparse(2*(NElem+2), floor(tEnd/dt)); 

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
    uRRKr = sparse(2*(NElem+2), floor(tEnd/dt)); 

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
    uEx = 15/4*(cosh(3*sqrt(2/5)*(xGrid - 5/2*tEnd)) - 2).* ...
            sech(3/sqrt(10)*(xGrid - 5/2*tEnd)).^4; 
    vEx = 15/2*sech(3/sqrt(10)*(xGrid - 5/2*tEnd)).^2; 


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
figure
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
% hold on
% plot(xGrid, uEx, '-k', 'LineWidth', 1); % u0(1:NElem+2, 1)
title('u vs x'); legend('RRK', 'RK4', 'FRuth', 'PEFRL', 'COMP4', 'RRK-root') %
xlabel('x')
ylabel('u(x,t)')
set(gca,'FontSize',10)
% 

figure
plot(xGrid, uRRK(NElem+3:end, end));
hold on
plot(xGrid, uRK4(NElem+3:end, end));
hold on
plot(xGrid, uFRuth(NElem+3:end, end));
hold on
plot(xGrid, uPEFRL(NElem+3:end, end));
hold on
plot(xGrid, uCOMP4(NElem+3:end, end), '-.', 'LineWidth', 2);
hold on
plot(xGrid, uRRKr(NElem+3:end, end));
% hold on
% plot(xGrid, uEx, '-k', 'LineWidth', 1); % u0(1:NElem+2, 1)
title('v vs x'); legend('RRK', 'RK4', 'FRuth', 'PEFRL', 'COMP4', 'RRK-root' ) %
xlabel('x')
ylabel('v(x,t)')
set(gca,'FontSize',10)


pn1 = 1; pn2 = 50; 
figure
plot(tRRK(pn1:pn2:end), abs( (eRRK(pn1:pn2:end) - eRRK(1,1))/eRRK(1,1)), '-ok', 'MarkerSize', 3)
hold on
plot(tRK4(pn1:pn2:end), abs( (eRK4(pn1:pn2:end) - eRK4(1,1))/eRK4(1,1)), '-*k', 'MarkerSize', 3)
hold on
plot(tFRuth(pn1:pn2:end), abs( (EFRuth(pn1:pn2:end) - EFRuth(1,1))/EFRuth(1,1) ), '-^k', 'MarkerSize', 6)
hold on
plot(tPEFRL(pn1:pn2:end), abs( (EPEFRL(pn1:pn2:end) - EPEFRL(1,1))/EPEFRL(1,1)  ), '-.sr', 'MarkerSize', 3)
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
title('Computation Time Comparison, Shallow Water Equations')
xlabel('Mesh Spatial Size, N'); ylabel('Computational Time [s]')
set(gca,'FontSize',10); grid on; 
% % 
% % % Function Evals
% % 
% % figure
% % loglog(fEvalRRK, normRRK4, '-ok'); 
% % hold on; loglog(fEvalRK4, normRK4, '-*k');
% % hold on; loglog(fEvalFRuth, normFRuth, '-^k');
% % hold on; loglog(fEvalPEFRL, normPEFRL, '-.sk');
% % hold on; loglog(fEvalCOMP4, normCOMP4, '--xk');
% % hold on; loglog(fEvalRRKr, normRRK4r, '-ob'); 
% % hold on; loglog([1e7, 4e7], [1e-2, 1/16*1e-2], '-.r', 'LineWidth', 1.5); 
% % legend('RRK', 'RK4', 'FRuth', 'PEFRL', 'COMP4', 'RRK-root', ...
% %     'order = 4', 'location', 'best') %
% % title('Computation Work Comparison')
% % xlabel('# of Functional Evaluations'); ylabel('|| U - u_{exact} ||_{\infty}')
% % set(gca,'FontSize',10); grid on; 
% % 




function [uRRK, tRRK, gamRRK, eRRK] = RRK(RKFlag,NElem,dh,dt,tStart,tEnd,u0, cV)
    
%     uRRK = zeros(2*(NElem+2), (tEnd)/dt); 


% Relaxation RK4
     
    D = div(4, NElem, dh);
    IntD = interpC2N(NElem, 4); 

    DIntD = D*IntD; 

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
    EE = y(1:NElem+2); %eta
    UU = y(NElem+3:end); %u
    
    E0 = EE'*EE + UU'*UU + EE'*UU.^2; 
    eRRK(1,1) = E0; 
    
    iCount = 2; 
    t = tStart + dt; 
   
    while t <= tEnd+dt
        z1 = y; 
        [k1] = - [DIntD*z1(NElem+3:end,1) + DIntD*(z1(1:NElem+2,1).*z1(NElem+3:end,1));
                  DIntD*z1(1:NElem+2) +  z1(NElem+3:end).*DIntD*z1(NElem+3:end) ]; %
               

        % Step 2        
        z2 =  y + dt*A(2,1)*k1;       
        [k2] = - [DIntD*z2(NElem+3:end,1) + DIntD*(z2(1:NElem+2,1).*z2(NElem+3:end,1));
                  DIntD*z2(1:NElem+2) +  z2(NElem+3:end).*DIntD*z2(NElem+3:end) ]; %

        % Step 3    
        z3 = y + dt*(A(3,1)*k1 + A(3,2)*k2); 
        [k3] = - [DIntD*z3(NElem+3:end,1) + DIntD*(z3(1:NElem+2,1).*z3(NElem+3:end,1));
                  DIntD*z3(1:NElem+2) +  z3(NElem+3:end).*DIntD*z3(NElem+3:end) ]; %
        
        % Step 4    
        z4 = y + dt*(A(4,1)*k1 + A(4,2)*k2 + A(4,3)*k3);         
        [k4] = - [DIntD*z4(NElem+3:end,1) + DIntD*(z4(1:NElem+2,1).*z4(NElem+3:end,1));
                  DIntD*z4(1:NElem+2) +  z4(NElem+3:end).*DIntD*z4(NElem+3:end) ]; %
        
        BKsum = B(1,1)*k1 + B(2,1)*k2 + B(3,1)*k3 + ...
            B(4,1)*k4; % + B(5,1)*k5;   
    
        switch RKFlag
            case 1
                EE = y(1:NElem+2,1); UU = y(NElem+3:end,1);
                dEE = BKsum(1:NElem+2,1); dUU = BKsum(NElem+3:end,1);
                
                BB = 2*EE'*dEE + 2*UU'*dUU + 2*UU'*(EE.*dUU) + dEE'*(UU.*UU);                
                CC = dEE'*dEE + dUU'*dUU + EE'*(dUU.*dUU) + 2*dEE'*(UU.*dUU);
                DD = dEE'*(dUU.*dUU);

                gam = (-CC + (CC^2 - 4*DD*BB)^(1/2))/(2*DD*dt);

            case 0
                gam = 1;
        end
        uNew = y + gam*dt*BKsum;  
        uRRK(:,iCount) = uNew;   
        y = uRRK(:,iCount);
        tRRK(iCount,1) = t;
        gamRRK(iCount,1) = gam; 
        
        eRRK(iCount,1) = uNew(1:NElem+2,1)'*uNew(1:NElem+2,1) + ...
            uNew(NElem+3:end,1)'*uNew(NElem+3:end,1) + ...
            uNew(1:NElem+2,1)'*(uNew(NElem+3:end,1).*uNew(NElem+3:end,1));    
        
        iCount = iCount+1; 
        t = t + gam*dt;           
        
    end   

end


function [uFRuth, EFRuth, tFRuth] = FRuth(NElem,dh,dt,tStart,tEnd,u0, cV)

% Forest Ruth
% uFRuth = zeros(2*(NElem+2), (tEnd - tStart)/dt); 

U = u0(1:NElem+2); V = u0(NElem+3:end); 
uFRuth(:,1) = [U; V]; % initial value

D = div(4, NElem, dh);
IntD = interpC2N(NElem, 4); 
DIntD = D*IntD; 

rt = -(1-2^(1/3)-2^(-1/3))/6; % root from Ruth's paper
c1 = rt + 1/2; c2 = -rt; 
c3 = c2; c4 = rt + 1/2;

d1 = 2*rt + 1; d2 = -4*rt - 1;
d3 = d1; d4 = 0; 

EE = U; %eta
UU = V; %u

E0 = EE'*EE + UU'*UU + (EE'*UU).^2; 
EFRuth(1,1) = E0; 
tFRuth(1,1) = tStart; % initial time    
iCount = 2; 
t = tStart + dt; 

    while t <= tEnd+dt 
        p1 = U - c1*dt*( DIntD*((1 + U).*V) ); 
        q1 = V - d1*dt*( DIntD*p1 + V.*DIntD*V );         %
        
        % Step 2        
        p2 = p1 - c2*dt*( DIntD*((1 + p1).*q1) );
        q2 = q1 - d2*dt*( DIntD*p2 + q1.*DIntD*q1 ); %
    
    %
        % Step 3    
        p3 = p2 - c3*dt*( DIntD*((1 + p2).*q2) ); 
        q3 = q2 - d3*dt*( DIntD*p3 + q2.*DIntD*q2 );    % 
    
    
        % Step 4    
        p4 = p3 - c4*dt*( DIntD*((1 + p3).*q3) ); 
        q4 = q3; 
                
        uNew = [p4; q4]; 

        
        EFRuth(iCount,1) = p4'*p4 + q4'*q4 + (p4'*q4).^2;  
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

U = u0(1:NElem+2,1); V = u0(NElem+3:end,1); 
uPEFRL(:,1) = [U; V]; % initial value

D = div(4, NElem, dh);
IntD = interpC2N(NElem, 4);
DIntD = D*IntD; 

Zeta = 0.1786178958448091E+00;
Lmb = -0.2123418310626054E+00;
Xi = -0.6626458266981849E-1; 

E0 = U'*U + V'*V + (U'*V).^2;
EPEFRL(1,1) = E0; 
tPEFRL(1,1) = tStart; % initial time    
iCount = 2; 

t = tStart + dt; 

    while t <= tEnd+dt 
        p1 = U - Zeta*dt*( DIntD*( (1+U).*V ) ); 
        q1 = V - (1-2*Lmb)/2*dt*( DIntD*p1 + V.*DIntD*V ); 
        
        % Step 2        
        p2 = p1 - Xi*dt*( DIntD*( (1+p1).*q1 ) ); 
        q2 = q1 - Lmb*dt*( DIntD*p2 + q1.*DIntD*q1 ); 
    
    %
        % Step 3    
        p3 = p2 - (1-2*(Xi + Zeta))*dt*( DIntD*( (1+p2).*q2 ) ); 
        q3 = q2 - Lmb*dt*( DIntD*p3 + q2.*DIntD*q2 );  
    
    
        % Step 4    
        p4 = p3 - Xi*dt*( DIntD*( (1+p3).*q3 ) ); 
        q4 = q3 - (1-2*Lmb)/2*dt*( DIntD*p4 + q3.*DIntD*q3 );  
        
        p4 = p4 - Zeta*dt*( DIntD*( (1+p4).*q4 ) );  
        
        uNew = [p4; q4];                  
    
        EPEFRL(iCount,1) = p4'*p4 + q4'*q4 + (p4'*q4).^2; 
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
IntD = interpC2N(NElem, 4);
DIntD = D*IntD; 

alp1 = (146 + 5*sqrt(19))/540;  bet5 = alp1;
alp2 = (-2 + 10*sqrt(19))/135;  bet4 = alp2; 
alp3 = 1/5;                     bet3 = alp3;
alp4 = (-23 - 20*sqrt(19))/270; bet2 = alp4;
alp5 = (14 - sqrt(19))/108;     bet1 = alp5; 

E0 = U'*U + V'*V + (U'*V).^2;
ECOMP4(1,1) = E0; 
tCOMP4(1,1) = tStart; % initial time    
iCount = 2; 
t = tStart + dt; 

    while t <= tEnd+dt 
        p1 = U - bet1*dt*( DIntD*((1+U).*V) ); 
        q1 = V - (bet1 + alp1)*dt*( DIntD*p1 + V.*DIntD*V );         
        
        % Step 2        
        p2 = p1 - (bet2 + alp1)*dt*( DIntD*((1+p1).*q1) ); 
        q2 = q1 - (bet2 + alp2)*dt*( DIntD*p2 + q1.*DIntD*q1 );    
    
    %
        % Step 3    
        p3 = p2 - (bet3 + alp2)*dt*( DIntD*((1+p2).*q2) );  
        q3 = q2 - (bet3 + alp3)*dt*( DIntD*p3 + q2.*DIntD*q2 );
    
    
        % Step 4    
        p4 = p3 - (bet4 + alp3)*dt*( DIntD*((1+p3).*q3) );   
        q4 = q3 - (bet4 + alp4)*dt*( DIntD*p4 + q3.*DIntD*q3 ); 

        % Step 5    
        p5 = p4 - (bet5 + alp4)*dt*( DIntD*((1+p4).*q4) ); 
        q5 = q4 - (bet5 + alp5)*dt*( DIntD*p5 + q4.*DIntD*q4 ); 
        
        p5 = p5 - alp5*dt*( DIntD*((1+p5).*q5) );  

        uNew = [p5; q5];                  
        
        ECOMP4(iCount,1) = p5'*p5 + q5'*q5 + (p5'*q5).^2;   
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

D = div(4, NElem, dh);
IntD = interpC2N(NElem, 4);
DIntD = D*IntD; 


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
E0 = y(1:NElem+2,1)'*y(1:NElem+2,1) + ...
        y(NElem+3:end,1)'*y(NElem+3:end,1) + ...
        y(1:NElem+2,1)'*y(NElem+3:end,1).^2; 
eRRK(1,1) = E0; 

iCount = 2; 
t = tStart + dt; MaxIt = 50; MinErr = 1E-12; 

    while t <= tEnd+dt    
        x0 = 0.9999999; x1 = 1.000001;                 

        % Secant root finder for gamma
        for i = 1:MaxIt
            u_x0 = RKSolver(y, DIntD, dt, x0, NElem);
            u_x1 = RKSolver(y, DIntD, dt, x1, NElem);
            
            ut0 = u_x0(1:NElem+2,1); vt0 = u_x0(NElem+3:end,1); 
            ut1 = u_x1(1:NElem+2,1); vt1 = u_x1(NElem+3:end,1);
                        
            E_x0 = abs(ut0'*ut0 + vt0'*vt0 + ut0'*vt0.^2 - E0);
            E_x1 = abs(ut1'*ut1 + vt1'*vt1 + ut1'*vt1.^2 - E0); 

            x2 = x1 - E_x1*(x1 - x0)/(E_x1 - E_x0); 
            x0 = x1; x1 = x2;         
            if abs(E_x1 - E_x0) <= MinErr, break, end
    
        end

        MaxItCount(iCount, 1) = i; 
        gam = x0; 
        uNew = u_x0; 

        uRRK(:,iCount) = uNew;   
        y = uRRK(:,iCount);
        tRRK(iCount,1) = t;
        gamRRK(iCount,1) = gam; 
        
        ut2 = uNew(1:NElem+2,1); 
        eRRK(iCount,1) = ut2'*ut2 + ...
            uNew(NElem+3:end,1)'*uNew(NElem+3:end,1) + ...
            ut2'*uNew(NElem+3:end,1).^2;  
        
        iCount = iCount+1; 
        t = t + gam*dt;           
        
        
    end

end


function uNew = RKSolver(y, DIntD, dt, gam, NElem)

C(1,1) = 0;   C(2,1) = 1/2;
C(3,1) = 1/2; C(4,1) = 1;

A(2,1) = 1/2;
A(3,1) = 0;   A(3,2) = 1/2;
A(4,1) = 0;   A(4,2) = 0;     A(4,3) = 1;    

B(1,1) = 1/6; B(2,1) = 1/3;
B(3,1) = 1/3; B(4,1) = 1/6;

    % Step 1
    z1 = y; 
    [k1] = - [DIntD*z1(NElem+3:end,1) + DIntD*(z1(1:NElem+2,1).*z1(NElem+3:end,1));
              DIntD*z1(1:NElem+2) +  z1(NElem+3:end).*DIntD*z1(NElem+3:end) ]; %          

    % Step 2        
    z2 =  y + dt*A(2,1)*k1;       
    [k2] = - [DIntD*z2(NElem+3:end,1) + DIntD*(z2(1:NElem+2,1).*z2(NElem+3:end,1));
              DIntD*z2(1:NElem+2) +  z2(NElem+3:end).*DIntD*z2(NElem+3:end) ]; %

    % Step 3    
    z3 = y + dt*(A(3,1)*k1 + A(3,2)*k2); 
    [k3] = - [DIntD*z3(NElem+3:end,1) + DIntD*(z3(1:NElem+2,1).*z3(NElem+3:end,1));
              DIntD*z3(1:NElem+2) +  z3(NElem+3:end).*DIntD*z3(NElem+3:end) ]; %
    
    % Step 4    
    z4 = y + dt*(A(4,1)*k1 + A(4,2)*k2 + A(4,3)*k3);         
    [k4] = - [DIntD*z4(NElem+3:end,1) + DIntD*(z4(1:NElem+2,1).*z4(NElem+3:end,1));
              DIntD*z4(1:NElem+2) +  z4(NElem+3:end).*DIntD*z4(NElem+3:end) ]; %
    
    uNew = y  + gam*dt*(B(1,1)*k1 + B(2,1)*k2 + B(3,1)*k3 + ...
        B(4,1)*k4);          

end



function [yR] = tInterp(tR, uR, tEnd, NElem)
% interpolate for t = tEnd, used for convergence

    t1 = tR(end-1, 1); t2 = tR(end, 1);
    u1 = uR(1:NElem+2, end-1); u2 = uR(1:NElem+2, end); 
    yR = u1 + (u2 - u1)*(tEnd - t1)/(t2 - t1); 

end

