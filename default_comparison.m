% Bankruptcy choice - Comparing 3 methods
% Code made by Gustavo Mellior
% This code compares LCP, Splitting method (SM) and Nuno,
% Hurtado and Thomas' (2023) method (NHTM)

clear; clc; close all

s           = 2;                      % CRRA utility
rho         = 0.05;                   % Discount rate
psi         = 2;                      % Mental cost of bankruptcy
csubs       = 1e-6;                   % Subsistence level of consumption - helps keep slope of V positive, see below
plotconv    = 0;                      % Set to 1 if you want to see plots during loops - faster if set to zero
z1          = 0.75;                   % Low income flow
z2          = 1.25;                   % High income
z           = [z1,z2];              
la1         = 0.25;                   % Poisson rate low to high
la2         = 0.25;                   % Poisson high to low
la          = [la1,la2];              %
I           = 300;                    % Wealth grid points
amin        = -4;                     % Exogenous debt limit
amax        = 4;                      % Maximum level of wealth
a           = linspace(amin,amax,I)'; % Wealth grid
da          = (amax-amin)/(I-1);      % Wealth step
aa          = [a,a];
zz          = ones(I,1)*z;
r           = 0.035./((0.5 + exp(aa-0.15)).*(aa<0) + ones(I,2).*(aa>0)); % Interest rate
maxit       = 100;                    % Max number of loops in V (when using LCP)
compute_VSM = 0;                      % Set to 1 if you want to compute V with the SM method (it takes a long time)
maxitSM     = 5000000;
maxitNM     = 9000000;                % Max number of loops in V (when using Nuno, Hurtado and Thomas)
crit        = 10^(-6);                % Set very small so that LCP error clears - see further below
Delta       = 1000;                   % Time step
dVf         = zeros(I,2);             % Pre-allocate memory
dVb         = zeros(I,2);
c           = zeros(I,2);
Aswitch     = [-speye(I)*la(1),speye(I)*la(1);speye(I)*la(2),-speye(I)*la(2)];
% Initial guess of non-default V
v0(:,1)     = ((z(1) + r(:,1).*a).^(1-s) - 1)/(1-s)/rho;
v0(:,2)     = ((z(2) + r(:,1).*a).^(1-s) - 1)/(1-s)/rho;
v           = v0;

% Solver options and initial guess for c(amin)
options   = optimset('Display','off','MaxIter', 10000, 'MaxFunEvals', 5000);
x0        = 0.05*z(1,1);
% Bankruptcy value function
Vstar     = ((z(2)*1.15 + 0.025*a.*(a<0)).^(1-s) - 1)/(1-s)/rho;
Vstar     = [Vstar-psi;Vstar-5000];
Vstarplot       = reshape(Vstar,I,2); % Just for plotting purposes
Vstarplot(aa>0) = NaN;                % Just for plotting purposes

%% Default not allowed - this will be used as a starting guess for bankruptcy choice
for n=1:maxit
    V = v;
    % Forward difference
    dVf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVf(I,:)     = (z + r(end,:).*amax).^(-s);
    % Backward difference
    dVb(2:I,:)   = (V(2:I,:)-V(1:I-1,:))/da;
    dVb(1,:)     = (z + r(1,:).*amin).^(-s); % State boundary constraint
    
    dVf = max(dVf,csubs); % Helps with stability to guarantee slopes are non negative
    dVb = max(dVb,csubs);
        
    % Consumption and savings with forward difference
    cf  = dVf.^(-1/s);
    ssf = zz + r.*aa - cf;
    % Consumption and savings with backward difference
    cb  = dVb.^(-1/s);
    ssb = zz + r.*aa - cb;
    % Consumption and derivative of value function at the temporary steady state
    c0  = zz + r.*aa;
    dV0 = c0.^(-s);
    
    % Upwind
    If        = ssf > 0; %positive drift --> forward difference
    Ib        = ssb < 0; %negative drift --> backward difference
    I0        = (1-If-Ib); %at steady state
    c         = cf.*If + cb.*Ib + c0.*I0;
    u         = (c.^(1-s)-1)/(1-s);
    % Build the A matrix
    X         = - min(ssb,0)/da;
    Y         = - max(ssf,0)/da + min(ssb,0)/da;
    Z         =   max(ssf,0)/da;
    A1        = spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
    A2        = spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
    A         = [A1,sparse(I,I);sparse(I,I),A2] + Aswitch;
    B         = (rho + 1/Delta)*speye(2*I) - A;
    u_stacked = [u(:,1);u(:,2)];
    V_stacked = [V(:,1);V(:,2)];
    b         = u_stacked + V_stacked/Delta;
    V_stacked = B\b; 
    V         = [V_stacked(1:I),V_stacked(I+1:2*I)];
    Vchange   = V - v;
    v         = V;
    dist(n)   = max(abs(Vchange(:)));
    if dist(n)<crit
        % disp('No default value function Converged, Iteration = ')
        % disp(n)
        break
    end
end


%% LCP

LCPerror        = 1;                  % See related comments below
fct             = -1;                 % See comments below
vnodef          = v;
% V^C = value function no default - V^D = value function default
for n=1:maxit
    V = v;
    % Find c(amin)
    cT          = zeros(1,1);
    a3store     = cT;
    a2store     = cT;
    x0          = 0.05*z(1,1);            % Re-start this at every nth loop
    driftminusC = r(1,1)*amin+zz(1,1);    % Autarky cons
    psflow      = la(1)*V(1,2);           % Poisson flow in HJB
    dv          = (rho+la(1))*Vstar(1,1); % Multip denominator of V^C out, affects default value net of mental cost
    params      = [s driftminusC psflow 1 dv]; % Pass arguments to solver
    while fct<0                           % Exploit positive slope of VMatching
        x0 = x0 + 0.1;
        myf1         = @(x) cTsolver(x,params);
        [a1, a2, a3] = fsolve(myf1, x0, options);
        fct          = s*a1.^(-s)-s*a1.^(-s-1)*driftminusC; % Derivative of VMatching wrt c
    end
    fct          = -1;                 % Reset fct for the next loop
    cT           = a1;                 % c(amin)
    cT           = min(cT,zz(1,1)*10); % We don't want strange cT, so let's give a bound to potential results
    cT           = max(cT,csubs);      % Might be unnecessary, but still keep it just in case
    if a2>0.0001                       % If a2 is positive -> normal V with no bankruptcy
       cT = r(1,1)*amin +zz(1,1); 
    elseif (a3<0)&&(a2<0) % Rarely used
       cT = 1000;
    end
    
    % Forward difference
    dVf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVf(I,:)     = (z(end,:) + r(end,:).*amax).^(-s);
    % Backward difference
    dVb(2:I,:)   = (V(2:I,:)-V(1:I-1,:))/da;
    dVb(1,1)     = (max(cT,r(1,1)*amin + zz(1,1))).^(-s); % Only low income type can default
    dVb(1,2)     = (z(2) + r(1,2).*amin).^(-s);
    dVf          = max(dVf,csubs);
    dVb          = max(dVb,csubs);

    % Consumption and savings with forward difference
    cf = (dVf).^(-1/s);
    sf = r.*aa+zz-cf;
    % Consumption and savings with backward difference
    cb = (dVb).^(-1/s);
    sb = r.*aa+zz-cb;
    % Consumption and derivative of value function at the temporary steady state
    c0 = r.*aa + zz;
    c0 = max(c0,csubs); 

    Hf = cf.^(1-s)/(1-s) - 1/(1-s) + dVf.*sf;
    Hb = cb.^(1-s)/(1-s) - 1/(1-s) + dVb.*sb;
    H0 = c0.^(1-s)/(1-s) - 1/(1-s);
    
    % More flexible approach allowing for non concavity
    Iunique = (sb<0).*(1-(sf>0)) + (1-(sb<0)).*(sf>0);
    Iboth   = (sb<0).*(sf>0);
    Ib      = Iunique.*(sb<0).*(Hb>H0) + Iboth.*(Hb>Hf).*(Hb>H0);
    If      = Iunique.*(sf>0).*(Hf>H0) + Iboth.*(Hf>Hb).*(Hf>H0);
    I0      = 1-Ib-If;
    
    % Use upwind scheme to get consumption
    c       = cf.*If + cb.*Ib + c0.*I0;
    u       = (c.^(1-s)-1)/(1-s);
    adot    = sf.*If+sb.*Ib;
    % Build the A matrix
    X       = -Ib.*sb./da;
    Y       = -If.*sf./da + Ib.*sb./da;
    Z       = If.*sf./da;
    
    % ******** Fix entries on A matrix ***********
    Y(1)    = Y(1) - min(adot(1,1)/da,0); % Y(1) + X(1)
    u(1,1)  = u(1,1) + (cb(1,1))^(-s)*min(adot(1,1)/da,0);
    % ******** Fix entries on A matrix ***********
    
    % A matrix
    A1            =spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
    A2            =spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
    A             = [A1,sparse(I,I);sparse(I,I),A2] + Aswitch;
    B             = (1/Delta + rho)*speye(2*I) - A;
    
    u_stacked     = reshape(u,2*I,1);
    V_stacked     = reshape(V,2*I,1);
    % LCP
    vec           = u_stacked + V_stacked/Delta;
    qlcp          = -vec + B*Vstar; 
    z0            = V_stacked-Vstar;
    zans          = LCP(B,qlcp,zeros(2*I,1),Inf*ones(2*I,1),z0,0);
    LCP_error     = max(abs(zans.*(B*zans + qlcp)));
    LCPerror(n+1) = LCP_error; % Track error - sometimes it gets stuck and slows down things
    
    V_stacked     = zans+Vstar; % Get new value function
    V             = reshape(V_stacked,I,2);
    Vchange       = V - v;
    v             = V;
    % We need to set a laxer crit for LCPerror and exit threshold
    spindex = find(abs(V_stacked-Vstar)<crit*1000); % Indicator of bankruptcy
    if isempty(spindex)==1
        spindex = 1;
    end
    
    dist(n) = max(abs(Vchange(:)));
    % Convergence
    if (dist(n)<crit)&&(LCP_error<crit*1000) % Laxer to avoid cycles on LCPerror
        conv = 1;
         disp('****** HJBVI SOLVED ******');
        break
    else
        sprintf('Change in V = %0.10g ---- LCP_error = %0.10g',[dist(n) LCP_error])
        if (LCPerror(n+1)-LCPerror(n))==0 % When this equals zero the algo might get stuck
            disp('LCPerror getting stuck');
            break
        end
    end

end
VLCP = V;

%% Splitting method
v       = vnodef;
dt      = 0.0001;
if compute_VSM==1
    for n=1:maxitSM
        V = v;
    
        % Find c(amin)
        cT          = zeros(1,1);
        a3store     = cT;
        a2store     = cT;
        x0          = 0.05*z(1,1);            % Re-start this at every nth loop
        driftminusC = r(1,1)*amin+zz(1,1);    % Autarky cons
        psflow      = la(1)*V(1,2);           % Poisson flow in HJB
        dv          = (rho+la(1))*Vstar(1,1); % Multip denominator of V^C out, affects default value net of mental cost
        params      = [s driftminusC psflow 1 dv]; % Pass arguments to solver
        while fct<0                           % Exploit positive slope of VMatching
            x0 = x0 + 0.1;
            myf1         = @(x) cTsolver(x,params);
            [a1, a2, a3] = fsolve(myf1, x0, options);
            fct          = s*a1.^(-s)-s*a1.^(-s-1)*driftminusC; % Derivative of VMatching wrt c
        end
        fct          = -1; % Reset fct for the next loop
        cT           = a1; % c(amin)
        cT           = min(cT,zz(1,1)*10); % We don't want strange cT, so let's give a bound to potential results
        cT           = max(cT,csubs);      % Might be unnecessary, but still keep it just in case
        if a2>0.0001 % If a2 is positive -> normal V with no bankruptcy
           cT = r(1,1)*amin +zz(1,1); 
        elseif (a3<0)&&(a2<0) % Rarely used
           cT = 1000;
        end
    
        % Forward difference
        dVf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
        dVf(I,:)     = (z(end,:) + r(end,:).*amax).^(-s);
        % Backward difference
        dVb(2:I,:)   = (V(2:I,:)-V(1:I-1,:))/da;
        dVb(1,1)     = (max(cT,r(1,1)*amin + zz(1,1))).^(-s); % Only low income type can default
        dVb(1,2)     = (z(2) + r(1,2).*amin).^(-s);
    
        dVf          = max(dVf,csubs);     
        dVb          = max(dVb,csubs);
        cf           = dVf.^(-1/s);
        ssf          = zz + r.*aa - cf;
        cb           = dVb.^(-1/s);
        ssb          = zz + r.*aa - cb;
        c0           = zz + r.*aa;
        dV0          = c0.^(-s);
        % Hamiltonians
        Hf           = cf.^(1-s)/(1-s) - 1/(1-s) + dVf.*sf;
        Hb           = cb.^(1-s)/(1-s) - 1/(1-s) + dVb.*sb;
        H0           = c0.^(1-s)/(1-s) - 1/(1-s);
    
        % Indicator functions
        Iunique      = (sb<0).*(1-(sf>0)) + (1-(sb<0)).*(sf>0);
        Iboth        = (sb<0).*(sf>0);
        Ib           = Iunique.*(sb<0).*(Hb>H0) + Iboth.*(Hb>Hf).*(Hb>H0);
        If           = Iunique.*(sf>0).*(Hf>H0) + Iboth.*(Hf>Hb).*(Hf>H0);
        I0           = 1-Ib-If;
        c            = cf.*If + cb.*Ib + c0.*I0;
        u            = (c.^(1-s)-1)/(1-s);
        % Build the A matrix
        X = - min(ssb,0)/da;
        Y = - max(ssf,0)/da + min(ssb,0)/da;
        Z =   max(ssf,0)/da;
    
        % ******** Fix entries on A matrix ***********
        Y(1)   = Y(1) - min(adot(1,1)/da,0); % Y(1) + X(1)
        u(1,1) = u(1,1) + (cb(1,1))^(-s)*min(adot(1,1)/da,0);
        % ******** Fix entries on A matrix ***********
    
        A1 = spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
        A2 = spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
        A  = [A1,sparse(I,I);sparse(I,I),A2] + Aswitch;
    
        B         = (rho + 1/dt)*speye(2*I) - A;
        u_stacked = [u(:,1);u(:,2)];
        V_stacked = [V(:,1);V(:,2)];
        b         = u_stacked + V_stacked/dt;
        V_stacked = max(B\b,Vstar);
        V         = reshape(V_stacked,I,2);
        Vchange   = V - v;
        v         = V;
    
        if mod(n,10000)==0
            plot(a,V,'linewidth',3)
            grid on
            drawnow
        end
    
        dist(n) = max(max(abs(Vchange)));
        if dist(n)<crit
            disp('Splitting method value Function Converged, Iteration = ')
            disp(n)
            break
        end
    end

    VSM = V;
    nSM = n;
end

%% Nuno, Hurtado and Thomas' method

v     = vnodef;
lam   = 370;
Delta = 0.005;
clear('dist')
for n=1:maxitNM
    V            = v;
    % Forward difference
    dVf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVf(I,:)     = (z + r(end,:).*amax).^(-s);
    % Backward difference
    dVb(2:I,:)   = (V(2:I,:)-V(1:I-1,:))/da;
    dVb(1,:)     = (z + r(1,:).*amin).^(-s);

    dVf          = max(dVf,csubs);
    dVb          = max(dVb,csubs);
    % Consumption and savings with forward difference
    cf           = dVf.^(-1/s);
    ssf          = zz + r.*aa - cf;
    % Consumption and savings with backward difference
    cb           = dVb.^(-1/s);
    ssb          = zz + r.*aa - cb;
    % Consumption and derivative of value function at the temporary steady state
    c0           = zz + r.*aa;
    dV0          = c0.^(-s);
    % Hamiltonians
    Hf           = cf.^(1-s)/(1-s) - 1/(1-s) + dVf.*sf;
    Hb           = cb.^(1-s)/(1-s) - 1/(1-s) + dVb.*sb;
    H0           = c0.^(1-s)/(1-s) - 1/(1-s);
    % Indicator functions
    Iunique      = (sb<0).*(1-(sf>0)) + (1-(sb<0)).*(sf>0);
    Iboth        = (sb<0).*(sf>0);
    Ib           = Iunique.*(sb<0) + Iboth.*(Hb>Hf).*(Hb>H0);
    If           = Iunique.*(sf>0) + Iboth.*(Hf>Hb).*(Hf>H0);
    I0           = 1-Ib-If;
    c            = cf.*If + cb.*Ib + c0.*I0;
    u            = (c.^(1-s)-1)/(1-s);
    
    % Build the A matrix
    X            = - min(ssb,0)/da;
    Y            = - max(ssf,0)/da + min(ssb,0)/da;
    Z            =   max(ssf,0)/da;
    A1           = spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
    A2           = spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
    A            = [A1,sparse(I,I);sparse(I,I),A2] + Aswitch;
    B            = (rho + 1/Delta)*speye(2*I) - A;
    u_stacked    = u(:);
    V_stacked    = V(:);
    b            = u_stacked + V_stacked/Delta + lam*(Vstar>=V(:)).*(Vstar-V_stacked).*(aa(:)<0);
    V_stacked    = B\b;
    V            = reshape(V_stacked,I,2);
    Vchange      = V - v;
    v            = V;
    dist(n)      = max(abs(Vchange(:)));

    if mod(n,1000)==0
        subplot(121)
        plot(a,V(:,1),'linewidth',3)
        hold on
        plot(a,Vstar(1:I),'k--','linewidth',3)
        hold off
        grid on
        subplot(122)
        plot(log10(dist),'linewidth',3)
        grid on
        drawnow
    end
    
    if dist(n)<crit
        disp('NHT Value Function Converged, Iteration = ')
        disp(n)
        break
    end
end
VNHTM = V;
nNHTM = n;

%% Plots

if compute_VSM==0
    load('VSM')
end

% Get the default boundary
LCPindx          = abs(VLCP(1:I,1)-Vstar(1:I))<1e-6;
SMindx           = abs(VSM(1:I,1)-Vstar(1:I))<1e-6;
NHTMindx         = VNHTM(1:I,1)-Vstar(1:I)<1e-6;
bdry_LCP         = find(LCPindx,1,'last');
bdry_SM          = find(SMindx,1,'last');
bdry_NHTM        = find(NHTMindx,1,'last');
VLCP_Vstar_diff  = VLCP(:,1)-Vstarplot(:,1);
VSM_Vstar_diff   = VSM(:,1)-Vstarplot(:,1);
VNHTM_Vstar_diff = VNHTM(:,1)-Vstarplot(:,1);
% For plotting purposes
VLCP2            = VLCP;
VSM2             = VSM;
VNHTM2           = VNHTM;
VLCP2(LCPindx)   = NaN;
VSM2(SMindx)     = NaN;
VNHTM2(NHTMindx) = NaN;

subplot(121)
plot(a,Vstarplot(:,1),'linewidth',6)
hold on
plot(a,VLCP2(:,1),'linewidth',4)
plot(a,VSM2(:,1),'-.','linewidth',2)
plot(a,VNHTM2(:,1),'k--','linewidth',1)
scatter(a(bdry_LCP),Vstar(bdry_LCP),60,'k','filled')
hold off
grid on
xlabel('Wealth')
ylabel('Value functions')
legend('V^D','VLCP','VSM','VNHTM','Default boundary','location','Northwest')
legend('boxoff')
subplot(122)
plot(a,VLCP_Vstar_diff,'linewidth',4)
hold on
plot(a,VSM_Vstar_diff,'-.','linewidth',2)
plot(a,VNHTM_Vstar_diff,'k--','linewidth',1)
scatter(a(bdry_LCP),0,60,'k','filled')
hold off
grid on
xlabel('Wealth')
ylabel('V^C-V^D')
legend('VLCP','VSM','VNHTM','Default boundary','location','Northwest')
legend('boxoff')

clc
disp('Default boundary')
disp('      LCP      SM        NHTM')
disp([a(bdry_LCP) a(bdry_SM) a(bdry_NHTM)])
