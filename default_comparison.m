% Bankruptcy choice - Comparing 3 methods
% Gustavo Mellior g.mellior@liverpool.ac.uk
% This code compares LCP, Splitting method (SM) and Hurtado, Nuno and Thomas' (2023) method (NHTM)

set(groot, 'defaulttextInterpreter','latex') 
set(groot, 'defaultAxesTickLabelInterpreter','latex') 
set(groot, 'defaultAxesFontsize',14) 
set(groot, 'defaultLegendInterpreter','latex')

clear; clc; close all

saveresults = 0;                      % Save results
show_c_S    = 0;                      % Show consumption policy and drift
compute_VSM = 1;                      % Set to 1 if you want to compute V with the SM method (it can take a long time)
s           = 2;                      % CRRA utility
rho         = 0.05;                   % Discount rate
csubs       = 1e-6;                   % Subsistence level of consumption - helps keep slope of V positive, see below
plotconv    = 0;                      % Set to 1 if you want to see plots during loops - faster if set to zero
z1          = 0.75;                   % Low income
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
r0          = 0.035;                  % Interest rate (when a>0)
utility     = @(x) (x.^(1-s))/(1-s);
maxit       = 100;                    % Max number of loops in V (when using LCP)
maxitSM     = 2000000;
maxitHNT    = 200000;                 % Max number of loops in V (when using Hurtado, Nuno and Thomas)
crit        = 10^(-6);                % Stop criteria
dtLCP       = Inf;                    % Time step in LCP and no default cases
Delta       = dtLCP;                  % Update parameter in HNT
gamma       = 185/2;                  % Arrival rate of bankruptcy opportunities in HNT
dtHNT       = 0.02;                   % Update parameter in HNT
dVf         = zeros(I,2);             % Pre-allocate memory
dVb         = zeros(I,2);
c           = zeros(I,2);
Aswitch     = [-speye(I)*la(1),speye(I)*la(1);speye(I)*la(2),-speye(I)*la(2)]; % Income type transitions

% Bankruptcy value function
zd          = 0.9;
Vstartype   = 1;
if Vstartype==1                             % Interior solution
    psi     = 0.07;
elseif Vstartype==2                         % Corner solution - V not flat
    psi     = 0.001;
else                                        % Corner solution - V flat
    psi     = 0;
end
r           = 0.0075*exp(-2.7*(aa+3))+r0;   % Interest rate
VD          = [((((zd*1 + psi*r(:,1).*a).^(1-s))/(1-s))/rho).*(a<0)+(a>=0)*-22.24;-500*ones(I,1)];
if Vstartype>1
    dtSM    = 0.00001;
else
    dtSM    = 0.1;
end
implicitmethod = 1;                         % Set to 1 to solve SM and HNT with the implicit method.

Vstarplot       = reshape(VD,I,2);          % Just for plotting purposes
Vstarplot(aa>0) = NaN;                      % Just for plotting purposes
Cases           = ['A','B','C'];
% Initial guess of non-default V
v0(:,1)     = ((z(1) + r(:,1).*a).^(1-s) - 1)/(1-s)/rho;
v0(:,2)     = ((z(2) + r(:,1).*a).^(1-s) - 1)/(1-s)/rho;
v           = v0;

% Solver options and initial guess for c(amin)
options   = optimset('Display','off','MaxIter', 10000, 'MaxFunEvals', 5000);
x0        = 0.05*z(1,1);

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
    sf  = zz + r.*aa - cf;
    % Consumption and savings with backward difference
    cb  = dVb.^(-1/s);
    sb  = zz + r.*aa - cb;
    % Consumption and derivative of value function at the temporary steady state
    c0  = zz + r.*aa;
    
    % Upwind
    If        = sf > 0;
    Ib        = sb < 0;
    I0        = (1-If-Ib);
    c         = cf.*If + cb.*Ib + c0.*I0;
    u         = utility(c);
    % Build the A matrix
    X         = - min(sb,0)/da;
    Y         = - max(sf,0)/da + min(sb,0)/da;
    Z         =   max(sf,0)/da;
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
        break
    end
end

adotnodef = zz+r.*aa-c;
cnodef    = c;

%% LCP
LCPerror        = 1;                  % See related comments below
fct             = -1;                 % See comments below
vnodef          = v;                  % Save previous result where bankruptcy is not allowed
t1 = tic;
for n=1:maxit
    V = v;
    % Find c(amin)
    cT          = zeros(1,1);
    a3store     = cT;
    a2store     = cT;
    x0          = 0.05*z(1,1);            % Re-start this at every nth loop
    driftminusC = r(1,1)*amin+zz(1,1);    % Autarky cons
    psflow      = la(1)*V(1,2);           % Poisson flow in HJB
    dv          = (rho+la(1))*VD(1,1); % Multip denominator of V^C out, affects default value net of mental cost
    params      = [s driftminusC psflow 0 dv]; % Pass arguments to solver
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

    Hf = utility(cf) + dVf.*sf;
    Hb = utility(cb) + dVb.*sb;
    H0 = utility(c0);
    
    % More flexible approach allowing for non concavity
    Iunique = (sb<0).*(1-(sf>0)) + (1-(sb<0)).*(sf>0);
    Iboth   = (sb<0).*(sf>0);
    Ib      = Iunique.*(sb<0).*(Hb>H0) + Iboth.*(Hb>Hf).*(Hb>H0);
    If      = Iunique.*(sf>0).*(Hf>H0) + Iboth.*(Hf>Hb).*(Hf>H0);
    I0      = 1-Ib-If;
    
    % Use upwind scheme to get consumption
    c       = cf.*If + cb.*Ib + c0.*I0;
    u       = utility(c);
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
    B             = (1/dtLCP + rho)*speye(2*I) - A;
    
    u_stacked     = reshape(u,2*I,1);
    V_stacked     = reshape(V,2*I,1);
    % LCP
    vec           = u_stacked + V_stacked/dtLCP;
    qlcp          = -vec + B*VD; 
    z0            = V_stacked-VD;
    zans          = LCP(B,qlcp,zeros(2*I,1),Inf*ones(2*I,1),z0,0);
    LCP_error     = max(abs(zans.*(B*zans + qlcp)));
    LCPerror(n+1) = LCP_error; % Track error - sometimes it gets stuck and slows down things
    
    V_stacked     = zans+VD; % Get new value function
    V             = reshape(V_stacked,I,2);
    Vchange       = V - v;
    v             = V;
    dist(n) = max(abs(Vchange(:)));
    % Convergence
    if (dist(n)<crit)&&(LCP_error<crit*1000) % Laxer to avoid cycles on LCPerror
        conv = 1;
         disp('LCP method value function converged');
        break
    end

end
nLCP         = n;
timeLCP      = toc(t1);
VLCP         = V;
distLCP      = dist(n);
LCPindx      = abs(VLCP(1:I,1)-VD(1:I))<1e-6;
bdry_LCP     = find(LCPindx,1,'last');
VrelLCP      = max(abs(Vchange(:)./VLCP(:)));
temp0        = rho*V(:)- utility(c(:)) - A*V(:);
temp1        = V(:);
HJBerrLCP    = max(abs( temp0(bdry_LCP+1:end) ));
HJBerrRelLCP = max(abs( temp0(bdry_LCP+1:end)./temp1(bdry_LCP+1:end) ));
adotLCP      = zz+r.*aa-c;
cLCP         = c;

%% Splitting method

v = vnodef;
t2 = tic;
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
        dv          = (rho+la(1))*VD(1,1); % Multip denominator of V^C out, affects default value net of mental cost
        params      = [s driftminusC psflow 0 dv]; % Pass arguments to solver
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
        sf           = zz + r.*aa - cf;
        cb           = dVb.^(-1/s);
        sb           = zz + r.*aa - cb;
        c0           = zz + r.*aa;
        % Hamiltonians
        Hf           = utility(cf) + dVf.*sf;
        Hb           = utility(cb) + dVb.*sb;
        H0           = utility(c0);
    
        % Indicator functions
        Iunique      = (sb<0).*(1-(sf>0)) + (1-(sb<0)).*(sf>0);
        Iboth        = (sb<0).*(sf>0);
        Ib           = Iunique.*(sb<0).*(Hb>H0) + Iboth.*(Hb>Hf).*(Hb>H0);
        If           = Iunique.*(sf>0).*(Hf>H0) + Iboth.*(Hf>Hb).*(Hf>H0);
        I0           = 1-Ib-If;
        c            = cf.*If + cb.*Ib + c0.*I0;
        u            = utility(c);
        adot         = sf.*If+sb.*Ib;
        % Build the A matrix
        X = - min(sb,0)/da;
        Y = - max(sf,0)/da + min(sb,0)/da;
        Z =   max(sf,0)/da;
    
        % ******** Fix entries on A matrix ***********
        Y(1)   = Y(1) - min(adot(1,1)/da,0); % Y(1) + X(1)
        u(1,1) = u(1,1) + (cb(1,1))^(-s)*min(adot(1,1)/da,0);
        % ******** Fix entries on A matrix ***********
    
        A1        = spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
        A2        = spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
        A         = [A1,sparse(I,I);sparse(I,I),A2] + Aswitch;
        u_stacked = [u(:,1);u(:,2)];
        V_stacked = [V(:,1);V(:,2)];
        switch implicitmethod
            case 0
                V_stacked = max(V_stacked + dtSM*(u_stacked + A*V(:)-rho*V_stacked),VD);
            case 1
                B         = (rho + 1/dtSM)*speye(2*I) - A;
                b         = u_stacked + V_stacked/dtSM;
                V_stacked = max(B\b,VD); 
        end
        V         = reshape(V_stacked,I,2);
        Vchange   = V - v;
        v         = V;
        dist(n)   = max(max(abs(Vchange)));

        if dist(n)<crit
            disp('Splitting method value function converged')
            break
        end
    end
    timeSM      = toc(t2);
    VSM         = V;
    nSM         = n;
    SMindx      = abs(VSM(1:I,1)-VD(1:I))<1e-6;
    distSM      = dist(n);
    bdry_SM     = find(SMindx,1,'last');
    temp0       = rho*VSM(:)- utility(c(:)) - A*VSM(:);
    temp1       = VSM(:);
    HJBerrSM    = max(abs( temp0(bdry_SM+1:end) ));
    HJBerrRelSM = max(abs( temp0(bdry_SM+1:end)./temp1(bdry_SM+1:end) ));
else
    typeres     = strcat('resultsCase',Cases(Vstartype));
    VSMres      = load(typeres);
    timeSM      = VSMres.timeSM;
    VSM         = VSMres.VSM;
    nSM         = VSMres.nSM;
    SMindx      = VSMres.SMindx;
    distSM      = VSMres.distSM;
    bdry_SM     = VSMres.bdry_SM;
    HJBerrSM    = VSMres.HJBerrSM;
    HJBerrRelSM = VSMres.HJBerrRelSM;
end
adotSM          = zz+r.*aa-c;
cSM             = c;

%% Hurtado, Nuno and Thomas' method
v     = vnodef;
clear('dist')
t3    = tic;
for n=1:maxitHNT
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
    sf           = zz + r.*aa - cf;
    % Consumption and savings with backward difference
    cb           = dVb.^(-1/s);
    sb           = zz + r.*aa - cb;
    % Consumption and derivative of value function at the temporary steady state
    c0           = zz + r.*aa;
    % Hamiltonians
    Hf           = utility(cf) + dVf.*sf;
    Hb           = utility(cb) + dVb.*sb;
    H0           = utility(c0);
    % Indicator functions
    Iunique      = (sb<0).*(1-(sf>0)) + (1-(sb<0)).*(sf>0);
    Iboth        = (sb<0).*(sf>0);
    Ib           = Iunique.*(sb<0) + Iboth.*(Hb>Hf).*(Hb>H0);
    If           = Iunique.*(sf>0) + Iboth.*(Hf>Hb).*(Hf>H0);
    I0           = 1-Ib-If;
    c            = cf.*If + cb.*Ib + c0.*I0;
    u            = utility(c);
    
    % Build the A matrix
    X            = - min(sb,0)/da;
    Y            = - max(sf,0)/da + min(sb,0)/da;
    Z            =   max(sf,0)/da;
    A1           = spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
    A2           = spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
    A            = [A1,sparse(I,I);sparse(I,I),A2] + Aswitch;
    u_stacked    = u(:);
    V_stacked    = V(:);
    switch implicitmethod
        case 0
            V_stacked = max(V_stacked + dtHNT*(u_stacked + A*V(:)-rho*V_stacked + gamma*(VD>=V(:)).*(VD-V_stacked).*(aa(:)<0)),VD);
        case 1
            B            = (rho + 1/dtHNT)*speye(2*I) - A;
            b            = u_stacked + V_stacked/dtHNT + gamma*(VD>=V(:)).*(VD-V_stacked).*(aa(:)<0);
            V_stacked    = B\b;
    end
    V            = reshape(V_stacked,I,2);
    Vchange      = V - v;
    v            = V;
    dist(n)      = max(abs(Vchange(:)));
    
    if dist(n)<crit
        disp('Hurtado, Nuno and Thomas method value function converged')
        disp(' ')
        break
    end
end
timeHNT          = toc(t3);
VHNT             = V;
distHNT          = dist(n);
nHNT             = n;
NHTMindx         = VHNT(1:I,1)-VD(1:I)<1e-6;
bdry_HNT         = find(NHTMindx,1,'last');
temp0            = rho*V(:)- utility(c(:)) - A*V(:) - gamma*(VD>=V(:)).*(VD-V_stacked).*(aa(:)<0);
temp1            = V(:);
HJBerrHNT        = max(abs( temp0(bdry_HNT+1:end) ));
HJBerrRelHNT     = max(abs( temp0(bdry_HNT+1:end)./temp1(bdry_HNT+1:end) ));
adotHNT          = zz+r.*aa-c;
cHNT             = c;
%% Results
VLCP_Vstar_diff  = VLCP(:,1)-Vstarplot(:,1);
VSM_Vstar_diff   = VSM(:,1)-Vstarplot(:,1);
VNHTM_Vstar_diff = VHNT(:,1)-Vstarplot(:,1);

% For plotting purposes
VLCP2               = VLCP;
VSM2                = VSM;
VHNT2               = VHNT;
cLCP2               = cLCP;
cSM2                = cSM;
cHNT2               = cHNT;
adotLCP2            = adotLCP;
adotSM2             = adotSM;
adotHNT2            = adotHNT;
VLCP2(LCPindx)      = NaN;
VSM2(SMindx)        = NaN;
VHNT2(NHTMindx)     = NaN;
cLCP2(LCPindx)      = NaN;
cSM2(SMindx)        = NaN;
cHNT2(NHTMindx)     = NaN;
adotLCP2(LCPindx)   = NaN;
adotSM2(SMindx)     = NaN;
adotHNT2(NHTMindx)  = NaN;

disp('Default boundary')
disp('      LCP      SM        NHTM')
disp([a(bdry_LCP) a(bdry_SM) a(bdry_HNT)])
disp('')
disp('Convergence')
fprintf('      LCP \t      SM \t        NHTM \t \n\n')
disp([distLCP distSM distHNT])

disp('')
disp('HJB equation errors')
disp('')
disp('Absolute error')
fprintf('      LCP \t      SM \t        NHTM \t \n\n')
disp([HJBerrLCP HJBerrSM HJBerrHNT])
disp('')
disp('Relative error')
fprintf('      LCP \t      SM \t        NHTM \t \n\n')
fprintf('\t %f\t %f\t %f\n', [HJBerrRelLCP HJBerrRelSM HJBerrRelHNT])


%% Create figures in the paper

if Vstartype>1
    if Vstartype==3
        figure
        plot(a,vnodef(:,1),'color',0.6*ones(1,3),'linewidth',2)
        hold on
        plot(a,VLCP(:,1),'linewidth',2,'Color',[0 0 0.6])
        plot(a,Vstarplot(:,1),'k--','linewidth',2)
        plot(a(bdry_LCP)*ones(2,1),linspace(min(vnodef(:))*1.025,max(VLCP(:))*1.1,2),'k--')
        hold off
        grid on
        axis tight
        xlabel('$a$')
        ylabel('$V(a)$')
        legend('$V_L$','$V_{L}^N$','$V^D$','$a^*=\underline{a}$','location','Southeast','Box','off')
        xlim([amin*1.1 2])
        saveas(gcf,'Figure1a','epsc')
    else
        figure
        plot(a,VLCP(:,1),'linewidth',2,'Color',[0 0 0.6])
        hold on
        plot(a,Vstarplot(:,1),'k--','linewidth',2)
        plot(a(bdry_LCP)*ones(2,1),linspace(min(vnodef(:))*1.025,max(VLCP(:))*1.1,2),'k--')
        hold off
        grid on
        axis tight
        xlabel('$a$')
        ylabel('$V(a)$')
        legend('$V_{L}^N$','$V^D$','$a^*=\underline{a}$','location','Southeast','Box','off')
        xlim([amin*1.01 0.9375*amin])
        ylim([Vstarplot(1)*1.00015 Vstarplot(1)*0.9995])
        yticks(linspace(Vstarplot(1)*1.00015, Vstarplot(1)*0.9995, 4));
        saveas(gcf,'Figure2','epsc')
    end
else
    figure
    plot(a,vnodef(:,1),'color',0.6*ones(1,3),'linewidth',2)
    hold on
    plot(a,VLCP(:,1),'linewidth',2,'Color',[0 0 0.6])
    plot(a,Vstarplot(:,1),'k--','linewidth',2)
    plot(a(bdry_LCP)*ones(2,1),linspace(min(vnodef(:))*1.025,max(VLCP(:))*1.1,2),'k--')
    hold off
    grid on
    xlabel('$a$')
    ylabel('$V(a)$')
    legend('$V_L$','$V_{L}^N$','$V^D$','$a^*>\underline{a}$','location','Southeast','Box','off')
    axis tight
    xlim([amin*1.1 2])
    saveas(gcf,'Figure1b','epsc')
end

if saveresults==1
    eval(strcat('save resultsCase',Cases(Vstartype)));
end

if show_c_S==1
    figure
    subplot(121)
    plot(a,cnodef(:,1),'color',0.6*ones(1,3),'linewidth',2)
    hold on
    plot(a,cLCP2(:,1),'linewidth',4,'Color',[0 0 0.6])
    plot(a,cSM2(:,1),'linewidth',3 ,'Color',[0.8 0.898 1])
    plot(a,cHNT2(:,1),'o'          ,'Color',[0 0 0])
    plot(a(bdry_LCP)*ones(2,1),linspace(min(cnodef(:))*1.025,max(max(cnodef(:),cLCP2(:)))*1.1,2),'k--')
    hold off
    grid on
    xlabel('$a$');ylabel('$c(a)$')
    legend('$c_L$','$c_{L}^N$ LCP','$c_{L}^N$ SM','$c_{L}^N$ HNT','$a^*$','location','Southeast','Box','off')
    axis tight
    xlim([amin*1.1 2])
    subplot(122)
    plot(a,adotnodef(:,1),'color',0.6*ones(1,3),'linewidth',2)
    hold on
    plot(a,adotLCP2(:,1),'linewidth',2,'Color',[0 0 0.6])
    plot(a,adotSM2(:,1),'linewidth',2 ,'Color',[0.8 0.898 1])
    plot(a,adotHNT2(:,1),'o','Color',[0 0 0])
    plot(a(bdry_LCP)*ones(2,1),linspace(min(adotLCP2(:))*1.025,max(adotLCP2(:))*1.1,2),'k--')
    hold off
    grid on
    xlabel('$a$');ylabel('$S(a)$')
    legend('$S_L$','$S_{L}^N$ LCP','$S_{L}^N$ SM','$S_{L}^N$ HNT','$a^*$','location','Southeast','Box','off')
    axis tight
    xlim([amin*1.1 2])
end