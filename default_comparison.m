% Bankruptcy choice - Comparing 3 methods
% Gustavo Mellior g.mellior@liverpool.ac.uk
% This code compares LCP, Splitting method (SM) and Hurtado, Nuno and Thomas' (2023) method (NHTM)

% set(groot, 'defaulttextInterpreter','latex') 
% set(groot, 'defaultAxesTickLabelInterpreter','latex') 
% set(groot, 'defaultAxesFontsize',14) 
% set(groot, 'defaultLegendInterpreter','latex')

clear; clc; close all

s           = 2;                      % CRRA utility
rho         = 0.05;                   % Discount rate
psi         = 0;                      % Mental cost of bankruptcy
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
r0          = 0.035;                  % Interest rate
utility     = @(x) (x.^(1-s))/(1-s);
maxit       = 100;                    % Max number of loops in V (when using LCP)
compute_VSM = 0;                      % Set to 1 if you want to compute V with the SM method (it takes a long time)
maxitSM     = 200000;
maxitHNT    = 200000;                 % Max number of loops in V (when using Hurtado, Nuno and Thomas)
crit        = 10^(-6);                % Set very small so that LCP error clears - see further below
Delta       = Inf;                    % Time step in LCP and no default cases
dVf         = zeros(I,2);             % Pre-allocate memory
dVb         = zeros(I,2);
c           = zeros(I,2);
Aswitch     = [-speye(I)*la(1),speye(I)*la(1);speye(I)*la(2),-speye(I)*la(2)]; % Income type transitions

% Bankruptcy value function
Vstartype = 3;
if Vstartype==1                             % Not flat V^A - Smooth pasting
    vnodef = load('vnodef');
    Vstar = vnodef.vnodef(:,1);
    Vstar = [Vstar-psi;Vstar-5000];
    r     = (exp(-(aa-amin))*0.1+r0).*(aa<0) + (aa>=0)*r0;
elseif Vstartype==2                         % Corner solution - V not flat
   Vstar  = [-23+1/5*a-psi;(-5000-psi)*ones(I,1)];
   r      = r0*ones(I,2);
else                                        % Corner solution - V flat
    Vstar = [(-23-psi)*ones(I,1);(-5000-psi)*ones(I,1)];
    r     = r0*ones(I,2);
end
Vstarplot       = reshape(Vstar,I,2);          % Just for plotting purposes
Vstarplot(aa>0) = NaN;                         % Just for plotting purposes
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
    dv          = (rho+la(1))*Vstar(1,1); % Multip denominator of V^C out, affects default value net of mental cost
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
    dist(n) = max(abs(Vchange(:)));
    % Convergence
    if (dist(n)<crit)&&(LCP_error<crit*1000) % Laxer to avoid cycles on LCPerror
        conv = 1;
         % disp('****** HJBVI SOLVED ******');
        break
    end

end
nLCP         = n;
timeLCP      = toc(t1);
VLCP         = V;
distLCP      = dist(n);
LCPindx      = abs(VLCP(1:I,1)-Vstar(1:I))<1e-6;
bdry_LCP     = find(LCPindx,1,'last');
VrelLCP      = max(abs(Vchange(:)./VLCP(:)));
temp0        = rho*V(:)- utility(c(:)) - A*V(:);
temp1        = V(:);
HJBerrLCP    = max(abs( temp0(bdry_LCP+1:end) ));
HJBerrRelLCP = max(abs( temp0(bdry_LCP+1:end)./temp1(bdry_LCP+1:end) ));

%% Splitting method
if Vstartype>1
    dt = 0.01;
else
    dt = 0.05;
end

plot(a,vnodef(:,1),a,V(:,1),'LineWidth',2)
hold on
plot(a,V(1:I,1),'LineWidth',2,'color',[0 0 0.6])
plot(a,Vstar(1:I),'k--','LineWidth',2)
hold off
grid
legend('$V^H$','$V^D$','$V_{LCP}$','Box','off','Location','southeast')
xlim([amin 2])
drawnow

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
        dv          = (rho+la(1))*Vstar(1,1); % Multip denominator of V^C out, affects default value net of mental cost
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
    
        A1 = spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
        A2 = spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
        A  = [A1,sparse(I,I);sparse(I,I),A2] + Aswitch;
        
        switch Vstartype
            case {2,3}
            if n>800
                dt = 0.00005;
            end
        end
        B         = (rho + 1/dt)*speye(2*I) - A;
        u_stacked = [u(:,1);u(:,2)];
        V_stacked = [V(:,1);V(:,2)];
        b         = u_stacked + V_stacked/dt;
        V_stacked = max(B\b,Vstar);
        V         = reshape(V_stacked,I,2);
        Vchange   = V - v;
        v         = V;
        dist(n)   = max(max(abs(Vchange)));
        % if mod(n,50000)==0
        %     figure(1)
        %     subplot(121)
        %     plot(a,vnodef(:,1),a,V(:,1),'LineWidth',2)
        %     hold on
        %     plot(a,Vstar(1:I),'k--','LineWidth',2)
        %     grid
        %     xlim([amin 2])
        %     plot(a,V(:,1),'-.','linewidth',3,'color',[0 0.5 0])
        %     SMindx  = abs(V(1:I,1)-Vstar(1:I))<1e-6;
        %     bdry_SM = find(SMindx,1,'last');
        %     scatter(a(bdry_LCP),VLCP(bdry_LCP),200,'filled','MarkerFaceColor',[0 0 0.6])
        %     scatter(a(bdry_SM),V(bdry_SM),100,'filled','MarkerFaceColor',[0 0.5 0])
        %     hold off
        %     legend('$V_{H}$','$V^D$','$V_{SM}$','Box','off','Location','southeast')
        %     subplot(122)
        %     plot(log10(dist),'linewidth',3)
        %     grid
        %     drawnow
        % end
    
        if dist(n)<crit
            disp('Splitting method value Function Converged, Iteration = ')
            disp(n)
            break
        end
    end
    timeSM      = toc(t2);
    VSM         = V;
    nSM         = n;
    SMindx      = abs(VSM(1:I,1)-Vstar(1:I))<1e-6;
    distSM      = dist(n);
    bdry_SM     = find(SMindx,1,'last');
    temp0       = rho*VSM(:)- utility(c(:)) - A*VSM(:);
    temp1       = VSM(:);
    HJBerrSM    = max(abs( temp0(bdry_SM+1:end) ));
    HJBerrRelSM = max(abs( temp0(bdry_SM+1:end)./temp1(bdry_SM+1:end) ));
else
    typeres     = strcat('VSMres',num2str(Vstartype));
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

%% Nuno, Hurtado and Thomas' method

v     = vnodef;
lam   = 370/2;
Delta = 0.005;
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
    B            = (rho + 1/Delta)*speye(2*I) - A;
    u_stacked    = u(:);
    V_stacked    = V(:);
    b            = u_stacked + V_stacked/Delta + lam*(Vstar>=V(:)).*(Vstar-V_stacked).*(aa(:)<0);
    V_stacked    = B\b;
    V            = reshape(V_stacked,I,2);
    Vchange      = V - v;
    v            = V;
    dist(n)      = max(abs(Vchange(:)));

    % if mod(n,2000)==0
    %     figure(1)
    %     subplot(121)
    %     plot(a,vnodef(:,1),a,V(:,1),'LineWidth',2)
    %     hold on
    %     plot(a,Vstar(1:I),'k--','LineWidth',2)
    %     grid
    %     xlim([amin 2])
    %     plot(a,V(:,1),'-.','linewidth',3,'color',[1 0 0])
    %     NHTMindx  = V(1:I,1)-Vstar(1:I)<1e-6;
    %     bdry_NHTM = find(NHTMindx,1,'last');
    %     scatter(a(bdry_LCP),VLCP(bdry_LCP,1),200,'filled','MarkerFaceColor',[0 0 0.6])
    %     scatter(a(bdry_SM),VSM(bdry_SM,1),100,'filled','MarkerFaceColor',[0 0.5 0])
    %     scatter(a(bdry_NHTM),V(bdry_NHTM,1),50,'filled','MarkerFaceColor',[1 0 0])
    %     hold off
    %     legend('$V_{H}$','$V^D$','$V_{NT}$','Box','off','Location','southeast')
    %     subplot(122)
    %     plot(log10(dist),'linewidth',3)
    %     grid
    %     drawnow
    % end
    
    if dist(n)<crit
        break
    end
end
timeHNT = toc(t3);
VNHNT   = V;
distHNT = dist(n);
nHNT    = n;
VrelHNT = max(abs(Vchange(:)./VNHNT(:)));

%% Results
% Get the default boundary

NHTMindx         = VNHNT(1:I,1)-Vstar(1:I)<1e-6;
bdry_HNT         = find(NHTMindx,1,'last');
VLCP_Vstar_diff  = VLCP(:,1)-Vstarplot(:,1);
VSM_Vstar_diff   = VSM(:,1)-Vstarplot(:,1);
VNHTM_Vstar_diff = VNHNT(:,1)-Vstarplot(:,1);
temp0            = rho*V(:)- utility(c(:)) - A*V(:);
temp1            = V(:);
HJBerrHNT        = max(abs( temp0(bdry_HNT+1:end) ));
HJBerrRelHNT     = max(abs( temp0(bdry_HNT+1:end)./temp1(bdry_HNT+1:end) ));
% For plotting purposes
VLCP2            = VLCP;
VSM2             = VSM;
VHNTM2           = VNHNT;
VLCP2(LCPindx)   = NaN;
VSM2(SMindx)     = NaN;
VHNTM2(NHTMindx) = NaN;

clc
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
% fprintf('\t\t\t\t\t %i\t %i\t %i\n', la, lb, la+lb);
fprintf('      LCP \t      SM \t        NHTM \t \n\n')
fprintf('\t %f\t %f\t %f\n', abs([HJBerrRelLCP HJBerrRelSM HJBerrRelHNT]))


%% Create figures in the paper

if Vstartype>1
    figure
    plot(a,vnodef(:,1),'color',0.6*ones(1,3),'linewidth',2)
    hold on
    plot(a,VLCP(:,1),'linewidth',2,'Color',[0 0 0.6])
    plot(a,Vstarplot(:,1),'k--','linewidth',2)
    plot(a(bdry_LCP)*ones(2,1),linspace(min(vnodef(:))*1.1,max(VLCP(:))*1.1,2),'k--')
    hold off
    grid on
    axis tight
    xlabel('$a$')
    ylabel('$V(a)$')
    legend('$V^H_L$','$V_{L}^N$','$V^D$','$a*=\underline{a}$','Default boundary','location','Southeast','Box','off')
    xlim([amin*1.1 2])
    if Vstartype==3
        saveas(gcf,'Figure1a','epsc')
    else
        saveas(gcf,'Figure2','epsc')
    end
else
    figure
    plot(a,vnodef(:,1),'color',0.6*ones(1,3),'linewidth',2)
    hold on
    plot(a,VLCP(:,1),'linewidth',2,'Color',[0 0 0.6])
    plot(a,Vstarplot(:,1),'k--','linewidth',2)
    plot(a(bdry_LCP)*ones(2,1),linspace(min(vnodef(:))*1.1,max(VLCP(:))*1.1,2),'k--')
    hold off
    grid on
    xlabel('$a$')
    ylabel('$V(a)$')
    legend('$V^H_L$','$V_{L}^N$','$V^D$','$a*>\underline{a}$','Default boundary','location','Southeast','Box','off')
    axis tight
    xlim([amin*1.1 2])
    saveas(gcf,'Figure1b','epsc')
end