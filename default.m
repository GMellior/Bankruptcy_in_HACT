% Bankruptcy choice
% Gustavo Mellior and Katsuyuki Shibayama

clear all; clc; close all

s    = 2;    % CRRA utility
rho  = 0.05; % Discount rate
psi  = 2;    % Mental cost of bankruptcy
csubs= 1e-6; % Subsistence level of consumption - helps keep slope of V positive, see below
plotconv = 0;% Set to 1 if you want to see plots during loops - faster if set to zero
z1   = 0.75; % Low income flow
z2   = 1.25; % High income
z    = [z1,z2];
la1  = 0.25; % Poisson rate low to high
la2  = 0.25; % Poisson high to low
la   = [la1,la2];
I    = 800;  % Wealth grid points
amin = -4;   % Exogenous debt limit
amax = 5;    % Maximum level of wealth
a    = linspace(amin,amax,I)'; % Wealth grid
da   = (amax-amin)/(I-1);
aa   = [a,a];
zz   = ones(I,1)*z;

%***** Risk premium *****
% Choose zero for no risk premium - Affects whether bankruptcy takes place in the interior or at amin
risk_premium = 1; 
if risk_premium==1
    temp = (0.5 + exp(aa-0.15)).*(aa<0) + ones(I,2).*(aa>0);
    r = 0.035./temp;
else
    r = ones(I,2)*0.035;
end
%***** Slope of V^A *****
Vstarflat = 0;          % Set to 1 if you want autarky V^A to be flat
%
maxit     = 100;        % Max number of loops in V
crit      = 10^(-9);    % Set very small so that LCP error clears - see further below
Delta     = 1000;       % Time step
dVf       = zeros(I,2); % Pre-allocate memory
dVb       = zeros(I,2);
c         = zeros(I,2);
Aswitch   = [-speye(I)*la(1),speye(I)*la(1);speye(I)*la(2),-speye(I)*la(2)];
%INITIAL GUESS
v0(:,1)   = ((z(1) + r(:,1).*a).^(1-s) - 1)/(1-s)/rho;
v0(:,2)   = ((z(2) + r(:,1).*a).^(1-s) - 1)/(1-s)/rho;
v         = v0;

%% Basic Huggett - this will be used as a starting guess for HJBVI
for n=1:maxit
    V = v;
    % Forward difference
    dVf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVf(I,:)     = (z + r(end,:).*amax).^(-s);
    % Backward difference
    dVb(2:I,:)   = (V(2:I,:)-V(1:I-1,:))/da;
    dVb(1,:)     = (z + r(1,:).*amin).^(-s); % Huggett state constraint boundary condition
    
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
    
    % dV_upwind makes a choice of forward or backward differences based on
    % the sign of the drift    
    If = ssf > 0; %positive drift --> forward difference
    Ib = ssb < 0; %negative drift --> backward difference
    I0 = (1-If-Ib); %at steady state
    
    dV_Upwind = dVf.*If + dVb.*Ib + dV0.*I0;
    c = dV_Upwind.^(-1/s);
    u = (c.^(1-s)-1)/(1-s);
    
    %CONSTRUCT MATRIX
    X = - min(ssb,0)/da;
    Y = - max(ssf,0)/da + min(ssb,0)/da;
    Z =   max(ssf,0)/da;
    
    A1=spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
    A2=spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
    A = [A1,sparse(I,I);sparse(I,I),A2] + Aswitch;
 
    %Check transition matrix conditions 
    if max(abs(sum(A,2)))>10^(-12)
        disp('Improper Transition Matrix')
        break
    end    
    
    B = (rho + 1/Delta)*speye(2*I) - A;
    
    u_stacked = [u(:,1);u(:,2)];
    V_stacked = [V(:,1);V(:,2)];
    
    b = u_stacked + V_stacked/Delta;
    V_stacked = B\b; %SOLVE SYSTEM OF EQUATIONS
    
    V = [V_stacked(1:I),V_stacked(I+1:2*I)];
    
    Vchange = V - v;
    v = V;

    dist(n) = max(max(abs(Vchange)));
    if dist(n)<crit
        disp('Huggett value Function Converged, Iteration = ')
        disp(n)
        break
    end
end


%% HJBVI

% Solver options and initial guess for c(amin)
options   = optimset('Display','off','MaxIter', 10000, 'MaxFunEvals', 5000);
x0        = 0.05*z(1,1);
% Bankruptcy value
if Vstarflat==0 % Not flat V^A - arbitrarily generated value of bankruptcy
    Vstar = ((z(2)*1.15 + 0.025*a.*(a<0)).^(1-s) - 1)/(1-s)/rho;
    Vstar = [Vstar-psi;Vstar-5000];
else % Flat V^A
   Vstar = [(-0-psi)*ones(I,1);(-5000-psi)*ones(I,1)];
end
Vstarplot       = reshape(Vstar,I,2); % Just for plotting purposes
Vstarplot(aa>0) = NaN;                % Just for plotting purposes
LCPerror        = 1;                  % See related comments below
fct             = -1;                 % See comments below
%**** Start V^C vs V^A
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
    while fct<0 % Exploit positive slope of VMatching
        x0 = x0 + 0.1;
        myf1         = @(x) cTsolver(x,params);
        [a1, a2, a3] = fsolve(myf1, x0, options);
        fct          = s*a1.^(-s)-s*a1.^(-s-1)*driftminusC; % Derviative of VMatching wrt c
    end
    fct          = -1; % Reset fct for the next loop
    cT           = a1; % c(amin)
    cT           = min(cT,zz(1,1)*10); % We don't want strange cT, so let's give a bound to potential results
    cT           = max(cT,csubs);      % Might be unnecessary, but still keep it just in case
    if a2>0.0001 % If a2 is positive -> normal Huggett model with no bankruptcy
       cT = r(1,1)*amin +zz(1,1); 
    elseif (a3<0)&&(a2<0) % Rarely used
       cT = 1000;
    end
    
    % Forward difference
    dVf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVf(I,:)     = (z(end,:) + r(end,:).*amax).^(-s);
    % Backward difference
    dVb(2:I,:)   = (V(2:I,:)-V(1:I-1,:))/da;
    dVb(1,1)     = (max(cT,r(1,1)*amin + zz(1,1))).^(-s); % Only low income can default
    dVb(1,2)     = (z(2) + r(1,2).*amin).^(-s);

    dVf = max(dVf,csubs); % Helps with stability - guarantees slopes are non negative
    dVb = max(dVb,csubs);

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
    %CONSTRUCT THE A MATRIX
    X       = -Ib.*sb./da;
    Y       = -If.*sf./da + Ib.*sb./da;
    Z       = If.*sf./da;
    
    % ******** Fix entries on A matrix ***********
    Y(1)   = Y(1) - min(adot(1,1)/da,0); % Y(1) + X(1)
    u(1,1) = u(1,1) + (cb(1,1))^(-s)*min(adot(1,1)/da,0);
    % ******** Fix entries on A matrix ***********
    
    % A matrix
    A1=spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
    A2=spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
    A = [A1,sparse(I,I);sparse(I,I),A2] + Aswitch;
    
    %Check transition matrix conditions 
    if max(abs(sum(A,2)))>10^(-12)
        disp('Improper Transition Matrix')
        break
    end
    
    % B matrix
    B = (1/Delta + rho)*speye(2*I) - A;

    u_stacked = reshape(u,2*I,1);
    V_stacked = reshape(V,2*I,1);
    
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
    
    dist(n) = max(max(max(abs(Vchange))));
    % Plots
    if plotconv==1 
        subplot(1,2,1)
        plot(a,V,'linewidth',3)
        hold on
        plot(a,Vstarplot,'k-.','linewidth',3)
        hold off
        grid on
        title(sprintf('Vchange = %f',max(max(Vchange))))
        ylim([-4 7])
        xlabel('Wealth')
        ylabel('V^C and V^A')
        legend('V_L^C','V_H^C','V^A')
        legend('boxoff')
        subplot(1,2,2)
        plot(a,zeros(I,1),'k--')
        hold on
        adotnan = adot;
        adotnan(1:spindex(end),1) = NaN; 
        plot(a,adotnan,'linewidth',3)
        hold on
        scatter(a(spindex(end)),adot(spindex(end),1),100,'filled','MarkerFaceColor',[0 0 0.6])
        hold off
        xlim([amin amax])
        grid on
        title(sprintf('LCP error = %f',LCP_error))
        xlabel('Wealth')
        ylabel('Drift')
        legend('Zero drift','S_L^C','S_H^C','a^*')
        legend('boxoff')
        drawnow
    end
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

%% Plots
subplot(1,2,1)
plot(a,V,'linewidth',3)
hold on
plot(a,Vstarplot,'k-.','linewidth',3)
hold off
grid on
title(sprintf('Vchange = %f',max(max(Vchange))))
ylim([-4 7])
xlabel('Wealth')
ylabel('V^C and V^A')
legend('V_L^C','V_H^C','V^A','location','Southeast')
legend('boxoff')
subplot(1,2,2)
plot(a,zeros(I,1),'k--')
hold on
adotnan = adot;
adotnan(1:spindex(end),1) = NaN; 
plot(a,adotnan,'linewidth',3)
hold on
scatter(a(spindex(end)),adot(spindex(end),1),100,'filled','MarkerFaceColor',[0 0 0.6])
hold off
xlim([amin amax])
grid on
title(sprintf('LCP error = %f',LCP_error))
xlabel('Wealth')
ylabel('Drift')
legend('Zero drift','S_L^C','S_H^C','a^*','location','Southeast')
legend('boxoff')