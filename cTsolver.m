function F = cTsolver(x,params)
CRRA     = params(1,1); % CRRA
drift    = params(:,2); % Drift minus C
psflow   = params(:,3); % Transition poor to rich
constant = params(1,4); % Constant of CRRA utility function
defval   = params(:,5); % (rho+la1)*(V^A - mentalcost)

F(:,1) = (CRRA)/(1-CRRA)*(x.^(1-CRRA)) - constant/(1-CRRA) + (x.^(-CRRA)).*(drift) + psflow - defval;

end