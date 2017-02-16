clear

% problem setting
nObj = 20;  % # objects
nPts = 20;   % # points on each object
p_ob = 0.7;  % observation ratio
p_error = 0.5;  % input error (the ratio of observations corrupted)

% synthesize data
OPT = problem_generator2(nObj,nPts,p_ob,p_error); % OPT.X_gt is the ground truth
X_gt = OPT.X_gt(2:end,2:end);
X_in = OPT.W(2:end,2:end);
dimGroup = OPT.stateDims;

% run multiple matching
X = mmatch_CVX_ALS(X_in,dimGroup);

% recovery error
norm(X-X_gt)/norm(X_gt)






