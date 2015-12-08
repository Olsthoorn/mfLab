function txt = lsqExit(exitNr)
% get reason for exiting lsqnonlin function
% TO 141010

exitReasons = {
1, 'Function converged to a solution x.';
2, 'Change in x was less than the specified tolerance';
3, 'Change in the residual was less than the specified tolerance.';
4, 'Magnitude of search direction was smaller than the specified tolerance.';
0, 'Number of iterations exceeded options.MaxIter or number of function evaluations exceeded options.FunEvals.';
-1,'Output function terminated the algorithm.';
-2,'Problem is infeasible: the bounds lb and ub are inconsistent.';
-3,'Regularization parameter became too large (levenberg-marquardt algorithm).';
-4,'Line search could not sufficiently decrease the residual along the current search direction.'
};

L = ismember([exitReasons{:,1}],exitNr);

txt = exitReasons(L,2);
txt = txt{1};