function writeNWT(basename,nwt)
%WRITENWT writes input for MODFLOW-NWT input package
%
% Example
%    writeNWT(basename,nwt) --- write NWT file
%
% TO 130401

fid=fopen([basename,'.',nwt.ext],'wt');

%0.
fprintf(fid,'%s\n',['# MATLAB writeNWT ' datestr(now)]);
fprintf(    '%s\n',['# MATLAB writeNWT ' datestr(now)]);

fprintf(fid,'%g  %g  %d  %g  %d  %d  %d  ',...
    nwt.HEADTOL, nwt.FLUXTOL, nwt.MAXITEROUT, nwt.THICKFACT,...
    nwt.LINMETH, nwt.IPRNWT, nwt.IBOTAV);

switch nwt.OPTIONS
    case 1, fprintf(fid,' SIMPLE');
    case 2, fprintf(fid,' MODERATE');
    case 3, fprintf(fid,' COMPLEX');
    case 4, fprintf(fid,' SPECIFIED');
    otherwise
        error(['%s: unknown option <<%s>>, must be 1,2,3 or 4\n',...
               'for SIMPLE, MODERATE COMPLEX and SPECIFIED'],mfilename,nwt.OPTIONS);
end

%fprintf(fid,'\n');
%fprintf(fid,'     HEADTOL FLUXTOL MAXITEROUT THICKFACT LINMETH IPRNWT IBOTAV OPTIONS\n');

% HEADTOL (units of length)?is the maximum head change between outer iterations
%          for solution of the nonlinear problem (real).
% FLUXTOL (units of length cubed per time)?is the maximum root-mean-squared
%          flux difference between outer iterations for solution of the
%          nonlinear problem (real).
% MAXITEROUT is the maximum number of iterations to be allowed for solution
%          of the outer (nonlinear) problem (integer).
% THICKFACT is the portion of the cell thickness (length) used for smoothly
%          adjusting storage and conductance coefficients to zero (symbol ?
%          in equation 9; real).
% LINMETH is a flag that determines which matrix solver will be used.
%          A value of 1 indicates GMRES will be used and a value of 2
%          indicates ?MD will be used (integer).
% IPRNWT is a flag that indicates whether additional information about solver
%          convergence will be printed to the main listing file (integer).
% IBOTAV is a flag that indicates whether corrections will be made to
%          groundwater head relative to the cell-bottom altitude if the cell
%          is surrounded by dewatered cells (integer).
%          A value of 1 indicates that a correction will be made and a value
%          of 0 indicates no correction will be made. This input variable is
%          problem specific and both options (IBOTAV=0 or 1) should be tested.
% OPTIONS  are keywords that activate options:
%    SPECIFIED indicates that the optional solver input values listed for
%              items 1 and 2 will be specified in the NWT input file by the user.
%    SIMPLE    indicates that default solver input values will be defined
%              that work well for nearly linear models. This would be used
%              for models that do not include nonlinear stress packages, and
%              models that are either confined or consist of a single unconfined
%              layer that is thick enough to contain the water table within a
%              single layer. (See table 2 for the solver input values that will
%              be used for this option.)
%    MODERATE  indicates that default solver input values will be defined that
%              work well for moderately nonlinear models. This would be used
%              for models that include nonlinear stress packages, and models
%              that consist of one or more unconfined layers. The ?MODERATE?
%              option should be used when the ?SIMPLE? option does not result
%              in successful convergence. (See table 2 for the solver input
%              values that will be used for this option.)
%   COMPLEX    indicates that default solver input values will be defined that
%              work well for highly nonlinear models. This would be used for
%              models that include nonlinear stress packages, and models that
%              consist of one or more unconfined layers representing complex
%              geology and sw/gw interaction. The ?COMPLEX? option should be
%              used when the ?MODERATE? option does not result in successful
%              convergence. (See table 2 for the solver input values that will
%              be used for this option.)

if nwt.OPTIONS==4

    fprintf(fid,' %.g  %.g  %.g  %.g  %d',...
       nwt.DBDTHETA, nwt.DBDKAPPA, nwt.DBDGAMMA, nwt.MOMFACT, nwt.BACKFLAG);
   
   if nwt.BACKFLAG>0
       fprintf('  %d  %.g  %.g',...
           nwt.MAXBACKITER, nwt.BACKTOL, nwt.BACKREDUCE);
      %   fprintf(fid,'     DBDTHETA DBDKAPPA DBDGAMMA MOMFACT BACKFLAG MAXBACKITER BACKTOL BACKREDUCE\n');
   end
   fprintf(fid,'\n');
%% Read the following values if OPTIONS = ?SPECIFIED.?
% DBDTHETA    is a coefficient used to reduce the weight applied to the head change
%             between nonlinear iterations (symbol ? in equation 21).
%             DBDTHETA is used to control oscillations in head. Values range
%             between 0.0 and 1.0, and larger values increase the weight
%             (decrease under-relaxation) applied to the head change (real).
% DBDKAPPA    is a coefficient used to increase the weight applied to the head
%             change between nonlinear iterations (symbol ? in equation 22).
%             DBDKAPPA is used to control oscillations in head. Values range
%             between 0.0 and 1.0, and larger values increase the weight applied
%             to the head change (real).
% DBDGAMMA    is a factor (symbol ? in equation 19) used to weight the head 
%             change for iterations n?1 and n. Values range between 0.0 and 1.0,
%             and greater values apply more weight to the head change calculated
%             during iteration n (real).
% MOMFACT     is the momentum coefficient m of equation 20 and ranges between 0.0
%             and 1.0. Greater values apply more weight to the head change for
%             iteration n (real).
% BACKFLAG    is a flag used to specify whether residual control will be used.
%             A value of 1 indicates that residual control is active and a value
%             of 0 indicates residual control is inactive (integer).
% MAXBACKITER is the maximum number of reductions (backtracks) in the head change
%             between nonlinear iterations (integer). A value between 10 and 50 works well.
% BACKTOL     is the proportional decrease in the root-mean-squared error of the groundwater-flow equation used to determine if residual control is required at the end of a nonlinear iteration, as applied in equation 23 (real).
% BACKREDUCE  is a reduction factor (symbol Br in equation 23) used for residual
%             control that reduces the head change between nonlinear iterations (real).
%             Values should be between 0.0 and 1.0, where smaller values result
%             in smaller head-change values.
 
end

if nwt.OPTIONS==4 && nwt.LINMETH == 1
       fprintf(fid,'%d  %d  %d  %g  %d',...
           nwt.MAXITINNER, nwt.ILUMETHOD, nwt.LEVFILL, nwt.STOPTOL, nwt.MSDR);
       fprintf(fid,'\n');
       %fprintf(fid,'     MAXITINNER ILUMETHOD LEVFILL STOPTOL MSDR');
       
%% Read the following values if LINMETH = 1 and OPTIONS = ?SPECIFIED.?
% MAXITINNER is the maximum number of iterations for the linear solution (integer).
% ILUMETHOD is the index for selection of the method for incomplete factorization (ILU)
%           used as a preconditioner. See Kipp and others (2008) for further details (integer).
%      ILUMETHOD=1?ILU with drop tolerance and fill limit. Fill-in terms less than
%           drop tolerance times the diagonal are discarded. The number of fill-in terms
%           in each row of L and U is limited to the fill limit. The fill-limit largest
%           elements are kept in the L and U factors.
%      ILUMETHOD=2 ? ILU(k), Order k incomplete LU factorization. Fill-in terms of higher
%           order than k in the factorization are discarded.
% LEVFILL   is the fill limit for ILUMETHOD = 1 and is the level of fill for ILUMETHOD = 2.
%           Recommended values: 5-10 for method 1, 0-2 for method 2. See Kipp and others (2008)
%           for further details (integer).
% STOPTOL   is the tolerance for convergence of the linear solver. This is the residual
%           of the linear equations scaled by the norm of the root mean squared error.
%           Usually 10-8 to 10-12 works well. See Kipp and others (2008) for further
%           details (integer).
% MSDR      is the number of iterations between restarts of the GMRES Solver.
%           See Kipp and others (2008) for further details (integer).

end

if nwt.OPTIONS==4 && nwt.LINMETH == 2
    
       fprintf(fid,'%d  %d  %d  %d  %d  %g  %d  %g  %g  %d',...
           nwt.IACL, nwt.NORDER, nwt.LEVEL, nwt.NORTH, nwt.IREDSYS,...
           nwt.RRCTOLS, nwt.IDROPTOL, nwt.EPSRN, nwt.HCLOSEXMD, nwt.MXITERXMD);
       fprintf(fid,'\n');
%       fprintf(fid,'     IACL NORDER LEVEL NORTH IREDSYS RRCTOLS IDROPTOL EPSRN HCLOSEXMD MXITERXMD\n');

%% Read the following values if LINMETH = 2 and OPTIONS=?SPECIFIED.?
% IACL      is a flag for the acceleration method: 0 = conjugate gradient;
%           1 = ORTHOMIN; 2 = Bi-CGSTAB (integer).
% NORDER    is a flag for the scheme of ordering the unknowns: 0= original ordering;
%           1= RCM ordering; 2= Minimum Degree ordering (integer).
% LEVEL     is the level of fill for incomplete LU factorization (integer).
% NORTH     is the number of orthogonalization for the ORTHOMIN acceleration scheme.
%           A number between 4 and 10 is appropriate. Small values require less
%           storage but more iterations may be required. This number should
%           equal 2 for the other acceleration methods (integer).
% IREDSYS   is a flag for reduced system preconditioning: =1 apply reduced system
%           preconditioning; = 0 do not apply reduced system preconditioning (integer).
% RRCTOLS   is the residual reduction-convergence criteria (real).
% IDROPTOL  is a flag for using drop tolerance in the preconditioning (integer).
% EPSRN     is the drop tolerance for preconditioning (real).
% HCLOSEXMD is the head closure criteria for inner (linear) iterations (real).
% MXITERXMD is the maximum number of iterations for the linear solution.
end          

fclose(fid);
