function paths=mf_mybatchpaths
% MF_MYBATCHPATHS Sets paths when running script off line using batch()
%
% USAGE mf_mybatchpaths
%     job=batch('mf_adapt','PathDependencies',mf_mybachpaths);
%     wait(job); % to see when it ends, but then you loos control of your PC
%     job        % to see if it is ready.
%     load(job); % when ready
% When running a script in batch node, off line, it will run in its own
% workspace environment and does not know about paths set in the original
% one. Therefore we have to set the paths to the mfiles director of mflab
% and pass these path names onto the batch job's environment
%
% Notice that batch uses the parallel toolbox
%
% TO 110511

P='/Domain/tudelft.net/Users/tolsthoorn/GRWMODELS/mflab/';
paths=...
  {[P 'mfiles/read'];...
   [P 'mfiles/write'];...
   [P 'mfiles/etc'];...
   [P 'mfiles/visualization'];...
   [P 'mfiles/gridcoords'];...
   [P 'mfiles/fdm'];...
   [P 'mfiles/analytic'];...
   [P 'bin']};


