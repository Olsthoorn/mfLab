function [STRTHD,p] = envHead(o,STRTHD,STCONC,drhodc,rho0)
%% STRTHD = gr.envHead(STRTHD,STCONC,drhodc,rho0);
%
% Compute pointwater head if environmental head is given.
% As an alternative to CHDDENSOPT and as its benchmark.
% Only the top head of STRTHD is used and unaltered.
%
% TO 120628 120704

if nargin<5, rho0 = 1000; end

[STRTHD,p] = o.hydrostaticPntwHead(STRTHD,STCONC,drhodc,rho0);
