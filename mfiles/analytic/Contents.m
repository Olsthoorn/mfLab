% ANALYTIC
%
% Files
%   Birsoy_Summers         - Variable-discharge well tests and tests in well fields
%   BruggemanMultilayerObj - class deff for objects of analytical solutions by Bruggeman (1999)
%   E                      - computes Hantush's modification of the Theis method for partially penetrating wells (ds_pp)
%   EdelmanFunc            - computes E_n=ierfc(u,n)/ierfc(0,n) used in 1D Edelman solutions
%   EdelmanQ               - computes transient 1D flow, constant profile
%   EdelmanS               - computes transient 1D flow, constant profile (Edelman solutions)
%   hantush                - computes Hantush's well function computed by integration
%   hantushE               - computes hantush but for partially penetrating well + self-test
%   hantushn               - Analytical N-layer solution of flow in semiconf multiple aquifer system
%   nsecn                  - solves multi layer flat analytical model (steady state)
%   shownsecn              - shows results of nsecn
%   solution               - class definition for objects of anlytical solutions for dynamics in 1D top layer
%   stehfest               - numerical back transformation from Laplace space according to Stehfest
%   TestHantushnScript     - script that verifies the analytical mutilayer nantushn
%   theis                  - Drawdown according to Theis computed as hantush(u,0) + selftest compares integration with expint.
%   TypeCurvesHantush      - produces Hantush type curves through function Wh(u,rho)
%   Wh                     - W computes Hantush's well function for flow to a well in a semi-confined aquifer
%   Wh1                    - computes Hantush's well function (Maas C, Veling, E, (2010), Stromingen, 16(2010)59-69
