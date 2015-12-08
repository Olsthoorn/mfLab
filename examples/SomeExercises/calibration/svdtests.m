% SVD excercises (normally done in Mathematica)
% www.uwlax.edu/faculty/will/svd/
% see pointobj in mfiles/etc
% TO 130625

close all;
p(1) = pointObj();

p(1) = p(1).circle(1,36); p(1).plot; view(2);

% R is the hanger, R' the aligner and S the stretcher matrix
% That is, S is a diagonal matrix, like singular value matrix.

a= pi/3;
R = [cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 0]  % hanger, R'=aligner
S = [2 0 0;  0 -0.5 0; 0 0 0]                   % stretcher
P = [2,1,0]';                                   % a point

% align, stretch and hang
p(2) = p(1).hit(R * S * R'); p(2).plot; p(2).connect(p(1)) % align stretch hang unit circle
p(3) = pointObj(P(1),P(2),P(3)); p(3).plot;           % plot P
p(4) = p(3).hit(R * S * R'); p(4).plot;  % align stretch and hang P

axis equal

