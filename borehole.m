function [y] = Borehole_Model(xx)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% BOREHOLE FUNCTION
%
% Authors: Sonja Surjanovic, Simon Fraser University
%          Derek Bingham, Simon Fraser University
% Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
%
% Copyright 2013. Derek Bingham, Simon Fraser University.
%
% THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
% FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
% derivative works, such modified software should be clearly marked.
% Additionally, this program is free software; you can redistribute it 
% and/or modify it under the terms of the GNU General Public License as 
% published by the Free Software Foundation; version 2.0 of the License. 
% Accordingly, this program is distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
% of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
% General Public License for more details.
%
% For function details and reference information, see:
% http://www.sfu.ca/~ssurjano/
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OUTPUT AND INPUT:
%
% y  = water flow rate
% xx = [rw, r, Tu, Hu, Tl, Hl, L, Kw]
%For the purposes of uncertainty quantification, the distributions of the input random variables are: 
%rw ~ N(?=0.10, ?=0.0161812) 
%r ~ Lognormal(?=7.71, ?=1.0056)    
%Tu ~ Uniform[63070, 115600] 
%Hu ~ Uniform[990, 1110] 
%Tl ~ Uniform[63.1, 116] 
%Hl ~ Uniform[700, 820] 
%L ~ Uniform[1120, 1680] 
%Kw ~ Uniform[9855, 12045] even [1500, 15 000] to be more challenging

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rw = xx(1);
r  = xx(2);
Tu = xx(3);
Hu = xx(4);
Tl = xx(5);
Hl = xx(6);
L  = xx(7);
Kw = xx(8);

frac1 = 2 * pi * Tu * (Hu-Hl);

frac2a = 2*L*Tu / (log(r/rw)*rw^2*Kw);
frac2b = Tu / Tl;
frac2 = log(r/rw) * (1+frac2a+frac2b);

y = frac1 / frac2;

end
