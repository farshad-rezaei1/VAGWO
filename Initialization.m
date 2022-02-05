%_________________________________________________________________________%
% An Enhanced Grey Wolf Optimizer with a Velocity-Aided Global Search     %
% Mechanism                                                               %
%                                                                         %
% Developed in MATLAB R2018b                                              %
%                                                                         %
% Author and programmer: Farshad Rezaei, PhD                              %
%                                                                         %
% e-Mail: farshad.rezaei@gmail.com                                        %
%         f.rezaei@alumni.iut.ac.ir                                       %
%                                                                         %
% Homepage: https://www.linkedin.com/in/farshad-rezaei-5a92559a/          %
%                                                                         %
% Main paper: Rezaei, F.; Safavi, H.R.; Abd Elaziz, M.; El-Sappagh,       %
% S.H.A.; Al-Betar, M.A.; Abuhmed, T. An Enhanced Grey Wolf Optimizer     %
% with a Velocity-Aided Global Search Mechanism. Mathematics 2022, 10,    %
% 351.https://doi.org/10.3390/math10030351                                %
%_________________________________________________________________________%

% This function is to initialize the position and velocity of the wolves to start the optimization process
function [pp,gv_alpha,gv_beta,gv_delta] = Initialization(np,nx,varmax,varmin,velmax,velmin)
pp=zeros(np,nx); 
gv_alpha=zeros(np,nx);
gv_beta=zeros(np,nx);
gv_delta=zeros(np,nx);
for j=1:np
    pp(j,1:nx)=(varmax-varmin).*rand(1,nx)+varmin;
    gv_alpha(j,1:nx)=(velmax-velmin).*rand(1,nx)+velmin;
    gv_beta(j,1:nx)=(velmax-velmin).*rand(1,nx)+velmin;
    gv_delta(j,1:nx)=(velmax-velmin).*rand(1,nx)+velmin;
end