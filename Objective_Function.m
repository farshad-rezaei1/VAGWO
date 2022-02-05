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

% This function containts full information and implementations of the benchmark 
% functions used as the test bed in the main paper

% lb is the lower bound: lb=[lb_1,lb_2,...,lb_d]
% up is the upper bound: ub=[ub_1,ub_2,...,ub_d]
% dim is the number of the decision variables (dimensionality of the problem)
% so:shifted optimum

function [lb,ub,dim,fobj] = Objective_Function(F)

switch F
    case 'F1'
        fobj = @F1;
        lb=-100;
        ub=100;
        dim=100;
        
    case 'F2'
        fobj = @F2;
        lb=-10;
        ub=10;
        dim=100;
        
    case 'F3'
        fobj = @F3;
        lb=-100;
        ub=100;
        dim=100;
        
    case 'F4'
        fobj = @F4;
        lb=-100;
        ub=100;
        dim=100;
        
    case 'F5'
        fobj = @F5;
        lb=-30;
        ub=30;
        dim=100;
        
    case 'F6'
        fobj = @F6;
        lb=-100;
        ub=100;
        dim=100;
        
    case 'F7'
        fobj = @F7;
        lb=-1.28;
        ub=1.28;
        dim=100;
        
    case 'F8'
        fobj = @F8;
        lb=-500;
        ub=500;
        dim=100;
        
    case 'F9'
        fobj = @F9;
        lb=-5.12;
        ub=5.12;
        dim=100;
        
    case 'F10'
        fobj = @F10;
        lb=-32;
        ub=32;
        dim=100;
        
    case 'F11'
        fobj = @F11;
        lb=-600;
        ub=600;
        dim=100;
        
    case 'F12'
        fobj = @F12;
        lb=-50;
        ub=50;
        dim=100;
        
    case 'F13'
        fobj = @F13;
        lb=-50;
        ub=50;
        dim=100;
end

% F1

function o = F1(x)
so = -30;
x = x - so*ones(1,size(x,2));
o=sum(x.^2);
end

% F2

function o = F2(x)
so = -3;
x = x - so*ones(1,size(x,2));
o=sum(abs(x))+prod(abs(x));
end

% F3

function o = F3(x)
so = -30;
x = x - so*ones(1,size(x,2));    
dim = size(x,2);
o=0;
for i=1:dim
    o=o+sum(x(1:i))^2;
end
end

% F4

function o = F4(x)
so = -30;
x = x - so*ones(1,size(x,2));
o=max(abs(x));
end

% F5

function o = F5(x)
so = -15;
x = x - so*ones(1,size(x,2));
dim=size(x,2);
o=sum(100*(x(2:dim)-(x(1:dim-1).^2)).^2+(x(1:dim-1)-1).^2);
end

% F6

function o = F6(x)
so = -750;
x = x - so*ones(1,size(x,2));
o=sum(floor((x+.5)).^2);
end

% F7

function o = F7(x)
so = -0.25;
x = x - so*ones(1,size(x,2));
dim=size(x,2);
o=sum([1:dim].*(x.^4))+rand;
end

% F8

function o = F8(x)
so = -300;
x = x - so*ones(1,size(x,2));
o=sum(-x.*sin(sqrt(abs(x))));
end

% F9

function o = F9(x)
so = -2;
x = x - so*ones(1,size(x,2));
dim=size(x,2);
o=sum(x.^2-10*cos(2*pi.*x))+10*dim;
end

% F10

function o = F10(x)
so = -16;
x = x - so*ones(1,size(x,2));
dim=size(x,2);
o=-20*exp(-.2*sqrt(sum(x.^2)/dim))-exp(sum(cos(2*pi.*x))/dim)+20+exp(1);
end

% F11

function o = F11(x)
so = -400;
dim=size(x,2);
x = x - so*ones(1,dim);
o=sum(x.^2)/4000-prod(cos(x./sqrt([1:dim])))+1;
end

% F12

function o = F12(x)
so = -30;
x = x - so*ones(1,size(x,2));
dim=size(x,2);
o=(pi/dim)*(10*((sin(pi*(1+(x(1)+1)/4)))^2)+sum((((x(1:dim-1)+1)./4).^2).*...
(1+10.*((sin(pi.*(1+(x(2:dim)+1)./4)))).^2))+((x(dim)+1)/4)^2)+sum(Ufun(x,10,100,4));
end

% F13

function o = F13(x)
so = -100;
x = x - so*ones(1,size(x,2));
dim=size(x,2);
o=.1*((sin(3*pi*x(1)))^2+sum((x(1:dim-1)-1).^2.*(1+(sin(3.*pi.*x(2:dim))).^2))+...
((x(dim)-1)^2)*(1+(sin(2*pi*x(dim)))^2))+sum(Ufun(x,5,100,4));
end
end