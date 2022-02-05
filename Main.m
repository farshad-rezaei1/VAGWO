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

% The initial parameters that you need are:
%_________________________________________________________________________
% fobj=@YourCostFunction
% nx=number of your variables
% lb=the lower bound of variables which can generally be a fixed number or a vector
% ub=the upper bound of variables which can generally be a fixed number or a vector
% notice: if the lower nad upper bounds are not fixed for all variables, 
% they appear in the forms of the vectors "varmin" and "varmax", as illustrated in following

% To run VAGWO: [z_iter,z_final,pos_final]=VAGWO(np,nx,maxit,varmax,varmin,velmax,velmin,k_max,k_min,a_max,a_min,c_max,c_min,fobj,elitism);

%_________________________________________________________________________
% Set the required parameters to run the VAGWO algorithm

% This code is for solving the minimization problems. To maximize a desired 
% cost function,please implement this code upon inverting the sign of the cost function

clc
clear
close all
tic
run=1; % Maximum number of the algorithm runnings conducted
np=30; % Number of search agents
Function_name='F1'; % Name of the test function that can be from F1 to F13 
maxit=1000; % Maximum number of iterations
elitism=1; % elitism=1: Elitism is applied; elitism=any other number: Elitism is not applied 
a_max=sqrt(2); % Upper bound of the acceleration coefficient
a_min=0; % Lower bound of the acceleration coefficient
c_max=1; % Upper bound of the leading wolves multipliers
c_min=0; % Lower bound of the leading wolves multipliers
k_max=0.9; % Upper bound of the inertia weight
k_min=0.4; % Lower bound of the inertia weight
[lb,ub,nx,fobj]=Objective_Function(Function_name); % Load details of the selected benchmark function
varmax=ub*ones(1,nx); % Upper bound defined for the positions which can generally be a desired vector
varmin=lb*ones(1,nx); % Lower bound defined for the positions which can generally be a desired vector
limvel=0.1; % A ratio of the maximum distance in the search space to form the maximum velocity 
velmax=limvel*(varmax(1,1:nx)-varmin(1,1:nx)); % Upper bound defined for the velocities
velmin=-velmax; % Lower bound defined for the velocities
z_iter_main=zeros(run,maxit);
z_final_main=zeros(run);
pos_final_main=zeros(run,nx);
x1=zeros(maxit);
y1=zeros(maxit);

% Run the VAGWO algorithm for "run" times 
for nrun=1:run
    [z_iter,z_final,pos_final]=VAGWO(np,nx,maxit,varmax,varmin,velmax,velmin,k_max,k_min,a_max,a_min,c_max,c_min,fobj,elitism);
     z_iter_main(nrun,1:maxit)=z_iter(1:maxit);
     z_final_main(nrun)=z_final;
     pos_final_main(nrun,1:nx)=pos_final(1:nx);
end

% Display the comprehensive results
disp(['The final statistical results calculated when implementing the VAGWO algorithm for ',num2str(run),' times are as follows:']);
disp(['The average of the final objective function values calculated over ',num2str(run),' times = ',num2str(mean(z_final_main(1:run)))]);
disp(['The median of the final objective function values calculated over ',num2str(run),' times = ',num2str(median(z_final_main(1:run)))]);
disp(['The best of the final objective function values calculated over ',num2str(run),' times = ',num2str(min(z_final_main(1:run)))]);
disp(['The standard deviation of the final objective function values calculated over ',num2str(run),' times = ',num2str(std(z_final_main(1:run)))]);

% Plot the convergence curve of the VAGWO over the course of iterations
for i=1:maxit
    x1(i)=i;sum1=0;
    for j=1:run
        sum1=sum1+z_iter_main(j,i);
    end
    y1(i)=sum1/run;
end
semilogy(x1,y1,'-r')
xlabel('Iteration');
ylabel('Average best-so-far');
legend('VAGWO');
hold on
time_vagwo = toc;
disp(['Elapsed time of running the VAGWO for ',num2str(run),' times = ',num2str(time_vagwo),' seconds']);