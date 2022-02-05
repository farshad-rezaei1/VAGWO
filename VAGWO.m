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

% VAGWO algorithm                                                                  
function [z_iter,z_final,pos_final] = VAGWO(np,nx,maxit,varmax,varmin,velmax,velmin,k_max,k_min,a_max,a_min,c_max,c_min,fobj,elitism)
it=1;
% disp(['Number of Iterations = ',num2str(it)]);
pp_pbest=zeros(np,nx);
pv=zeros(np,nx);
optimal_pos=zeros(1,nx);
z=zeros(np);
z_pbest=zeros(np);
pos_final=zeros(nx);
z_iter=zeros(maxit);
z_alpha=inf;
z_beta=inf;
z_delta=inf;

% Initialization process of the algorithm
[pp,gv_alpha,gv_beta,gv_delta]=Initialization(np,nx,varmax,varmin,velmax,velmin);

% Start the optimization process

% Objective function evaluations and determine the personal best solutions
% and objectives, if elitism is applied
for j=1:np
    z(j)=fobj(pp(j,1:nx));
    
    % Elitism
    if elitism==1
        z_pbest(j)=z(j);
        pp_pbest(j,1:nx)=pp(j,1:nx);
    end
end
for j=1:np
    if z(j)<=z_alpha
        z_alpha=z(j);
        alpha(1,1:nx)=pp(j,1:nx);
    elseif z(j)>z_alpha && z(j)<=z_beta
        z_beta=z(j);
        beta(1,1:nx)=pp(j,1:nx);
    elseif z(j)>z_beta && z(j)<=z_delta
        z_delta=z(j);
        delta(1,1:nx)=pp(j,1:nx);
    end
end
z_optimal(it)=z_alpha;
optimal_pos(it,:)=alpha(1,:);

% Save the best-so-far objective value in the current run
z_iter(it)=z_optimal(it);

% The Main Loop
while it<maxit
    it=it + 1;
    aa=a_max-(a_max-a_min)*(it-1)/(maxit-1); % Eq.(16)
    cc=c_max-(c_max-c_min)*(it-1)/(maxit-1); % Eq.(23)
    k=k_max-(k_max-k_min)*(it-1)/(maxit-1); % Eq.(24)
%     disp(['Number of Iterations= ',num2str(it)]);
    for j=1:np        
        a_alpha(1,1:nx)=(2*rand(1,nx)-ones(1,nx))*(aa^2); % Eq.(13)
        c_alpha(1,1:nx)=ones(1,nx)+(2*rand(1,nx)-ones(1,nx))*(cc^2); % Eq.(20)
        gv_alpha(j,1:nx)=k*(sign(a_alpha(1,1:nx)).*abs(gv_alpha(j,1:nx)))+...
            a_alpha(1,1:nx).*abs(c_alpha(1,1:nx).*alpha(1,1:nx)-pp(j,1:nx)); % Eq.(10)
        
        a_beta(1,1:nx)=(2*rand(1,nx)-ones(1,nx))*(aa^2); % Eq.(14)
        c_beta(1,1:nx)=ones(1,nx)+(2*rand(1,nx)-ones(1,nx))*(cc^2); % Eq.(21)
        gv_beta(j,1:nx)=k*(sign(a_beta(1,1:nx)).*abs(gv_beta(j,1:nx)))+...
            a_beta(1,1:nx).*abs(c_beta(1,1:nx).*beta(1,1:nx)-pp(j,1:nx)); % Eq.(11)
        
        a_delta(1,1:nx)=(2*rand(1,nx)-ones(1,nx))*(aa^2); % Eq.(15)
        c_delta(1,1:nx)=ones(1,nx)+(2*rand(1,nx)-ones(1,nx))*(cc^2); % Eq.(22);
        gv_delta(j,1:nx)=k*(sign(a_delta(1,1:nx)).*abs(gv_delta(j,1:nx)))+...
            a_delta(1,1:nx).*abs(c_delta(1,1:nx).*delta(1,1:nx)-pp(j,1:nx)); % Eq.(12)
        
        sum1=gv_alpha(j,:)+gv_beta(j,:)+gv_delta(j,:);
        sum2=alpha+beta+delta;
        pv(j,1:nx)=sum1/3;
        
        % Return back the velocity of the particles if going beyond the velocity boundaries
        flag4lbv=pv(j,:)<velmin(1,:);
        flag4ubv=pv(j,:)>velmax(1,:);
        pv(j,:)=(pv(j,:)).*(~(flag4lbv+flag4ubv))+velmin.*flag4lbv+velmax.*flag4ubv;
        
        pp(j,1:nx)=sum2/3-pv(j,:); % Eq.(28)
            
        % Return back the position of the particles if going beyond the position boundaries
        flag4lbp=pp(j,:)<varmin(1,:);
        flag4ubp=pp(j,:)>varmax(1,:);
        pp(j,:)=(pp(j,:)).*(~(flag4lbp+flag4ubp))+varmin.*flag4lbp+varmax.*flag4ubp; 
    
    % Objective function evaluations and determining of the personal best solutions and objectives
        z(j)=fobj(pp(j,:));
    end
    if elitism==1
        for j=1:np
            if z_pbest(j)<z(j)
                z(j)=z_pbest(j);
                pp(j,:)=pp_pbest(j,:);
            else
                z_pbest(j)=z(j);
                pp_pbest(j,:)=pp(j,:);
            end
        end
    end
    for j=1:np
        if z(j)<=z_alpha
            z_alpha=z(j);
            alpha(1,1:nx)=pp(j,1:nx);
        elseif z(j)>z_alpha && z(j)<=z_beta
            z_beta=z(j);
            beta(1,1:nx)=pp(j,1:nx);
        elseif z(j)>z_beta && z(j)<=z_delta
            z_delta=z(j);
            delta(1,1:nx)=pp(j,1:nx);
        end
    end
    z_optimal(it)=z_alpha;
    optimal_pos(it,:)=alpha(1,:);
    
    % Save the best-so-far objective value in the current run
    z_iter(it)=z_optimal(it);
end

% Save the final best solution and objective revealed upon the end of the optimization process
z_final=z_optimal(maxit);
pos_final(1:nx)=optimal_pos(maxit,1:nx);
end