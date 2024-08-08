close all
clc
clearvars


% We consider a mixed-autonomy vehicle chain, 
%   with one head vehicle, one CAV, and n following HVs

n = 2; % number of following vehicles

%% Settings of safepolicy and CBF

safepolicy.name = 'SDH';
% We provide three safe policies: 
%  constant time headway (CTH)
%  time to collision (TTC)
%  stoppinng distance headway (SDH)

safepolicy.tau = [0.8,0.8,0.8,0.8,0.8]; 
% the safe time headway of each vehicle: [CAV, HV1, HV2, HV3, HV4]


CBF.K = [10,10,10,10,10]; 
% The extended K_infinity function $alpha(x)$ in CBF 
% We use linear function  a(x) = Kx, 
% But Other functions can also be used

CBF.penalty = [100,100,100,100];
% penalty in the CBF-QP for each following vehicles

CBF.eta = [0.6,0.6,0.6,0.6];


%% run simulations

% We use two safety-critical traffic scenario examples: 
%   head vehicle suddenely decelerates: ACCIDENT.type = 1;
%   follow vehicle suddenly decelerates: ACCIDENT.type = 2;

ACCIDENT.type = 1;
ACCIDENT.begin = 5; % the time when accident happens
ACCIDENT.duration = 4; % duration of the accident deceleration
ACCIDENT.acc = 5; % acceleration


useCBF = 0;
checkissafe(n,ACCIDENT,useCBF,CBF,safepolicy);
useCBF = 1;
checkissafe(n,ACCIDENT,useCBF,CBF,safepolicy);



ACCIDENT.type = 2;
ACCIDENT.begin = 5; % the time when accident happens
ACCIDENT.duration = 3; % duration of the accident deceleration
ACCIDENT.acc = 5; % acceleration
useCBF = 0;
checkissafe(n,ACCIDENT,useCBF,CBF,safepolicy);
useCBF = 1;
checkissafe(n,ACCIDENT,useCBF,CBF,safepolicy);



function [] = checkissafe(n,ACCIDENT,useCBF,CBF,safepolicy)
SimTime = 40;
Tstep = 0.05;
NumStep = SimTime/Tstep;

switch ACCIDENT.type
    case 1
        savename = 'HeadDec_';
    case 2
        savename = 'FollowAcc_';
end

if useCBF
   savename = [savename 'CBF'];
else
   savename = [savename 'nominal'];
end

disp(savename);

if ~exist('result/','dir')
    mkdir ('result/')
end



v_star = 20;   % Equilibrium speed

acel_max = 7; % maximum accelearion
dcel_max = -7; % maximum deceleration


s_star = HV_sstar(v_star);
s_star_CAV = s_star;

tic


S = zeros(NumStep,n+2,3);
% S stores state of each vehicle
% It is a three-dimension matrix
% S(k,i,j) k:time-step; i:vehicle index; j=x/v/a

% vehicle index
% 1: head vehicle
% 2: CAV
% 3-(n+2): following vehicles


S(1,:,2) = v_star * ones(n+2,1); 
S(1,1,1) = 0;
S(1,2,1) = -s_star_CAV;
for vi=1:n
    S(1,vi+2,1) = S(1,vi+1,1) - s_star;
end

sigma_record = zeros(NumStep-1,n);

%% Simulation

u = zeros(NumStep,1);

for k = 1:NumStep - 1
    % following HV model
    acel = HV_model(S(k,1:(end-1),1) - S(k,2:end,1),S(k,:,2));
    acel(acel>acel_max) = acel_max;
    acel(acel<dcel_max) = dcel_max;
    
    
    S(k,2:end,3) = acel;
    S(k,1,3) = 0; % the preceding vehicle
    
    switch ACCIDENT.type
        case 1  % head vehicle suddenly brakes
            if (k*Tstep>ACCIDENT.begin)&&(k*Tstep<=ACCIDENT.begin+ACCIDENT.duration)
                 S(k,1,3) = -ACCIDENT.acc;
            end
            if (k*Tstep>ACCIDENT.begin+ACCIDENT.duration)&&(k*Tstep<=ACCIDENT.begin+2*ACCIDENT.duration)
                 S(k,1,3) = ACCIDENT.acc;
            end
        case 2  % following vehicles suddenly accelerates
            if (k*Tstep>ACCIDENT.begin)&&(k*Tstep<=ACCIDENT.begin+ACCIDENT.duration)
                S(k,2+n,3) = ACCIDENT.acc; 
            end
    end
    
    
    
    Sk = squeeze(S(k,:,:));
    % nominal controller, 
    % The designed STC has a flexibility in choosing this nominal controller,
    %  any pre-designed controllers can be used here
    u(k) = nominal_controller(Sk,s_star,v_star,s_star_CAV);

    if useCBF
        if strcmp(safepolicy.name,'CTH')
            [u(k),sigma_k] = CBF_CTH(Sk,u(k),CBF,safepolicy);
        end
        if strcmp(safepolicy.name,'TTC')
            [u(k),sigma_k] = CBF_TTC(Sk,u(k),CBF,safepolicy);
        end
        if strcmp(safepolicy.name,'SDH')
            [u(k),sigma_k] = CBF_SDH(Sk,u(k),CBF,safepolicy);
        end
        sigma_record(k,:) = sigma_k;
    end

    u(k) = max(min(u(k),acel_max), dcel_max);
    

    S(k,2,3) = u(k);
    
    S(k+1,:,2) = S(k,:,2) + Tstep*S(k,:,3);
    S(k+1,:,1) = S(k,:,1) + Tstep*S(k,:,2);
end


%% Save

save(['result/' savename '.mat'],'S');

tsim = toc;

fprintf('  ends at %6.4f seconds \n', tsim);
plot_trajectory(savename,Tstep,safepolicy)
end






function [u0] = nominal_controller(S,s_star,v_star,s_star_CAV)
    K = [1.26,-1.5,-2,0.2,-2,0.2,0.9];
    
    [n,~] = size(S);
    n = n-2;
    X = zeros(2*n+3,1);
    %X = [s0,v0, s1,v1, ..., sn,vn, vhead]

    X(1) = S(1,1) - S(2,1) - s_star_CAV;
    X(3:2:2*n+2) = reshape(S(2:n+1,1)-S(3:n+2,1), n, 1) - s_star;
    X(2:2:2*n+2) = reshape(S(2:n+2,2), n+1, 1) - v_star;
    X(2*n+3) = S(1,2)-v_star;
    u0 = K*X;

end


function [h] = safe_index(s,v,v_l,policy,tau)
    switch policy
        case 'CTH'
            h = s - tau*v;
        case 'TTC'
            h = s - (v-v_l)*tau;
        case 'SDH'
            a_max = 7;
            dv = v-v_l;
            h = s - tau*dv-dv.*dv/(2*a_max);
    end
end


function [s_star] = HV_sstar(v_star)
    s_st  = 5;
    s_go  = 35;
    v_max = 40;
    
    s_star   = acos(1-v_star/v_max*2)/pi*(s_go-s_st)+s_st; % Equilibrium spacing
end

function [acc] = HV_model(s,v)
    alpha = 0.6;
    beta  = 0.9;
    s_st  = 5;
    s_go  = 35;
    v_max  = 40;
    
    cal_D = min(s_go,max(s_st,s));
    acc = alpha*(v_max/2*(1-cos(pi*(cal_D-s_st)/(s_go-s_st))) -  v(2:end))+beta*(v(1:(end-1)) - v(2:end));
end


function [u,sigma] = CBF_SDH(S,u0,CBF,policy)
    [n,~] = size(S);
    
    n = n-2;
    A = zeros(n+1,n+1);
    b = zeros(n+1,1);

    a_max = 7;

    
    dv = S(2,2) - S(1,2);

    Lfh0 = S(1,2)-S(2,2) + policy.tau(1)*S(1,3)+dv/a_max*S(1,3);
    Lgh0 = -policy.tau(1) - dv/a_max;
    h0 = safe_index(S(1,1)-S(2,1),S(2,2),S(1,2),policy.name,policy.tau(1));
    A(1,1) = -Lgh0;
    b(1) = Lfh0 + CBF.K(1)*h0;


    dv = S(3,2) - S(2,2);
    Lfh1 = S(2,2) - S(3,2) - policy.tau(2)*S(3,3)-dv/a_max*S(3,3);
    Lgh1 = policy.tau(2) + dv/a_max;
    h1 = safe_index(S(2,1)-S(3,1),S(3,2),S(2,2), policy.name,policy.tau(2));
    Lf_hbar = Lfh1 - CBF.eta(1)*Lfh0;
    Lg_hbar = Lgh1 -CBF.eta(1)*Lgh0;
    hbar = h1 - h0;
    A(2,1) = -Lg_hbar;
    A(2,2) = -1;
    b(2) = Lf_hbar + CBF.K(2)*hbar;

    for i=2:n
        dv = S(i+2,2) - S(i+1,2);
        hi = safe_index(S(i+1,1)-S(i+2,1),S(i+2,2),S(i+1,2), policy.name,policy.tau(i+1));
        Lfhi = S(i+1,2) - S(i+2,2) - policy.tau(i+1)*(S(i+2,3)-S(i+1,3))-dv/a_max*(S(i+2,3)-S(i+1,3));
        Lghi = 0;
        Lf_hbar = Lfhi - CBF.eta(i)*Lfh0;
        Lg_hbar = Lghi -CBF.eta(i)*Lgh0;
        hbar = hi - h0;
        A(i+1,1) = -Lg_hbar;
        A(i+1,i+1) = -1;
        b(i+1) = Lf_hbar + CBF.K(i+1)*hbar;
    end

    H=eye(n+1);
    H(1,1)=1;
    for i=2:n+1
        H(i,i) = CBF.penalty(i-1);
    end

    f = zeros(n+1,1);
    f(1) = -u0;

    options =  optimset('Display','off');
    u_slacks = quadprog(H, f, A, b, [], [], [], [], [], options);
    u=u_slacks(1);
    sigma = u_slacks(2:end);
end



function [u,sigma] = CBF_TTC(S,u0,CBF,policy)
    [n,~] = size(S);
    
    n = n-2;
    A = zeros(n+1,n+1);
    b = zeros(n+1,1);

    Lfh0 = S(1,2) - S(2,2) + policy.tau(1)*S(1,3);
    Lgh0 = -policy.tau(1);
    h0 = safe_index(S(1,1)-S(2,1),S(2,2),S(1,2),policy.name,policy.tau(1));
    A(1,1) = -Lgh0;
    b(1) = Lfh0 + CBF.K(1)*h0;

    Lfh1 = S(2,2) - S(3,2) - policy.tau(2)*S(3,3);
    Lgh1 = policy.tau(2);
    h1 = safe_index(S(2,1)-S(3,1),S(3,2),S(2,2), policy.name,policy.tau(2));
    Lf_hbar = Lfh1 - CBF.eta(1)*Lfh0;
    Lg_hbar = Lgh1 -CBF.eta(1)*Lgh0;
    hbar = h1 - h0;
    A(2,1) = -Lg_hbar;
    A(2,2) = -1;
    b(2) = Lf_hbar + CBF.K(2)*hbar;

    for i=2:n
        hi = safe_index(S(i+1,1)-S(i+2,1),S(i+2,2),S(i+1,2), policy.name,policy.tau(i+1));
        Lfhi = S(i+1,2) - S(i+2,2) - policy.tau(i+1)*(S(i+2,3)-S(i+1,3));
        Lghi = 0;
        Lf_hbar = Lfhi - CBF.eta(i)*Lfh0;
        Lg_hbar = Lghi -CBF.eta(i)*Lgh0;
        hbar = hi - h0;
        A(i+1,1) = -Lg_hbar;
        A(i+1,i+1) = -1;
        b(i+1) = Lf_hbar + CBF.K(i+1)*hbar;
    end

    H=eye(n+1);
    H(1,1)=1;
    for i=2:n+1
        H(i,i) = CBF.penalty(i-1);
    end

    f = zeros(n+1,1);
    f(1) = -u0;

    options =  optimset('Display','off');
    u_slacks = quadprog(H, f, A, b, [], [], [], [], [], options);
    u=u_slacks(1);
    sigma = u_slacks(2:end);
end



function [u,sigma] = CBF_CTH(S,u0,CBF,policy)
    [n,~] = size(S);
    
    n = n-2;
    A = zeros(n+1,n+1);
    b = zeros(n+1,1);

    Lfh0 = S(1,2) - S(2,2);
    Lgh0 = -policy.tau(1);
    h0 = safe_index(S(1,1)-S(2,1),S(2,2),S(1,2),policy.name,policy.tau(1));
    
    A(1,1) = -Lgh0;
    b(1) = Lfh0 + CBF.K(1)*h0;

    for i=1:n
        hi = safe_index(S(i+1,1)-S(i+2,1),S(i+2,2),S(i+1,2), policy.name,policy.tau(i+1));
        Lfhi = S(i+1,2) - S(i+2,2) - policy.tau(i+1)*S(i+2,3);
        Lghi = 0;
        Lf_hbar = Lfhi - CBF.eta(i)*Lfh0;
        Lg_hbar = Lghi -CBF.eta(i)*Lgh0;
        hbar = hi - h0;
        A(i+1,1) = -Lg_hbar;
        A(i+1,i+1) = -1;
        b(i+1) = Lf_hbar + CBF.K(i+1)*hbar;
    end

    H=eye(n+1);
    H(1,1)=1;
    for i=2:n+1
        H(i,i) = CBF.penalty(i-1);
    end

    f = zeros(n+1,1);
    f(1) = -u0;

    options =  optimset('Display','off');
    u_slacks = quadprog(H, f, A, b, [], [], [], [], [], options);
    u=u_slacks(1);
    sigma = u_slacks(2:end);
end


function [] = plot_trajectory(resultfile,Tstep,policy)
    load(['result/' resultfile '.mat'],"S");
    
    [Num,n,~] = size(S);
    n = n-2;
    t = 1:Num;
    t = t*Tstep;



    figure()
    set(gcf,"unit","centimeters","Position",[2,2,30,16])
    

    subplot(1,4,1)
    hold on
    plot(t,S(:,1,1) - S(:,2,1))
    for i=1:n
        plot(t,S(:,i+1,1) - S(:,i+2,1))
    end
    legend('CAV','HV1','HV2')
    xlabel('time (s)')
    ylabel('gap (m)')


    subplot(1,4,2)
    hold on
    
    plot(t,S(:,2,2))
    for i=1:n
        plot(t,S(:,i+2,2))
    end
    plot(t,S(:,1,2))
    legend('CAV','HV1','HV2','Head')
    xlabel('time (s)')
    ylabel('speed (m/s)')
    
    
    subplot(1,4,3)
    hold on

    plot(t,S(:,2,3))
    for i=1:n
        plot(t,S(:,i+2,3))
    end
    plot(t,S(:,1,3))
    legend('CAV','HV1','HV2','Head')
    xlabel('time (s)')
    ylabel('acc (m/s^2)')
    
    subplot(1,4,4)
    hold on
    plot(t,safe_index(S(:,1,1)-S(:,2,1),S(:,2,2),S(:,1,2),policy.name,policy.tau(1)))
    for i=1:n
        plot(t,safe_index(S(:,i+1,1)-S(:,i+2,1),S(:,i+2,2),S(:,i+1,2),policy.name,policy.tau(i+1)))
    end
    legend('CAV','HV1','HV2')
    xlabel('time (s)')
    ylabel('h (m)')

    sgtitle(resultfile,Interpreter='none')
end
