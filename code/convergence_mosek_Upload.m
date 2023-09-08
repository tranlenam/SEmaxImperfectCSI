%%%% robust design WSR
%%%% Author Quang-Doanh Vu
%%%% Note: Need Mosek installed to run
%%
clear all
clc
%% system
B = 1; % number of BS
U = 4; % number of users per BS
M = 20; % number of antenna at a BS
%%%%
%% parameter
noisevar = 1; %(ALWAYS KEEP VALUE 1, IF CHANGE LOOK OTHER SCHEMES FOR SURE)
Pow = 100; % total transmit power at a BS (20 dB)
w = ones(U*B,1); %% weighted rate
%% Scaling factor: to scale value of f function based on number of uer and power
%%%%%% use for proposed stochastic model
scale_factor = (1/Pow)*(U/6);
noisevar_scale = noisevar*scale_factor;
%% error parameter (ellipsoilds)
%%%%%%%% Gaussian estimation error
sigmaerror = 0.01;%%%% sigmaerror
errorprob = 0.1; %%% error probability (90%)(fix)
d_ellipsoid = sqrt(sigmaerror*chi2inv(1-errorprob,2*M)/2);%%inverse CDF Chi square
Imax=5000;
%%%%%%%% ellipsoild
paraep = 1;
E = repmat(paraep*eye(M),[1 1 B U*B]);
Eh = zeros(M,M,B,U*B);%%
for iu = 1:U*B
    for ik = 1:B
        Eh(:,:,ik,iu)= E(:,:,ik,iu)^(-0.5);
    end
end
delta = (d_ellipsoid^2)*ones(B,U*B);
%% (random) channel (estimate)
% chan = sqrt((1-sigmaerror)/2)*(randn(1, M, B, U*B)+1i*randn(1, M, B, U*B));
%% real (random) channel (for test feasibility)
%chan_real = chan + sqrt(sigmaerror/2)*(randn(1, M, B, U*B)+1i*randn(1, M, B, U*B));
%%
load ../data/channelconvergeU4_revised.mat
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% setting parameters for stochastic

%%%%%%%%%%%%% for general
kappa = 0.85;%%% keep this for all
tau1 = 0.9; tau2 = 0.91; %% parameter for stepsize
tausquad = 10^-4; %%%% parameter of quadratic term (rho in paper), same for all constraint
bigM = 5000;
gamma_pa_0=1;
mu_ini = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initial
%       mu_old = mu_ini*ones(U,1);
%       F_old = zeros(2*U*M+3,U);
%       v_old_temp = randn(2*M,U);
%       v_old = sqrt(0.9*Pow)*v_old_temp/norm(myvec(v_old_temp));
load ../data/initial_point_U4.mat

%%%%%%%%%%%%%%%%%
stopflag =0;
loop = 0;
objsequence = 0;
timesequence = 0;
while(stopflag ==0)
    time1 = 0;time2=0;

    rho_pa = (1/(1+loop))^tau1; %%%
    gamma_pa= gamma_pa_0*(1/(1+loop))^tau2;
    %%%% new error sample
    chan_sample = chan + sqrt(sigmaerror/2)*(randn(1, M, B, U*B)+1i*randn(1, M, B, U*B));
    %%%%% variables
    clear prob;
    [r, res] = mosekopt('echo(0) symbcon');
    lengthva = 2*M*U + 4*U + U-2 + U-2 + 4;

    lengthva2 = 2*M*U + 4*U  + 4;

    %%%% linear constraint
    acon = [];
    lowcon = [];
    upcon= [];

    acon2 = [];
    lowcon2 = [];
    upcon2 = [];

    lowboundx = [-inf*ones(1,2*M*U),zeros(1,U),ones(1,U),zeros(1,U),ones(1,U),zeros(1,U-2),zeros(1,U-2),1,0,0,0];
    upboundx =  [inf*ones(1,2*M*U),bigM*ones(1,U),bigM*ones(1,U)+1,inf*ones(1,U),ones(1,U),inf*ones(1,U-2),inf*ones(1,U-2),1,Pow,inf,inf];

    lowboundx2 = [-inf*ones(1,2*M*U),zeros(1,U),ones(1,U),zeros(1,U),ones(1,U),1,0,0,0];
    upboundx2 =  [inf*ones(1,2*M*U),bigM*ones(1,U),bigM*ones(1,U)+1,inf*ones(1,U),ones(1,U),1,Pow,inf,inf];

    %%% quadratic constraint
    contype = [];
    cones_sub = [];
    cones_subptr = [];
    cone_par = [];

    contype2 = [];
    cones_sub2 = [];
    cones_subptr2 = [];
    cone_par2 = [];
    %%%%

    %%%%% geomean function by power
    contype = [contype,res.symbcon.MSK_CT_PPOW];
    cone_par = [cone_par,0.5];
    cones_subptr = [cones_subptr,length(cones_sub)+1];
    cones_sub = [cones_sub,2*M*U+U+1,2*M*U+U+2,2*M*U+4*U+1];

    temp0 = zeros(1,lengthva);temp0(2*M*U+4*U+1)=1;temp0(2*M*U+4*U+1+U-2)=-1;
    acon = [acon;temp0]; lowcon = [lowcon,0]; upcon = [upcon,0];

    for ipoint = 4:U
        contype = [contype,res.symbcon.MSK_CT_PPOW];
        cone_par = [cone_par,1-1/(ipoint-1)];
        cones_subptr = [cones_subptr,length(cones_sub)+1];
        cones_sub = [cones_sub,2*M*U+4*U+ipoint-3+U-2,2*M*U+U+ipoint-1,2*M*U+4*U+ipoint-2];

        temp0 = zeros(1,lengthva);temp0(2*M*U+4*U+ipoint-2)=1;temp0(2*M*U+4*U+ipoint-2+U-2)=-1;
        acon = [acon;temp0]; lowcon = [lowcon,0]; upcon = [upcon,0]; %%% 1+mu=mutide
    end
    contype = [contype,res.symbcon.MSK_CT_PPOW];
    cone_par = [cone_par,1-1/U];
    cones_subptr = [cones_subptr,length(cones_sub)+1];
    cones_sub = [cones_sub,2*M*U+4*U+U-2+U-2,2*M*U+U+U,lengthva];
    %                %%%%%% end of geomean

    for iuser = 1:U
        Hchan = chan_sample(:,:,1,iuser)'*chan_sample(:,:,1,iuser);
        Htide = scale_factor*[real(Hchan),-imag(Hchan);imag(Hchan),real(Hchan)];

        valx = v_old(:,iuser)'*Htide*v_old(:,iuser)-...
            mu_old(iuser)*(myvec((v_old(:,[1:iuser-1,iuser+1:U])'*Htide )')'*myvec(v_old(:,[1:iuser-1,iuser+1:U]))+noisevar_scale);
        g_0 = 0.5+0.5*tanh(kappa*valx);
        ghat = 0.5*kappa*(1-(tanh(kappa*valx))^2);

        grad_g = zeros(2*U*M+1,1); %%% there is U constraints

        for juser = 1:U
            if (juser~=iuser)
                grad_temp = -2*ghat*mu_old(iuser)*Htide*v_old(:,juser);
            else
                grad_temp = 2*ghat*Htide*v_old(:,juser);
            end
            grad_g(1+2*(juser-1)*M:2*juser*M) = grad_temp;
        end
        grad_g(2*U*M+1)= -ghat*(myvec((v_old(:,[1:iuser-1,iuser+1:U])'*Htide )')'*myvec(v_old(:,[1:iuser-1,iuser+1:U]))+noisevar_scale);%%gradient for mu variable
        %%%%%%%%%% construct surrogate function for constraint iuser
        %%% gtide
        newconst = g_0-grad_g'*[myvec(v_old);mu_old(iuser)]-...
            tausquad*sum(([myvec(v_old);mu_old(iuser)]).^2);
        para_order1 = grad_g+2*tausquad*[myvec(v_old);mu_old(iuser)];
        para_order2 = -tausquad;

        const =   (1-rho_pa)*F_old(2*U*M+2,iuser) +  rho_pa*newconst;
        const_order1 = (1-rho_pa)*F_old(1:2*U*M+1,iuser)+ rho_pa* para_order1;
        const_order2 = (1-rho_pa)*F_old(end,iuser)+ rho_pa* para_order2;


        %%% update F_old for recursive
        F_old(2*U*M+2,iuser) = const;F_old(1:2*U*M+1,iuser)=const_order1;
        F_old(end,iuser)=const_order2;
        %%%% constraint

        contype = [contype,res.symbcon.MSK_CT_PPOW];
        cones_subptr = [cones_subptr,length(cones_sub)+1];
        cones_sub = [cones_sub,2*M*U+U+U+iuser,2*M*U+U+U+U+iuser,2*M*U+iuser];
        cone_par = [cone_par,0.5];

        contype2 = [contype2,res.symbcon.MSK_CT_PPOW];
        cones_subptr2 = [cones_subptr2,length(cones_sub2)+1];
        cones_sub2 = [cones_sub2,2*M*U+U+U+iuser,2*M*U+U+U+U+iuser,2*M*U+iuser];
        cone_par2 = [cone_par2,0.5];

        %%% linear
        temp0 = zeros(1,lengthva);temp0(2*M*U+iuser)=1;temp0(2*M*U+U+iuser)=-1;
        acon = [acon;temp0]; lowcon = [lowcon,-1]; upcon = [upcon,-1];

        temp0 = zeros(1,lengthva2);temp0(2*M*U+iuser)=1;temp0(2*M*U+U+iuser)=-1;
        acon2 = [acon2;temp0]; lowcon2 = [lowcon2,-1]; upcon2 = [upcon2,-1];


        temp0 = zeros(1,lengthva);
        temp0(1:2*M*U)= const_order1(1:end-1);temp0(2*M*U+iuser)=const_order1(end);
        temp0(2*M*U+U+U+iuser)=const_order2;
        temp0(end-1)=-1;
        acon = [acon;temp0];lowcon = [lowcon,1-errorprob-const]; upcon = [upcon,inf];


        temp0 = zeros(1,lengthva2);
        temp0(1:2*M*U)= const_order1(1:end-1);temp0(2*M*U+iuser)=const_order1(end);
        temp0(2*M*U+U+U+iuser)=const_order2;
        temp0(end-1)=-1;
        temp0(end)=1;
        acon2 = [acon2;temp0];lowcon2 = [lowcon2,1-errorprob-const]; upcon2 = [upcon2,inf];

    end


    contype = [contype,res.symbcon.MSK_CT_PPOW];
    cones_subptr = [cones_subptr,length(cones_sub)+1];
    cones_sub = [cones_sub,lengthva-2,lengthva-3,[1:2*M*U]];
    cone_par = [cone_par,0.5];

    contype2 = [contype2,res.symbcon.MSK_CT_PPOW];
    cones_subptr2 = [cones_subptr2,length(cones_sub2)+1];
    cones_sub2 = [cones_sub2,lengthva2-2,lengthva2-3,[1:2*M*U]];
    cone_par2 = [cone_par2,0.5];

    %%%%additional variabels
    temp0 = zeros(1,lengthva);

    temp0(lengthva-2)=-const_order2;
    temp0(lengthva-1)=-1;
    acon = [acon;temp0]; lowcon = [lowcon,0]; upcon = [upcon,0]; %%%

    temp0 = zeros(1,lengthva2);
    temp0(lengthva2-2)=-const_order2;
    temp0(lengthva2-1)=-1;
    acon2 = [acon2;temp0]; lowcon2 = [lowcon2,0]; upcon2 = [upcon2,0]; %%%


    %%% objective
    ccon =  [zeros(1,lengthva-1),1];
    ccon2 =  [zeros(1,lengthva2-1),1];
    %%%%%%%% problem1

    prob.c   = ccon;
    prob.a   = sparse(acon);
    prob.blc = lowcon;
    prob.buc = upcon;
    prob.blx = lowboundx;
    prob.bux = upboundx;
    %%%%% cone
    prob.cones.type   = contype;
    prob.cones.sub    = cones_sub;
    prob.cones.subptr = cones_subptr;
    prob.cones.conepar = cone_par;

    tic
    [r,res]=mosekopt('maximize echo(0)',prob);
    time1 = toc;
    XX = res.sol.itr.xx;



    if strcmp(res.sol.itr.solsta,'PRIMAL_INFEASIBLE_CER') %%%% infeasible
        clear prob;
        [r, res] = mosekopt('echo(0) symbcon');



        prob.c   = ccon2;
        prob.a   = sparse(acon2);
        prob.blc = lowcon2;
        prob.buc = upcon2;
        prob.blx = lowboundx2;
        prob.bux = upboundx2;
        %%%%% cone
        prob.cones.type   = contype2;
        prob.cones.sub    = cones_sub2;
        prob.cones.subptr = cones_subptr2;
        prob.cones.conepar = cone_par2;

        tic
        [r,res]=mosekopt('minimize echo(0)',prob);
        time2 = toc;

        XX = res.sol.itr.xx;


    end
    vnew_temp = XX(1:2*M*U);vnew = reshape(vnew_temp,2*M,U);
    mu_new = XX(2*M*U+1:2*M*U+U);
    v_old = (1-gamma_pa)*v_old + gamma_pa*vnew;
    mu_old = (1-gamma_pa)*mu_old + gamma_pa*mu_new;
    loop = loop+1;
    objsequence = [objsequence,real(sum(log(1+mu_old)))];
    toltime = time1+time2;
    timepoint = timesequence(end) + toltime;
    timesequence = [timesequence, timepoint];
    if((loop>=Imax))
        stopflag = 1;
    end


end
%%%%%%%%%%% check feasibility for each user with a new (i.e. real) channel

%%%%%%%%%%%%
feasibleflag_stoch = zeros(1,U);
for iuser = 1:U
    checkval=(abs(chan_real(:,:,1,iuser)* (v_old(1:M,iuser)+1j*v_old(1+M:end,iuser))   ))^2 - mu_old(iuser)*((norm(chan_real(:,:,1,iuser)*(v_old(1:M,[1:iuser-1,iuser+1:U])+1j*v_old(1+M:end,[1:iuser-1,iuser+1:U]))))^2+noisevar);
    if(checkval<0)
        feasibleflag_stoch(iuser) = 1;
    end
end
realobj_stochastic = abs(1-feasibleflag_stoch)*log(1+mu_old);
plot(objsequence)
saveas(gcf, '../../results/ConvergencePlot.png')



