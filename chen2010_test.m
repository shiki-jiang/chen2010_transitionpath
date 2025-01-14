n_d=0; 
n_a=[151,11];
n_z=5;


N_j=66; 
Params.agejshifter=19; 
Params.J=N_j;


useTauchen=1; 


%% Parameters


% Discount factor
Params.beta=0.9622;


% Preferences
Params.sigma=2; 
Params.upsilon=0; 
Params.theta=0.8954; 


% Demographics
Params.Jr=46; 
Params.agej=1:1:Params.J; 
Params.n=0.01; 


% Deterministic productivity growth
Params.g=0.015;


% Production function
Params.alpha=0.2743; 
% Depreciation rates
Params.delta_r=0.0254; 
Params.delta_o=0.013; 
Params.delta_k=0.0951;


% Housing
Params.phi=0.05; 
Params.gamma=0.2;


% Earnings
Params.rho_z=0.96; 
Params.sigma_z_epsilon=sqrt(0.045); 
Params.sigma_z1=sqrt(0.38); 


% Pension
Params.vartheta=0.483; % Not used for anything (Chen2010 reports this a the 'replacement rate', but is not used here)


% Taxes
Params.tau_p=0.107;


Params.kappaj=[1.0000, 1.0719, 1.1438, 1.2158, 1.2842, 1.3527, 1.4212, 1.4897, 1.5582, 1.6267, 1.6952, 1.7217, 1.7438, 1.7748, 1.8014, 1.8279, 1.8545, 1.8810, 1.9075, 1.9341, 1.9606, 1.9623, 1.9640, 1.9658, 1.9675, 1.9692, 1.9709, 1.9726, 1.9743, 1.9760, 1.9777, 1.9700, 1.9623, 1.9546, 1.9469, 1.9392, 1.9315, 1.9238, 1.9161, 1.9084, 1.9007, 1.8354, 1.7701, 1.7048, 1.6396];
Params.kappaj=[Params.kappaj,zeros(1,Params.J-Params.Jr+1)]; % zeros for retirement
Params.sj=[0.9985995, 0.9985575, 0.998498, 0.9985075, 0.9986075, 0.998641, 0.998743, 0.9987235, 0.9987385, 0.998705, 0.9986825, 0.9986935, 0.9986795, 0.9986615, 0.998573, 0.9984555, 0.998305, 0.9981125, 0.9979675, 0.997858, 0.9977755, 0.997652, 0.997384, 0.99714, 0.9969315, 0.9967195, 0.99666, 0.996468, 0.9962415, 0.9959035, 0.995637, 0.9953435, 0.9950225, 0.9946615, 0.9942525, 0.9937795, 0.993245, 0.9926615, 0.992031, 0.9913415, 0.9906035, 0.9897745, 0.9887855, 0.987602, 0.9862635, 0.9847635, 0.9832185, 0.981763, 0.9804655, 0.979231, 0.977836, 0.976205, 0.974443, 0.9725355, 0.9704165, 0.967929, 0.9650565, 0.9619165, 0.9585285, 0.954784, 0.9506335, 0.945818, 0.940017, 0.933017, 0.9249375, 0.9160355];
Params.sj(end)=0; % I need this for the accidental bequests calculations. Not sure how Chen (2010) treated it.


%% Grids and exogenous shocks
% Grid for housing, from 0 to 20 (which are min and max values used by Chen2010 codes)
minh=0; % min for housing; note that h'=0 means renter
maxh=20; % maximum value of housing: 
h_grid=minh + (maxh-minh)*linspace(0,1,n_a(2))'.^3; 
minassets=-(1-Params.gamma)*maxh;
maxassets=100; % maximum value of assets (Chen 2010 codes uses 100)
n_a_neg = round(0.1*101); % use 1/10 of points in the negative assets range
asset_grid_neg = linspace(minassets,0,n_a_neg)';
asset_grid_pos=maxassets*linspace(0,1,101-n_a_neg+1)'.^3;
asset_grid = [asset_grid_neg(1:n_a_neg-1); asset_grid_pos]; % Note: both contain zero, so omit it from asset_grid_neg
Params.gamma1=0.1;
minassets1=-(1-Params.gamma1)*maxh;
asset_grid1=linspace(minassets1,99,50)';
asset_grid=unique([asset_grid;asset_grid1]);
a_grid=[asset_grid; h_grid]; % stacked column vector


if useTauchen==0
    kfttoptions.initialj1sigmaz=Params.sigma_z1;
    kfttoptions.nSigmas=sqrt(7/2); 
    [z_grid_J, pi_z_J,jequaloneDistz,otheroutputs] = discretizeLifeCycleAR1_KFTT(0,Params.rho_z,Params.sigma_z_epsilon,n_z,N_j,kfttoptions);
    z_grid_J=exp(z_grid_J);
elseif useTauchen==1
    tauchenoptions=struct();
    Tauchen_q=sqrt(7/2); % plus/minus sqrt(7/2) std deviations as the max/min grid points % This is exactly what Chen (2010) codes use
    [z_grid,pi_z]=discretizeAR1_Tauchen(0,Params.rho_z,Params.sigma_z_epsilon,n_z,Tauchen_q,tauchenoptions);
    z_grid=exp(z_grid); % Ranges from 0.7 to 1.4
    [z_mean,~,~,~]=MarkovChainMoments(z_grid,pi_z);
    z_grid=z_grid/z_mean; % renormalize grid to be mean 1 (Chen2010 codes do this)
    jequaloneDistz=MVNormal_ProbabilitiesOnGrid(log(z_grid),1,Params.sigma_z1,n_z); % Chen2010 codes do initial dist on z as log-normal
end


d_grid=[]; % no d variable


%% Return fn
DiscountFactorParamNames={'beta','sj'};


ReturnFn=@(aprime,hprime,a,h,z,kappaj,r,tau_p,theta,upsilon,sigma,gamma,phi,alpha,delta_k,delta_o,delta_r,agej,Jr,b,Tr)...
    Chen2010_ReturnFn(aprime,hprime,a,h,z,kappaj,r,tau_p,theta,upsilon,sigma,gamma,phi,alpha,delta_k,delta_o,delta_r,agej,Jr,b,Tr);


%% General eqm parameters
Params.r=0.0635; 
Params.Tr=0.05; 
Params.p=(Params.r+Params.delta_r)/(1+Params.r); 
Params.w=(1-Params.alpha)*((Params.r+Params.delta_k)/Params.alpha)^(Params.alpha/(Params.alpha-1)); 
% Because labor supply is exogenous you can actually figure out of b is without having to
% solve general eqm (point 4 in Definition 1 of his Appendix A.1).
Params.b=0.1;
% This is an initial guess, and then we set b to clear the general eqm
% below. It is just below the first creation of AggVars.




%% Solve the value function
%vfoptions.divideandconquer=1;
vfoptions=struct();%vfoptions.level1n=[9,n_a(2)];
tic;
if useTauchen==0
    [V,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid, a_grid, z_grid_J, pi_z_J, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
elseif useTauchen==1
    [V,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
end
vftime=toc




%% Age distribution
AgeWeightParamNames={'mewj'};
Params.mewj=ones(1,N_j);
for jj=1:N_j-1
    Params.mewj(jj+1)=Params.mewj(jj)*Params.sj(jj)/(1+Params.n); 
end
Params.mewj=Params.mewj/sum(Params.mewj); % Normalize total mass to one


%% Initial age 20 (period 1) distribution
jequaloneDist=zeros([n_a,n_z],'gpuArray');
[~,zeroassetindex]=min(abs(asset_grid));
jequaloneDist(zeroassetindex,1,:)=shiftdim(jequaloneDistz,-2); % initial dist of z, with zero assets and zero housing


%% Agent distribution
simoptions=struct(); % defaults
if useTauchen==0
    StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z_J,Params,simoptions);
elseif useTauchen==1
    StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
end


%% AggVars
FnsToEvaluate.A=@(aprime,hprime,a,h,z) a;
FnsToEvaluate.N=@(aprime,hprime,a,h,z,kappaj) kappaj*z;
FnsToEvaluate.H=@(aprime,hprime,a,h,z) h; % housing
FnsToEvaluate.Hr=@(aprime,hprime,a,h,z,kappaj,r,tau_p,theta,upsilon,phi,alpha,delta_k,delta_o,delta_r,agej,Jr,b,Tr) (hprime==0)*Chen2010_HousingServicesFn(aprime,hprime,a,h,z,kappaj,r,tau_p,theta,upsilon,phi,alpha,delta_k,delta_o,delta_r,agej,Jr,b,Tr); % rental housing=housing services used by renters
FnsToEvaluate.pensiontaxrevenue=@(aprime,hprime,a,h,z,tau_p,w,kappaj) tau_p*w*kappaj*z;
FnsToEvaluate.pensionspend=@(aprime,hprime,a,h,z,agej,Jr,b) (agej>=Jr)*b;
FnsToEvaluate.accidentalbeqleft=@(aprime,hprime,a,h,z,r,sj,delta_o) (1+r)*aprime*(1-sj)+(1-delta_o)*hprime*(1-sj);
FnsToEvaluate.Homeownership=@(aprime,hprime,a,h,z) (hprime>0);
if useTauchen==0
    AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist,Policy, FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_J,[],simoptions);
elseif useTauchen==1
    AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist,Policy, FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,[],simoptions);
end


Params.b=Params.b*(AggVars.pensiontaxrevenue.Mean/AggVars.pensionspend.Mean);


GEPriceParamNames={'r','Tr'};


GeneralEqmEqns.capitalmarkets=@(r,A,N,alpha,delta_k,Hr,delta_r) r-(alpha*((A-Hr*(1-((r+delta_r)/(1+r))))^(alpha-1))*(N^(1-alpha))-delta_k); % r=marginal product of capital, minus depreciation; with K'=A'-Hr'*(1-p), and p=(r+delta_r)/(1+r);
GeneralEqmEqns.accidentalbequests=@(Tr,accidentalbeqleft,n) Tr-accidentalbeqleft/(1+n); % Eqn A.2 from Appendix of Chen2010


%% Initial eqm: gamma=0.2
Params.gamma=0.2;


Params.r=0.0635;
Params.Tr=.05;
heteroagentoptions.verbose=1; 
if useTauchen==0
    p_eqm_init=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames, n_d, n_a, n_z, N_j, [], pi_z_J, d_grid, a_grid, z_grid_J, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
elseif useTauchen==1
    p_eqm_init=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames, n_d, n_a, n_z, N_j, [], pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
end


% Update Params based on the general eqm
Params.r=p_eqm_init.r;
Params.Tr=p_eqm_init.Tr;


if useTauchen==0
    [V_init,Policy_init]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid, a_grid, z_grid_J, pi_z_J, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
    StationaryDist_init=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_init,n_d,n_a,n_z,N_j,pi_z_J,Params,simoptions);
    AggVars_init=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist_init,Policy_init, FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_J,[],simoptions);
elseif useTauchen==1
    [V_init,Policy_init]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
    StationaryDist_init=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_init,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
    AggVars_init=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist_init,Policy_init, FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,[],simoptions);
end


%% Final eqm: gamma=0.1
Params.gamma=0.1;


Params.r=0.0635;
Params.Tr=.005;
if useTauchen==0
    p_eqm_final=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames, n_d, n_a, n_z, N_j, [], pi_z_J, d_grid, a_grid, z_grid_J, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
elseif useTauchen==1
    p_eqm_final=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames, n_d, n_a, n_z, N_j, [], pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
end


% Update Params based on the general eqm
Params.r=p_eqm_final.r;
Params.Tr=p_eqm_final.Tr;


if useTauchen==0
    [V_final,Policy_final]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid, a_grid, z_grid_J, pi_z_J, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
    StationaryDist_final=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_final,n_d,n_a,n_z,N_j,pi_z_J,Params,simoptions);
    AggVars_final=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist_final,Policy_final, FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_J,[],simoptions);
elseif useTauchen==1
    [V_final,Policy_final]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
    StationaryDist_final=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_final,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
    AggVars_final=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist_final,Policy_final, FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,[],simoptions);
end


save pretpath.mat



%% Quick look at initial and final before we go to transition path
fprintf(' \n')
fprintf('Initial vs final: exogenous change \n')
fprintf('Initial vs final: gamma= %1.2f vs %1.2f \n',[0.2,0.1])
fprintf('Initial vs final: general eqm params \n')
fprintf('Initial vs final: r= %1.3f vs %1.3f \n',[p_eqm_init.r, p_eqm_final.r])
fprintf('Initial vs final: Tr= %1.3f vs %1.3f \n',[p_eqm_init.Tr, p_eqm_final.Tr])
fprintf('Initial vs final: aggregates \n')
fprintf('Initial vs final: A= %2.3f vs %2.3f \n',[AggVars_init.A.Mean, AggVars_final.A.Mean])
fprintf('Initial vs final: Hr= %2.3f vs %2.3f \n',[AggVars_init.Hr.Mean, AggVars_final.Hr.Mean])
fprintf('Initial vs final: N= %2.3f vs %2.3f \n',[AggVars_init.N.Mean, AggVars_final.N.Mean])
fprintf(' \n')


max(abs(V_init(:)-V_final(:)))
max(abs(Policy_init(:)-Policy_final(:)))
max(abs(StationaryDist_init(:)-StationaryDist_final(:)))


%% Transition path
T=50;
ParamPath.gamma=0.1*ones(T,1);
PricePath0.r=[linspace(p_eqm_init.r, p_eqm_final.r, floor(T/2))'; p_eqm_final.r*ones(T-floor(T/2),1)];
PricePath0.Tr=[linspace(p_eqm_init.Tr, p_eqm_final.Tr, floor(T/2))'; p_eqm_final.Tr*ones(T-floor(T/2),1)];
TransPathGeneralEqmEqns.capitalmarkets=@(r,A,N,alpha,delta_k,Hr,delta_r) r-(alpha*((A-Hr*(1-((r+delta_r)/(1+r))))^(alpha-1))*(N^(1-alpha))-delta_k); 
TransPathGeneralEqmEqns.accidentalbequests=@(Tr,accidentalbeqleft,n) Tr-accidentalbeqleft/(1+n);
transpathoptions.weightscheme=1;
transpathoptions.verbose=1;
transpathoptions.GEnewprice=3;
transpathoptions.GEnewprice3.howtoupdate={... % a row is: GEcondn, price, add, factor
    'capitalmarkets','r',0,0.001;...  
    'accidentalbequests','Tr',0,0.001;... 
    };


transpathoptions.graphpricepath=1;
transpathoptions.graphaggvarspath=1;




PricePath=TransitionPath_Case1_FHorz(PricePath0, ParamPath, T, V_final, StationaryDist_init, jequaloneDist,n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, TransPathGeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions, simoptions, vfoptions);


save tpath.mat


PricePath=PricePath0; %DEBUG


[VPath,PolicyPath]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);


AgentDistPath=AgentDistOnTransPath_Case1_FHorz(StationaryDist_init, jequaloneDist, PricePath, ParamPath, PolicyPath, AgeWeightParamNames,n_d,n_a,n_z,N_j,pi_z, T,Params, transpathoptions, simoptions);


AggVarsPath=EvalFnOnTransPath_AggVars_Case1_FHorz(FnsToEvaluate,AgentDistPath,PolicyPath,PricePath,ParamPath, Params, T, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid,transpathoptions, simoptions);


save tpath2.mat




figure(3)
plot(0:1:T, [AggVars_init.Homeownership.Mean, AggVarsPath.Homeownership.Mean])
title('path of transtion')




%%
%load temp2.mat AgentDistPath2 AggVarsPath2
%N_a=prod(n_a);
%N_z=prod(n_z);
%AgentDistPath=reshape(AgentDistPath,[N_a,N_z,N_j,T]);
%AgentDistPath3=reshape(AgentDistPath2,[N_a,N_j,N_z,T]);
%AgentDistPath3=permute(AgentDistPath3,[1,3,2,4]);


%max(abs(AgentDistPath(:)-AgentDistPath3(:)))


% squeeze(sum(sum(sum(AgentDistPath,1),2),3))' % Should all be ones. They are
% squeeze(sum(sum(sum(AgentDistPath2,1),2),3))'


% max(max(max(abs(AgentDistPath(:,:,:,1)-reshape(StationaryDist_init,[N_a,N_z,N_j])))))
% max(max(max(abs(AgentDistPath3(:,:,:,1)-reshape(StationaryDist_init,[N_a,N_z,N_j]))))) 
% max(max(max(abs(AgentDistPath3(:,:,:,1)-permute(reshape(AgentDist_initial,[N_a,N_j,N_z]),[1,3,2])))))
% max(max(max(abs(permute(reshape(AgentDist_initial,[N_a,N_j,N_z]),[1,3,2]))-reshape(StationaryDist_init,[N_a,N_z,N_j]))))


%max(max(max(abs(AgentDistPath(:,:,:,1)-AgentDistPath3(:,:,:,1))))) % period 1 is good
%max(max(max(abs(AgentDistPath(:,:,:,2)-AgentDistPath3(:,:,:,2)))))


% squeeze(sum(sum(AgentDistPath(:,:,:,2),1),2))- squeeze(sum(sum(AgentDistPath3(:,:,:,2),1),2)) % age masses seem fine
% squeeze(sum(sum(AgentDistPath(:,:,:,2),1),3))- squeeze(sum(sum(AgentDistPath3(:,:,:,2),1),3)) % shocks are at 1e-4
% squeeze(sum(sum(AgentDistPath(:,:,:,2),2),3))- squeeze(sum(sum(AgentDistPath3(:,:,:,2),2),3)) % assets, most are zeros, but then therre is a 0.0152 in the early values




%% Policy is same in both
%load temp.mat PolicyIndexesPath
%size(PolicyPath)
%size(PolicyIndexesPath)


%PolicyPath2=PolicyPath(1,:,:,:,:,1:end-1)+n_a(1)*(PolicyPath(2,:,:,:,:,1:end-1)-1);
%PolicyPath2=reshape(PolicyPath2,[N_a,N_z,N_j,T-1]);
%max(max(max(max(abs(PolicyPath2-PolicyIndexesPath))))) % zero


%%
%load temptt1.mat AgentDist AgeWeights AgeWeightsOld optaprime pi_z_J_sim exceptlastj exceptfirstj justfirstj jequalOneDist
%load temptt.mat AgentDist2 AgeWeights2 AgeWeightsOld2 optaprime2 pi_z_J_sim2 exceptlastj2 exceptfirstj2 justfirstj2 jequalOneDist2


%max(abs(AgentDist(:)-AgentDist2(:)))
%max(abs(AgeWeights(:)-AgeWeights2(:)))
%max(abs(AgeWeightsOld(:)-AgeWeightsOld2(:)))
%max(abs(optaprime(:)-optaprime2(:)))
%max(abs(pi_z_J_sim(:)-pi_z_J_sim2(:)))
%max(abs(exceptlastj(:)-exceptlastj2(:)))
%max(abs(exceptfirstj(:)-exceptfirstj2(:)))
%max(abs(justfirstj(:)-justfirstj2(:)))
%max(abs(jequalOneDist(:)-jequalOneDist2(:)))
