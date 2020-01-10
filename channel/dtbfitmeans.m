function [D,W,T_opt,T_opt_se,grad,hess,nlogl,Wpred]=dtb_fit_means(D,varargin)
% [D,W,T_opt,,T_opt_se,grad,hess]=dtb_fit_means(D,varargin)
% optimizes dtb (diffusion to bound) model to mean data with possible parameters:
%   B           bounds at �B 
%   kappa       multiplier on strenghth to get drift
%   mu          strength bias
%   y0          starting point offset as a multiplier on B  [-1 to +1]
%   tnd         non decision time
%   tnd_delta   non decision time offset for positive strengths 
%   var_coeff   changes  variance noise to be  1+var_coeff.*var_data 
%
% D is a  data structutre that contains
%           strength  signed  [e.g. -0.256]
%           choice   choice  [0 1]
%           rt   in seconds
%           include_for_rt - flag for which of the data to use for rt part only
%           [optional] var_data    see var_coeff above
%           [optional] theta  used for prediction when opt=0
% 
%  
%  flags
%   rt_only          whether to use rt for fitting only ([0] 1)
%   opt              whether to optimize or use parameters in D.theta (0 [1])
%   mu_opt           whether to optimize mu ([0] 1)
%   tnd_delta_opt    whether to optimize tnd_delta which is for pos ([0] 1)
%   var_coeff_opt    whether to optimize var_coeff ([0] 1)
%   y0_opt           whether to optimize y0 ([0] 1)
% 
%   B,mu,kappa,tnd,tnd_delta,var_coeff,y0 lets you fix these parameters by setting them 
%                    so set things to zero not to fit them
% 
% returns
%
%    D  original D appended with optimal theta and other usesful things
%    W  predictions of the data and other useful things
%    T_opt -  fit parameters as a structure
%    T_opt_se - s.e. of fit parameters as a structure
%    grad, hess and nlogl  - grdient, hessian and nlogl of optimal solution

pSet = inputParser;

addParameter(pSet,"B",NaN);
addParameter(pSet,"kappa",NaN);
addParameter(pSet,"tnd",NaN);

addParameter(pSet,"mu",0);
addParameter(pSet,"tnd_delta",0);
addParameter(pSet,"var_coeff",0);
addParameter(pSet,"y0",0);

addParameter(pSet,"tnd_delta_opt",0);
addParameter(pSet,"var_coeff_opt",0);
addParameter(pSet,"y0_opt",0);
addParameter(pSet,"mu_opt",0);

addParameter(pSet,"rt_only",0);
addParameter(pSet,"opt",1); 

addParameter(pSet,"printGrad",true,@islogical);

parse(pSet,varargin{:});
v2struct(pSet.Results);
%done unpacking

D.corr=(D.strength>0) .* (D.choice==1)  + (D.strength<0) .* (D.choice==0) + (D.strength==0) .* 0.5;

%initial parameter settings and lower and upper bounds
clear T
if isnan(B)
    T{1}.B = 0.8;         T{2}.B=0.0;        T{3}.B=2;
else
    T{1}.B = B;         T{2}.B=B;        T{3}.B=B;
end

if isnan(kappa)
    T{1}.kappa =1;      T{2}.kappa=0.1;    T{3}.kappa=40;
else
    T{1}.kappa =kappa;      T{2}.kappa=kappa;    T{3}.kappa=kappa;
end

if isnan(tnd)
    T{1}.tnd =0.35;      T{2}.tnd=0.020;     T{3}.tnd=max(D.rt);
else
    T{1}.tns =tnd;      T{2}.tnd=tnd;     T{3}.tnd=tnd;
end

if tnd_delta_opt
    T{1}.tnd_delta =0.05;      T{2}.tnd_delta=-0.5;     T{3}.tnd_delta=+0.5;
else
    T{1}.tnd_delta =tnd_delta;      T{2}.tnd_delta=tnd_delta;     T{3}.tnd_delta=tnd_delta;
end

if mu_opt
    T{1}.mu =0.1;        T{2}.mu=-1;      T{3}.mu=1;
else
    T{1}.mu =mu;        T{2}.mu=mu;      T{3}.mu=mu;
end

if  y0_opt
    T{1}.y0 =0.1;        T{2}.y0=-1;      T{3}.y0=1;
else
    T{1}.y0 =y0;        T{2}.y0=y0;      T{3}.y0=y0;
end

if  var_coeff_opt
    T{1}.var_coeff =0.1;        T{2}.var_coeff=-1;      T{3}.var_coeff=10;
else
    if var_coeff==0
        D.var_data=ones(size(D.strength));
    end
    T{1}.var_coeff =var_coeff;    T{2}.var_coeff=var_coeff;      T{3}.var_coeff=var_coeff;
end

[theta,theta_lo,theta_hi,P]=opt_pack(T,0.1);
options = optimoptions("fmincon","GradObj","on","Display","notify");

if opt
  hess_flag=0;
    [theta_opt,nlogl] = fmincon(@(theta) dtb_cost_means(theta,P,D,rt_only,hess_flag),theta,[],[],[],[],theta_lo,theta_hi,[],options);
    T_opt=opt_unpack(theta_opt,P);
else
    T_opt=opt_unpack(theta,P);
    theta_opt=D.theta_opt;
end

hess_flag=1;
[nlogl,grad,hess, W]=dtb_cost_means(theta_opt,P,D,rt_only,hess_flag);
theta_opt_se=sqrt(diag(inv(hess)))";
T_opt_se=opt_unpack(theta_opt_se,P);

if pSet.Results.printGrad  % mns added conditional
    fprintf("Norm of final gradient = %f\n",norm(grad))
end

%this generates predictions over and interpolated range
Dopt.strength=sort([-logspace(log10(0.032) ,log10(0.512) ,33) 0 logspace(log10(0.032) ,log10(0.512) ,33)])";
Dopt.choice=Dopt.strength>=0;
Dopt.include_for_rt=Dopt.choice==Dopt.choice;
Dopt.rt=rand(size(Dopt.choice));
Dopt.rt_std=rand(size(Dopt.choice));
Dopt.var_data=abs(Dopt.strength);

hess_flag=0;
[~,~,~, Wpred]=dtb_cost_means(theta_opt,P,Dopt,rt_only,hess_flag);

D.theta_opt=theta_opt;
D.theta_opt_se=theta_opt_se;

%%
function [nlogl grad hess R]=dtb_cost_means(theta,S,D,rt_only,hess_flag)

T=opt_unpack(theta,S);

drift=T.kappa*(D.strength+T.mu);
if hess_flag  %we only need to do this once at the end of the optimization to get s.e. on parameters
    [d1,d2,e1,e2,f1,f2,g1,g2] =  dtb_grad(T.B,T.kappa,T.mu,T.tnd,T.tnd_delta,T.var_coeff,T.y0,D.strength,D.choice,D.rt,D.rt_std,D.var_data,D.rt*0);
    
    d1=[d1 f1];
    d2=[d2 f2];
    e1=[e1 g1];
    e2=[e2 g2];
else
    [d1,d2,e1,e2] =  dtb_grad(T.B,T.kappa,T.mu,T.tnd,T.tnd_delta,T.var_coeff,T.y0,D.strength,D.choice,D.rt,D.rt_std,D.var_data,D.rt*0);
end

%ds are RT stuff
d1(isnan(d1))=0;
d2(isnan(d2))=0;
%es are prob rightward choice 
e1(isnan(e1))=0;
e2(isnan(e2))=0;

s=abs(drift)<1e-4;
d=(bsxfun(@times,d1,(1-s))+bsxfun(@times,d2,s));
e=(bsxfun(@times,e1,(1-s))+bsxfun(@times,e2,s));

R.strength=D.strength;
R.choice=d(:,2);
R.corr=(D.strength>0) .* R.choice  + (D.strength<0) .* (1-R.choice) + (D.strength==0) .* 0.5;
R.rt=d(:,3);
R.rt_pos=d(:,4);
R.rt_neg=d(:,5);

d=bsxfun(@times,d,D.include_for_rt);
R.rt_nlogl=sum(d(:,1));
R.p_nlogl=sum(e(:,1));
R.nlogl_item_rt = d(:,1);
R.nlogl_item_p = e(:,1);

if ~rt_only
    d=d+ e;
end

ds=sum(d);
grad=ds(6:12);

if hess_flag
    hess=reshape(ds(13:end),7,7);
else
    hess=[];
end
%remove the parameters that we have fixed
s=find(row_sum(S.tied));
p=cumsum(S.n);
s=p(s);
grad(s-1)=grad(s-1)+grad(s);

s=find(row_sum(S.fixed|S.tied));
p=cumsum(S.n);
s=p(s);
grad(s)=[];

if hess_flag
    hess(s,:)=[];
    hess(:,s)=[];
end
nlogl=ds(1);
%fprintf("nlogl = %f\n",nlogl)
R.bic= 2*nlogl + length(theta)*log(sum(D.strength));
R.nlogl=nlogl;



