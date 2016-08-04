function [m,K,Kxsx,Kxsxs,Kxxs]=generateFromGP(gp,Xstar,predictWithNoise)
% generates data from a GP at points Xstar
%
% gp is a GP created using createGP
% Xstar the points at which to predict the output
% predictWithNoise [optional] if true then noise will be added to the
% output according to the sigmaN stored in gp

if nargin<3 || isempty(predictWithNoise)
    predictWithNoise=0;
end
if isempty(gp.X)    
    Kxsxs=gp.covFun(Xstar,Xstar,gp.covFunParams);
    if predictWithNoise
        Kxsxs=Kxsxs+gp.sigmaN^2.*eye(size(Kxsxs));
    end
    m=gp.fMean(Xstar)';
    K=Kxsxs;
else
    Kxsx=gp.covFun(Xstar,gp.X,gp.covFunParams);
    Kxsxs=gp.covFun(Xstar,Xstar,gp.covFunParams);
    if predictWithNoise
        Kxsxs=Kxsxs+gp.sigmaN^2.*eye(size(Kxsxs));
    end
    Kxxs=Kxsx';
    Ktmp=Kxsx*gp.Kxxinv;
    m=real(Ktmp*(gp.y'-gp.fMeanX'))+gp.fMean(Xstar)';
    K=Kxsxs-Ktmp*Kxxs;
end
%fstar=grand(m,K)';
end