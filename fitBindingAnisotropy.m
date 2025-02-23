function [fr,xForPlot,fittedAni]=fitBindingAnisotropy(cFluorescentPeptide,anisotropy,cLargeBinder)
% The program fits the anisotropy increase of a small peptide as it binds
% to a large target site.
% cFluorescentPeptide - concentration of the peptide whose anisotropy is measured
% anisotropy - anisotropy to be fitted
% cLargeBinder - concentration of the target the peptide binds to
%
% Model equation:
% r = rUnbound * (1-fBound) + rBound * fBound
% rUnbound (ru in the program), rBound (rb in the program) - anisotropies of the free and bound peptide, respectively
% fBound (fb in the program) - fraction of bound peptides
% The bound fraction is calculated by the 'simpleBindingWithDepletion'
% MATLAB program that calculates binding with ligand depletion.
% Output:
% fr - fit results
% xForPlot - concentration range of the large target for plotting
% fittedAni - calculated anisotropy at xForPlot
% Nov 22, 2023, email: peter.v.nagy@gmail.com, https://peternagyweb.hu
fo=fitoptions('method','nonlinearleastsquares');
ft=fittype('(1-simpleBindingWithDepletion(kd,cLargeBinder,cFluorescentPeptide)/cFluorescentPeptide)*ru + simpleBindingWithDepletion(kd,cLargeBinder,cFluorescentPeptide)/cFluorescentPeptide*rb',...
    'independent','cLargeBinder','problem','cFluorescentPeptide','options',fo);
% coeffnames(ft)
% coefficients: kd, rb, ru
fr=fit(cLargeBinder,anisotropy,ft,'problem',cFluorescentPeptide,'startpoint',[1 0.25 0.15],'lower',[0 0 0]);
figure;
plot(fr)
hold on;
plot(cLargeBinder,anisotropy,'bo');
xForPlot=linspace(min(cLargeBinder),max(cLargeBinder),100)';
fittedAni=feval(fr,xForPlot);