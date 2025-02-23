function result=fitManyAniSatCurves(anidata,elementForBoundAni,ccPeptide)
% The program can be used for analyzing anisotropy measurements in which a
% the binding of a small fluorescent binder to a large target is measured.
% The program determines the min and max anisotropies and the dissociation
% constants for each curve by global fitting.
% anidata - a structure array with each element having the following fields:
%       - data: an n x 2 array, 1st column - concentration, 2nd column: anisotropy
%       - id: name of the sample
%       - remark: what its name suggests
% elementForBoundAni - this dataset will be used for finding the max
%       anisotropy, i.e., the anisotropy of the bound peptide (if other datasets don't reach saturation)
%       If empty, anisotropy of the free and bound peptide will also be fit globally.
% ccPeptide - concentration of the fluorescent peptide that binds to the antibody
%
% result.fittedKds
% result.freeAni
% result.boundAni
% result.fittedAnis
% result.ccForPlot
%
% Written by Peter Nagy
% Nov 21, 2023, email: peter.v.nagy@gmail.com, https://peternagyweb.hu
if ~isempty(elementForBoundAni) % anisotropy of the free and bound peptide is found from one of the datasets
    % Find max anisotropy
    fr=fitBindingAnisotropy(ccPeptide,anidata(elementForBoundAni).data(:,2),anidata(elementForBoundAni).data(:,1));
    % fr parameters: kd, rb, ru
    ci=confint(fr);
    result.freeAni=fr.ru;
    result.freeAni_conf95=[ci(1,3) ci(2,3)];
    result.boundAni=fr.rb;
    result.boundAni_conf95=[ci(1,2) ci(2,2)];
    % global fitting of Kds using the anistropy of the free and bound peptide from the previous fit
    fitFunc=@(Kds) squaredDev(Kds,anidata,ccPeptide,result.freeAni,result.boundAni);
    [result.fittedKds,~,~,~,~,hessian]=fminunc(fitFunc,ones(numel(anidata),1));
    numOfPoints=sum(arrayfun(@(x) size(x.data,1),anidata));
    errorVariance=squaredDev(result.fittedKds,anidata,ccPeptide,result.freeAni,result.boundAni)/(numOfPoints-numel(result.fittedKds));
    result.fittedKds_SD=diag(real(sqrt(errorVariance*inv(hessian))));
else % global fitting of the anisotropy of the free and bound peptide
    fitFuncTotalGlobal=@(params) squaredDevTotalGlobal(params,anidata,ccPeptide);
    % params: 1 ... n : Kds
    % n+1: anisotropy of the free peptide
    % n+2: anisotropy of the bound peptide
    minAni=min(cell2mat(arrayfun(@(x) x.data(:,2),anidata,'UniformOutput',false)'));
    maxAni=max(cell2mat(arrayfun(@(x) x.data(:,2),anidata,'UniformOutput',false)'));
    [fittedParams,~,~,~,~,hessian]=fminunc(fitFuncTotalGlobal,[ones(numel(anidata),1);minAni;maxAni]);
    numOfPoints=sum(arrayfun(@(x) size(x.data,1),anidata));
    errorVariance=squaredDevTotalGlobal(fittedParams,anidata,ccPeptide)/(numOfPoints-numel(fittedParams));
    parameterSDs=diag(real(sqrt(errorVariance*inv(hessian))));
    result.fittedKds=fittedParams(1:numel(anidata));
    result.fittedKds_SD=parameterSDs(1:numel(anidata));
    result.freeAni=fittedParams(end-1);
    result.freeAni_SD=parameterSDs(end-1);
    result.boundAni=fittedParams(end);
    result.boundAni_SD=parameterSDs(end);
end
% plot the results
result.fittedAnis=cell(numel(anidata),1);
figure;
minCc=min(cell2mat(arrayfun(@(x) x.data(:,1),anidata,'UniformOutput',false)'));
maxCc=max(cell2mat(arrayfun(@(x) x.data(:,1),anidata,'UniformOutput',false)'));
result.ccForPlot=linspace(minCc,maxCc,100)';
colorCodes='bgrcmk';
for i=1:numel(anidata)
    currentColorCode=colorCodes(mod(i-1,6)+1);
    plot(anidata(i).data(:,1),anidata(i).data(:,2),[currentColorCode,'o']);
    hold on;
    result.fittedAnis{i}=(1-simpleBindingWithDepletion(result.fittedKds(i),result.ccForPlot,ccPeptide)/ccPeptide)*result.freeAni + simpleBindingWithDepletion(result.fittedKds(i),result.ccForPlot,ccPeptide)/ccPeptide*result.boundAni;
    plot(result.ccForPlot,result.fittedAnis{i},[currentColorCode,'-']);
end
xlabel('Concentration of antibody');
ylabel('Anisotropy');

function ss=squaredDev(Kds,anidata,ccPeptide,freeAni,boundAni)
ss=0;
for i=1:numel(Kds)
    calcAni=(1-simpleBindingWithDepletion(Kds(i),anidata(i).data(:,1),ccPeptide)/ccPeptide)*freeAni + simpleBindingWithDepletion(Kds(i),anidata(i).data(:,1),ccPeptide)/ccPeptide*boundAni;
    ss=ss+sum((calcAni-anidata(i).data(:,2)).^2);
end

function ss=squaredDevTotalGlobal(params,anidata,ccPeptide)
ss=0;
numOfSamples=numel(anidata);
Kds=params(1:numOfSamples);
freeAni=params(numOfSamples+1);
boundAni=params(numOfSamples+2);
for i=1:numOfSamples
    calcAni=(1-simpleBindingWithDepletion(Kds(i),anidata(i).data(:,1),ccPeptide)/ccPeptide)*freeAni + simpleBindingWithDepletion(Kds(i),anidata(i).data(:,1),ccPeptide)/ccPeptide*boundAni;
    ss=ss+sum((calcAni-anidata(i).data(:,2)).^2);
end