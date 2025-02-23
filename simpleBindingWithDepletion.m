function [TL,T,L]=simpleBindingWithDepletion(kd,Ttot,Ltot)
% Ligand (L) binds to a target (T) with a dissociation constant of Kd.
TL=(1/2).*(kd+Ltot+Ttot+(-1).*((kd+Ltot).^2+2.*(kd+(-1).*Ltot).*Ttot+Ttot.^2).^(1/2));
T=(1/2).*((-1).*kd+(-1).*Ltot+Ttot+((kd+Ltot).^2+2.*(kd+(-1).*Ltot).*Ttot+Ttot.^2).^(1/2));
L=(1/2).*((-1).*kd+Ltot+(-1).*Ttot+((kd+Ltot).^2+2.*(kd+(-1).*Ltot).*Ttot+Ttot.^2).^(1/2));