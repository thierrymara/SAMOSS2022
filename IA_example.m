if exist ('OCTAVE_VERSION', 'builtin')
    warning('off')
    pkg load statistics
    %pkg load gnuplot
    struct_levels_to_print(0)
end
close all
clear all

N=1000;
d=3;

A=rand(N,d);yA=A(:,1)+10*(A(:,2)-0.5).*(A(:,3)-0.5);%Model run: x1+10(x2-0.5)(x3-0.5)
B=rand(N,d);yB=B(:,1)+10*(B(:,2)-0.5).*(B(:,3)-0.5);
for i=1:d
  Ab=A;Ab(:,i)=B(:,i);
  Ba=B;Ba(:,i)=A(:,i);
  yAb=Ab(:,1)+10*(Ab(:,2)-0.5).*(Ab(:,3)-0.5);
  yBa=Ba(:,1)+10*(Ba(:,2)-0.5).*(Ba(:,3)-0.5);
  S(i) = 2*sum((yA-yAb).*(yBa-yB))/sum((yBa-yAb).^2+(yA-yB).^2);
  T(i) = sum((yA-yAb).^2+(yBa-yB).^2)/sum((yBa-yAb).^2+(yA-yB).^2);
endfor
S
T
