% Timing of gene regulated by fixed TF
% Plot both fold-change and Response time heatmap
%

clear all
close all
%hold off
tic

T12 = 40; % Cell division
ga = log(2)/(60*T12);
TEnd = 3*T12*60;
r1 = 0.1;
P = 0.62;

kpoff = 1;
kpon  = P*kpoff;

kTFoff= 0.001;
kTFon = 0.0002; % kTFon is the effective on rate for all the TFs
V2 = r1/kpoff;

M0 = r1*P/((1+P+V2)*ga); % Constitutive expression

be = logspace(log10(0.01),log10(1000),50); %[0.001 0.01 0.1 1 100];
al = logspace(log10(0.01),log10(1000),50);%[0.001 0.01 0.1 1 100];
%Pa = logspace(log10(0.01),log10(1000),100);
%kArr = logspace(log10(0.0000001),log10(1000),100);

cnt = 0;
k=0;
for kk=1:1%length(kArr)
    %kTFon = kArr(kk);
for i=1:length(be)
	for j=1:length(al)
		%j=i;
	
		V2 = r1/kpoff;
		R = kTFon/kTFoff;
		P = kpon/kpoff;
		%P = Pa(kk);
		kpon = P*kpoff;
		V1 = kTFoff/kpoff;
		K = [kTFon kTFoff kpon kpoff];

		N1 = be(i)*R*P*( V1+V2+1+P+R*V1);
		N2 = (1+P+V2)*(1+R)*(1+V1+al(j)*be(i)*V2) + (1+P+V2)*P*(1+be(i)*R) + R*V1*(1+R)*( 1 + P*be(i) +al(j)*be(i)*V2 ); 
		P11 = N1/N2;
		P01 = R/(1+R) - (1+be(i)*R)*P11/(be(i)*(1+R));
		P10 = P/(1+P+V2) - (1+be(i)*(al(j)*V2+P))*P11/(be(i)*(1+P+V2));

		FC = 1 + (al(j)*be(i)*(1+P)-P*be(i)-1)*P11/(P*be(i));

		Mss = M0 + (r1/ga)*P11*(al(j) - (1+be(i)*(P+al(j)*V2))/(be(i)*(1+P+V2))); % mrna expression

		%[Y(end,1)  Mss Y(end,1)/M0 FC]
		%[Mss Y(end,1)/M0 Y(end,2:4)]
		%[Mss/M0 FC Y(end,1)/M0]

		fun = @(t,y) ODEs(t,y,K,al(j),be(i),r1,ga);

%		[T Y]=ode45(fun,[0 4000],[0 0 0 0]);

		[T Y]=ode45(fun,[0 TEnd],[M0 0 P/(1+P+V2) 0]); % Start at constitutive level
		%plot(T,Y(:,1))
		%pause(1)

		if (Mss>=M0)
			ind = find(Y(:,1)> M0 + 0.5*(Mss-M0));
		else
			ind = find(Y(:,1) < M0 - 0.5*(M0-Mss));
		end
		if (length(ind)>0)
			%FPT(j,kk) = T(ind(1))/(T12*60);
			FPT(i,j) = T(ind(1))/(T12*60);
		else
			disp("No Sol")
		end
		%FCT(j,kk) = FC;
		FCT(i,j) = FC; % Fold-change in expression
		Pol(j,kk) = P10;
		TFb(j,kk) = P01;
		Col(j,kk) = P11;

		%subplot(1,2,1)
		%plot(T/(60*T12),Y(:,1)/Mss,'--')
		%hold on
		%subplot(1,2,2)
		%plot(T/(60*T12),Y(:,1),'--')
		%pause(1)
		%hold on
		%cnt = cnt +1
%[i j]
    end
end
end
subplot(1,2,1)
imagesc(log10(FCT))
axis xy
colorbar
title("Fold-change")
subplot(1,2,2)
imagesc(FPT)
axis xy
colorbar
title("Response time")
clc
toc

%%%%%%%%%%%%%%%%%%%

function [dy] = ODEs(t,y,K,al,be,r1,ga)

dy = zeros(4,1);
% Rates for the reactions

r2 = al*r1;

kTFon = K(1,1); % Binding
kTFoff = K(1,2); % Unbinding

kpon = K(1,3); % Binding to decoy site
kpoff = K(1,4); % Unbinding decoy


f1 = 1-y(2)-y(3)-y(4); %Free

dy(1) = r1*y(3) + r2*y(4) - ga*y(1); %

dy(2) = kTFon*f1 - (kTFoff + kpon)*y(2) + (kpoff/be + r2)*y(4); % Occupied TF promoter

dy(3) = kpon*f1 + kTFoff*y(4)/be - kpoff*y(3) - kTFon*y(3) - r1*y(3); % Occupied polymerase

dy(4) = kpon*y(2) + kTFon*y(3) - (kpoff + kTFoff)*y(4)/be  - r2*y(4); % cobound

end

%%%%%%%%%%%%


%numPlots = 9
%ramp = linspace(0, 0.75, numPlots);
%listOfGrayColors = [ramp; ramp; ramp]';
%for i=1:9
%plot(al,FPT(i,:,2),'color',listOfGrayColors(i, :),'LineWidth',1.5)
%hold on
%end
