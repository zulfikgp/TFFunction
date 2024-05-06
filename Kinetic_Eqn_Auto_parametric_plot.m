% Kinetic equations for autoregulation
% Compute fold-change for extreme values of P and use that to find region
% where the TF switcesh from activator to repressor and vice versa as
% promoter strength is tuned
% Also find other 4 regimes in the alpha-beta space

clear

tic

T12 = 40; % Cell division
ga = log(2)/(60*T12); % Degradation rate
r1 = 0.1; % basal rate
kpoff = 1;
kpon = 1;
kTFoff= 0.001;
kTFon = 0.0002; % sinlge TF binding rate


al = round(logspace(log10(0.01),log10(100),100),4);% alpha values
be = round(logspace(log10(0.01),log10(100),100),4); % beta values
P1 = 0.001;  % low P value
P2 = 1000; % high P value

cnt = 0;

for i=1:length(be)
	for j=1:length(al)
		%Steady state FC for low P
		kpon  = P1*kpoff; % Tune kpon
     	%kpoff  = kpon/P1; % Tune kpoff
        V2 = r1/kpoff;
		M0 = r1*P1/((1+P1+V2)*ga); % Constitutive expression		
        K = [kTFon kTFoff kpon kpoff];
	    [Mmon sol1 P01 P10 P11 nsol] = YSS(r1,al(j),be(i),ga,K);
        Mss = Mmon + P01 + P11; % Total TF free+bound
        FCl = Mss/M0;
        
        %Steady state expression and FC for high P
        kpon  = P2*kpoff; % Tune kpon
     	%kpoff  = kpon/P2; % Tune kpoff
        V2 = r1/kpoff;
		M0 = r1*P2/((1+P2+V2)*ga); % Constitutive expression
        K = [kTFon kTFoff kpon kpoff];
	    [Mmon sol1 P01 P10 P11 nsol] = YSS(r1,al(j),be(i),ga,K);
        Mss = Mmon + P01 + P11; % Total TF free+bound
        FCh = Mss/M0;

        % Find regions of interest
        if FCl < 1 & FCh < FCl
            ParPlot(i,j)=1;
        elseif FCl >1 & FCh > FCl
            ParPlot(i,j)=2;
        elseif FCl < 1 & FCh > 1
            ParPlot(i,j)=3;
        elseif FCl > 1 & FCh < 1
            ParPlot(i,j)=4;
        elseif FCl < 1 & FCh > FCl
            ParPlot(i,j)=5;
        else
            ParPlot(i,j)=6;
        end


		%FCT(i,j) = (Mss+P01+P11)/M0; % Fold-change of total TF
		%subplot(1,2,1)
		%plot(T/(60*T12),Y(:,1)/Mss,'--')
		%hold on
		%subplot(1,2,2)
		%plot(T/(60*T12),Y(:,1),'--')
		%pause(1)
		%hold on
		cnt = cnt +1;
    end
end
%return
% generate heatmap for Fold-change and response time
imagesc(log10(ParPlot))
axis xy
clc
toc


%%% Steady state for auto-regulation

function [Mss sol1 sol2 sol3 sol4 nsol] = YSS(r1, al, be, ga, K)
	kTFon = K(1,1); kTFoff = K(1,2); kpon = K(1,3); kpoff = K(1,4);
	r2 = r1*al;
	syms m P01 P10 P11;
	
	eqn1 = r1*P10 + r2*P11 - ga*m  - m*kTFon*(1-P01-P11) + kTFoff*(P11/be + P01);% m
	eqn2 = m*kTFon*(1-P01-P10-P11) - (kTFoff + kpon)*P01 + (kpoff/be +r2)*P11 - ga*P01; % TF bound
	eqn3 = kpon*(1-P01-P10-P11) - (kpoff + m*kTFon)*P10 + kTFoff*P11/be - r1*P10 + ga*P11; %Pol bound
	eqn4 = kpon*P01 + m*kTFon*P10 - ((kTFoff + kpoff)/be + r2 + ga)*P11; % Co-bound

	%sol = solve([eqn1, eqn2, eqn3, eqn4], [m, P01, P10, P11], 'Real',true);
	sol = vpasolve([eqn1, eqn2, eqn3, eqn4], [m, P01, P10, P11]);
	sol1 = double(vpa(sol.m));
	sol2 = double(vpa(sol.P01));
	sol3 = double(vpa(sol.P10));
	sol4 = double(vpa(sol.P11));

	Mss = 0;
	nsol = length(find(real(sol1)>0));
	if nsol>1
		[al be K]
	end
	if (length(find(real(sol1)>0))>1)
		disp('Multiple solution found')
	end
	%disp('# of Soluton = ')
	%[length(sol1) sol1(1) sol1(2) sol1(3)]
	for i=1:length(sol1)
		if (sol1(i))>=0
			Mss = real(sol1(i));
			sol2 = real(sol2(i));
			sol3 = real(sol3(i));
			sol4 = real(sol4(i));
		end
	end
end


%numPlots = 9
%ramp = linspace(0, 0.75, numPlots);
%listOfGrayColors = [ramp; ramp; ramp]';
%for i=1:9
%plot(al,FPT(i,:,2),'color',listOfGrayColors(i, :),'LineWidth',1.5)
%hold on
%end
