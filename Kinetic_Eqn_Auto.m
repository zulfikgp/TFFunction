% Kinetic equations for autoregulation
% Compute fold-change and response time
% For timing expression starts at constitutive level

clear

tic

T12 = 40; % Cell division
ga = log(2)/(60*T12); % Degradation rate
r1 = 0.1; % basal rate
P = 0.1; 
t_end =2/ga; % end time for solving ODEs to get response time
% vary t_end based on alpha, beta values. It get very slow for some values 

kpoff = 1;
kpon = 1;
kTFoff= 0.001;
kTFon = 0.0002; % sinlge TF binding rate


al = round(logspace(log10(0.001),log10(1000),10),4);% alpha values
be = round(logspace(log10(0.001),log10(1000),10),4); % beta values
Pa = P;  % round(logspace(log10(0.0001),log10(1000),20),4); % P values
FPT=zeros(length(be),length(al));

cnt = 0;

for i=1:length(be)
	for j=1:length(al)
	for k=1:length(Pa)		
		%k
		P = Pa(k);
		%kpon  = P*kpoff; % Tune kpon
		kpoff  = kpon/P; % Tune kpoff
		K = [kTFon kTFoff kpon kpoff];
		V2 = r1/kpoff;
		M0 = r1*P/((1+P+V2)*ga); % Constitutive expression

		%Steady state expression
		[Mmon sol1 P01 P10 P11 nsol] = YSS(r1,al(j),be(i),ga,K);
        Mss = Mmon + P01 + P11; % Total TF free+bound 
		fun = @(t,y) ODEs(t,y,K,al(j),be(i),r1,ga);

        % If expression starts at 0 use this
		% [T Y]=ode45(fun,[0 t_end],[0 0 0 0]);
		% ind = find(Y(:,1)>= 0.5*Mss);

        % or Start at constitutive level
		[T Y]=ode45(fun,[0 t_end],[M0 0 P/(1+P+V2) 0]); 
		if (Mss>=M0)
			ind = find(Y(:,1)> M0 + 0.5*(Mss-M0));
		else
			ind = find(Y(:,1) < M0 - 0.5*(M0-Mss));
		end

		if (length(ind)>0)
			FPT(i,j) = T(ind(1))/(T12*60); % Reponse time
		else
			%disp("No Sol")
		end
		FCT(i,j) = (Mss+P01+P11)/M0; % Fold-change of total TF
        %FCT(i,j) = Mss/M0; % Fold-change of free TF
        %Eff_TF(i,j) = Mss; % Equilibrium TF number

		%FCT(i,j,k) = Mss/M0;

		%Pol(i,k) = P01;
		%Col(i,k) = P11;
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
end
%return
% generate heatmap for Fold-change and response time
subplot(1,2,1)
imagesc(log10(FCT))
axis xy
subplot(1,2,2)
imagesc(FPT)
%imagesc(Eff_TF)
axis xy
clc
toc
%%%%%%%%%%%%%%ODes for auto-regulation %%%%%

function [dy] = ODEs(t,y,K,al,be,r1,ga)

dy = zeros(4,1);
% Rates for the reactions

r2 = al*r1;

kTFon = K(1,1); % Binding of single TF
kTFoff = K(1,2); % Unbinding

kpon = K(1,3); % Binding to decoy site
kpoff = K(1,4); % Unbinding decoy


f1 = 1-y(2)-y(3)-y(4); % Unbound state

dy(1) = r1*y(3) + r2*y(4) - ga*y(1) - y(1)*kTFon*(f1+y(3)) + kTFoff*(y(2)+y(4)/be); % Eqn for m

dy(2) = y(1)*kTFon*f1 - (kTFoff + kpon)*y(2) + (kpoff/be + r2)*y(4) - ga*y(2); % Occupied TF promoter

dy(3) = kpon*f1 + kTFoff*y(4)/be - kpoff*y(3) - y(1)*kTFon*y(3) - r1*y(3) + ga*y(4); % Occupied polymerase

dy(4) = kpon*y(2) + y(1)*kTFon*y(3) - (kpoff + kTFoff)*y(4)/be  - (r2+ga)*y(4); % cobound

end

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
