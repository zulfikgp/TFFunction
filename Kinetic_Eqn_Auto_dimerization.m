% Kinetic equations for autoregulation where TF dimerizes
% Compute fold-change and response time
% For timing expression starts at constitutive level

clear
clc

tic

T12 = 40; % Cell division
ga = log(2)/(60*T12); % Degradation rate
r1 = 0.1; % basal rate
P = 0.06; 
t_end = 1*T12*60; % end time for solving ODEs to get response time
kdim = 0.1;% Dimerization rate
kmon = 1; % monomerization rate
% vary t_end based on alpha, beta values. It get very slow for some values 

kpoff = 1;
kpon = 1;
kTFoff= 0.001;
kTFon = 0.0002; % sinlge TF binding rate


al = round(logspace(log10(0.001),log10(1000),10),4);% alpha values
be = round(logspace(log10(0.001),log10(1000),10),4); % beta values
Pa = 0.1;%round(logspace(log10(0.0001),log10(1000),50),4); % P values
%FPT=zeros(length(be),length(al));

cnt = 0;

for i=1:length(be)
	for j=1:length(al)
	for k=1:length(Pa)		
		%k
		P = Pa(k);
		%kpon  = P*kpoff; % Tune kpon
		kpoff  = kpon/P; % Tune kpoff
		K = [kTFon kTFoff kpon kpoff kdim kmon];
		V2 = r1/kpoff;
		M0 = r1*P/((1+P+V2)*ga); % Constitutive expression

		%Steady state expression
		[TFm P01 P10 P11 TFd nsol] = YSS(r1,al(j),be(i),ga,K);
        Mss = 2*(TFd + P01 + P11 ) + TFm; % Steady state monomer+dimers
		%Mss = 0;
        fun = @(t,y) ODEs(t,y,K,al(j),be(i),r1,ga);

        % If expression starts at 0 use this
		%[T Y]=ode45(fun,[0 t_end],[0 0 0 0 0]);
		% ind = find(Y(:,1)>= 0.5*Mss);

        % or Start at constitutive level
		[T Y]=ode45(fun,[0 t_end],[M0 0 P/(1+P+V2) 0 0]); 
		if (Mss>=M0)
			ind = find(Y(:,1)+2*(Y(:,5)+Y(:,2)+Y(:,4))> M0 + 0.5*(Mss-M0));
		else
			ind = find(Y(:,1)+2*(Y(:,5)+Y(:,2)+Y(:,4)) < M0 - 0.5*(M0-Mss));
		end

		if (length(ind)>0)
			FPT(i,j) = T(ind(1))/(T12*60); % Reponse time
        else
            FPT(i,j) = 0;
			disp("No Sol")
		end
		FCT(i,j) = Mss/M0; % Fold-change
        %FCT(k) = Mss/M0; % Fold-change
        %plot(T,(Y(:,1)+2*(Y(:,5)+Y(:,2)+Y(:,4)))/Mss)
        %hold on
        %title("alpha = "+al(j)+" beta = "+be(i))
        %pause(1)
		cnt = cnt +1;
    end
    %loglog(Pa,FCT)
    %hold on
    %pause(1)
end
end
%return
% generate heatmap for Fold-change and response time
figure
subplot(1,2,1)
imagesc(log10(FCT))
axis xy
subplot(1,2,2)
imagesc(FPT)
%imagesc(Eff_TF)
axis xy
toc
%%%%%%%%%%%%%%ODes for auto-regulation %%%%%

function [dy] = ODEs(t,y,K,al,be,r1,ga)

dy = zeros(5,1);
% Rates for the reactions

r2 = al*r1;

kTFon = K(1,1); % Binding of single TF
kTFoff = K(1,2); % Unbinding

kpon = K(1,3); % Binding to decoy site
kpoff = K(1,4); % Unbinding decoy

kdim = K(1,5);
kmon = K(1,6);


f1 = 1-y(2)-y(3)-y(4); % Unbound state

dy(1) = r1*y(3) + r2*y(4) - ga*y(1) + 2*kmon*y(5)-2*kdim*y(1)^2; % Eqn for m

dy(2) = y(5)*kTFon*f1 - (kTFoff + kpon)*y(2) + (kpoff/be + r2)*y(4) - ga*y(2); % Occupied TF promoter

dy(3) = kpon*f1 + kTFoff*y(4)/be - kpoff*y(3) - y(5)*kTFon*y(3) - r1*y(3) + ga*y(4); % Occupied polymerase

dy(4) = kpon*y(2) + y(5)*kTFon*y(3) - (kpoff + kTFoff)*y(4)/be  - (r2+ga)*y(4); % cobound

dy(5) = -kmon*y(5) + kdim*y(1)^2 - ga*y(5) - y(5)*kTFon*(f1+y(3)) + kTFoff*(y(2)+y(4)/be); %Dimerized TF

end

%%% Steady state for auto-regulation

function [TFm P01 P10 P11 TFd nsol] = YSS(r1, al, be, ga, K)
	kTFon = K(1,1); kTFoff = K(1,2); kpon = K(1,3); kpoff = K(1,4); 
    kmon = K(1,5); kdim = K(1,6);
	r2 = r1*al;
	syms m P01 P10 P11 mdim;
	
	eqn1 = r1*P10 + r2*P11 - ga*m + 2*kmon*mdim - 2*kdim*m^2 ;
	eqn2 = mdim*kTFon*(1-P01-P10-P11) - (kTFoff + kpon)*P01 + (kpoff/be +r2)*P11 - ga*P01;
	eqn3 = kpon*(1-P01-P10-P11) + kTFoff*P11/be - (kpoff + mdim*kTFon)*P10 - r1*P10 + ga*P11;
	eqn4 = kpon*P01 + mdim*kTFon*P10 - ((kTFoff + kpoff)/be + r2 + ga)*P11;
    eqn5 = -kmon*mdim + kdim*m^2 - ga*mdim - mdim*kTFon*(1-P01-P11) + kTFoff*(P11/be + P01); % Dimerized TF


	%sol = solve([eqn1, eqn2, eqn3, eqn4 eqn5], [m, P01, P10, P11, mdim], 'Real',true);
    %sol1 = double(sol.m) % Monomer
    %ind = find(sol1>0);
    %sol1 = sol1(ind);
	%sol2 = double(sol.P01(ind)); % TF bound
	%sol3 = double(sol.P10(ind)); % pol bound
	%ol4 = double(sol.P11(ind)); % Cobound
    %sol5 = double(sol.mdim(ind)); % Dimer

	sol = vpasolve([eqn1, eqn2, eqn3, eqn4 eqn5], [m, P01, P10, P11, mdim]);
	%double(vpa(sol.m))
    sol1 = double(vpa(sol.m)); % Monomer
	sol2 = double(vpa(sol.P01)); % TF bound
	sol3 = double(vpa(sol.P10)); % pol bound
	sol4 = double(vpa(sol.P11)); % Cobound
    sol5 = double(vpa(sol.mdim)); % Dimer

	nsol = find(sol1 == real(sol1) & real(sol1)>0);
    %sol1(nsol)
	
	if (length(nsol)>1)
		disp('Multiple solution found');
        TFm = 0;
        P01 = 0;
        P10 = 0;
        P11 = 0;
        TFd = 0;
    else
        TFm = sol1(nsol); % Free TF monomer
        P01 = sol2(nsol); % TF bound
        P10 = sol3(nsol); %POl bound
        P11 = sol4(nsol); % Cobound
        TFd = sol5(nsol); % TF dimer
    end
	%disp('# of Soluton = ')
	%[length(sol1) sol1(1) sol1(2) sol1(3)]
	
end


%numPlots = 9
%ramp = linspace(0, 0.75, numPlots);
%listOfGrayColors = [ramp; ramp; ramp]';
%for i=1:9
%plot(al,FPT(i,:,2),'color',listOfGrayColors(i, :),'LineWidth',1.5)
%hold on
%end
