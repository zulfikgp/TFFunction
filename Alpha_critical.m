% Find condition for alpha to switch from increasing to decreasing curve when beta is tuned
% 
clear all

tic
kpoff = 1;
r1 = 0.1;
T12 = 40; % Cell division in minutes
ga = log(2)/(60*T12); % Degradation rate


kTFoff= 0.001;
kTFon = 0.0002*kTFoff*1000;
be1 = 0.01; % low beta value
be2 = 100; % High beta value

% Compare Thermodynamic and Full model for a given P and varying R

P = 1; % Choose a P value
RArr = logspace(log10(0.0001),log10(1000),100); % Use an array of R
al = P./(1+P+RArr); % alpha critical from thermodynamic model

% Full model

alf = logspace(log10(0.00001),log10(2),1000); % Array for alpha values

for j=1:length(RArr)
check=0;
i=0;
while check==0 | i==length(alf)
	i=i+1;
	FC1 = FullModel(alf(i),be1,r1,kpoff,P,kTFoff,RArr(j)); % FC for beta 1
	FC2 = FullModel(alf(i),be2,r1,kpoff,P,kTFoff,RArr(j)); % FC for beta 2
	%FC1 = FullModel(alf(i),be1,r1,kpoff,PArr(j),kTFoff,R);
	%FC2 = FullModel(alf(i),be2,r1,kpoff,PArr(j),kTFoff,R);
	DFC = FC2-FC1;
	if DFC>0
		check=1;
		alc(j,1) = alf(i); % alpha critical for switching
	end
end
end
subplot(1,2,1)
loglog(RArr(1:5:end),al(1:5:end),'o')
hold on
loglog(RArr,alc)
xlabel("R");
ylabel("alpha critical");

%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%Compare thermodynamic model and autoregulation %%%%%%%%%%

R1 = 0.01; R2 = 100; % Choose some value for effective TF concentration
PArr = logspace(log10(0.001),log10(1),10); % Array of P
al1 = PArr./(1+PArr+R1); % alpha critical from thermodynamic model
al2 = PArr./(1+PArr+R2); % alpha critical from thermodynamic model
subplot(1,2,2)
loglog(PArr,al1,'r--')
hold on
loglog(PArr,al2,'r')
xlabel("R");
ylabel("alpha critical");

% Autoregulation
kTFon = 0.01; % Choose some value for kTFon per TF
for j=1:length(PArr)
check=0;
i=0;
while check==0 | i==length(alf)
	i=i+1;
	%FC1 = AutoModel(alf(i),be1,r1,kpoff,P,kTFoff,R(j),ga);% Use if varying R
	%FC2 = AutoModel(alf(i),be2,r1,kpoff,P,kTFoff,R(j),ga);
	FC1 = AutoModel(alf(i),be1,r1,kpoff,PArr(j),kTFoff,kTFon,ga);
	FC2 = AutoModel(alf(i),be2,r1,kpoff,PArr(j),kTFoff,kTFon,ga);
	DFC = FC2-FC1;
	if DFC>0
		check=1;
		alc_auto(j,1) = alf(i); % alpha critical for auto-gene
	end
end
end

loglog(PArr,alc_auto,'b')
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to compute fold-change from full model

function FC = FullModel(al,be,r1,kpoff,P,kTFoff,R)

	V2 = r1/kpoff;
	kpon = P*kpoff;
	V1 = kTFoff/kpoff;
	
	N1 = be*R*P*( V1+V2+1+P+R*V1);
	N2 = (1+P+V2)*(1+R)*(1+V1+al*be*V2) + (1+P+V2)*P*(1+be*R) + R*V1*(1+R)*( 1 + P*be +al*be*V2 ); 
	P11 = N1/N2;
	P01 = R/(1+R) - (1+be*R)*P11/(be*(1+R));
	P10 = P/(1+P+V2) - (1+be*(al*V2+P))*P11/(be*(1+P+V2));

	FC = 1 + (al*be*(1+P)-P*be-1)*P11/(P*be);

end

% Function to compute Fold-change in auto-regulation
function FC = AutoModel(al,be,r1,kpoff,P,kTFoff,kTFon,ga)

	kpon = P*kpoff;
	r2 = r1*al;
	V2 = r1/kpoff;
	syms m P01 P10 P11;
	
	eqn1 = r1*P10 + r2*P11 - ga*m;
	eqn2 = m*kTFon*(1-P01-P10-P11) - (kTFoff + kpon)*P01 + (kpoff/be +r2)*P11;
	eqn3 = kpon*(1-P01-P10-P11) - (kpoff + m*kTFon)*P10 + kTFoff*P11/be - r1*P10 + ga*P11;
	eqn4 = kpon*P01 + m*kTFon*P10 - ((kTFoff + kpoff)/be + r2 + ga)*P11;

	%sol = solve([eqn1, eqn2, eqn3, eqn4], [m, P01, P10, P11], 'Real',true);
	sol = vpasolve([eqn1, eqn2, eqn3, eqn4], [m, P01, P10, P11]);
	sol1 = double(vpa(sol.m));
	sol2 = double(vpa(sol.P01));
	sol3 = double(vpa(sol.P10));
	sol4 = double(vpa(sol.P11));

	Mss = 0;
	Cons = r1*P/((1+P+V2)*ga); % Constitutive expression

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
	FC = Mss/Cons;
end
