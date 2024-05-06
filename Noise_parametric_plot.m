% Code to generate Figure 5. Noise for different alpha/beta values

clear all
close all
tic

bea = logspace(log10(0.001),log10(1000),30); % Array for beta
ala = logspace(log10(0.001),log10(1000),30);% Array for alpha

T12 = 40; % Cell division
ga = log(2)/(60*T12); % Degradation rate
r1 = 0.1; % Basal production rate
kTFoff= 0.001;
kTFon = 0.02;

Parr = [0.001, 0.01, 0.1, 1];

for lP=1:4
    subplot(2,2,lP)
    P=Parr(lP);
    kpoff = 1;
    kpon  = P*kpoff;

    V2 = r1/kpoff;
    NTF = 1;
    TFonArr = logspace(log10(0.0000001),log10(1000000),1000); % TF on rate

    cnt = 0;
    for i=1:length(bea)
    for j=1:length(ala)
	    cnt = cnt+1;
	    kTFon = 1;%*kTFoff*10^(k-1);

	    kpon  = P*kpoff; % Tune kpon
	    %kpoff  = kpon/P; % Tune kpoff

	    be = bea(i);
	    al = ala(j);

	    [M0 Mss FC FCmax(cnt,1) CVexact CVmax(cnt,1) CVcon] = Find_CV(al, be, r1, kpon, kpoff, kTFon, kTFoff, ga);
	    % Vary TF on rate
        for kk=1:length(TFonArr)
		    kTFon = TFonArr(kk);
		    [M0k Mssk FCk(cnt,kk) FCmaxk(cnt,1) CVexactk(cnt,kk) CVmaxk(cnt,1) CVcon] = Find_CV(al, be, r1, kpon, kpoff, kTFon, kTFoff, ga);
	    end
	    alarr(cnt,1) = al;
	    bearr(cnt,1) = be;
	    [y l w p] = findpeaks(CVexactk(cnt,:),'MinPeakProminence',0.001);

	    v = find(islocalmin(CVexactk(cnt,:),'MinProminence',0.01));
	    if length(y)>0
		    locs(i,j) = TFonArr(l);
		    %if y< CVmaxk(cnt,1)
		    if  length(v)==0
			    %subplot(1,2,2) 
			    loglog(al,be,'r*')
		    	hold on
		    	%subplot(1,2,1)
	   	    	%loglog(FCk(cnt,:),CVexactk(cnt,:)/CVmaxk(cnt,1),'m')
		    	%loglog(FCk(cnt,:),CVexactk(cnt,:)/CVcon,'r')
		    	hold on
		    	%waitforbuttonpress
		    	%pause(10)
	    	else
	    		%subplot(1,2,2) 
		    	loglog(al,be,'m*')
			    hold on
			    %subplot(1,2,1)
		    	%loglog(FCk(cnt,:),CVexactk(cnt,:)/CVmaxk(cnt,1),'r')
		    	%loglog(FCk(cnt,:),CVexactk(cnt,:)/CVcon,'m')
		    	hold on
	    		%waitforbuttonpress
	    		%pause(10)
	    	end	
    	else
    		locs(i,j) = 0;
	    	%subplot(1,2,2) 
	    	loglog(al,be,'b*')
		    hold on
	    	%subplot(1,2,1)
	    	%loglog(FCk(cnt,:),CVexactk(cnt,:)/CVmaxk(cnt,1),'b')
    		%loglog(FCk(cnt,:),CVexactk(cnt,:)/CVcon,'b')
    		hold on
    		%waitforbuttonpress		
    		%pause(5)
	    end
    end
    end
    title("P="+P)
end
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%function to cluster data based on noise monotonicity%%%%%%%%%%%

function [FCk_m,CVexactk_m,alarr_m,bearr_m,indc]=cluster_data(t_s,t_e,FCk,CVexactk,alarr,bearr,CVcon,col)
	cnt = 1;
	cnt2 = 0;
	x = logspace(log10(t_s),log10(t_e),1000);
	y1 = interp1(FCk(1,:),CVexactk(1,:),x);
	alc = alarr(1,1);
	blc = bearr(1,1);
	indc(cnt) = 1; % counter for the row selected

	for i=2:length(FCk(:,1))
		y = interp1(FCk(i,:),CVexactk(i,:),x);
		check = 1;
		for j=1:length(alc)
			y1 = interp1(FCk(indc(j),:),CVexactk(indc(j),:),x);
			ind1 = find(y>y1);
			ind2 = find(y<y1);
			if length(ind1)*length(ind2)>0
				check = 0;
			end
		end
	
		if check == 1
			cnt = cnt+1;
			indc(cnt) = i; % counter for the row selected
			alc(cnt,1) = alarr(i,1);
			blc(cnt,1) = bearr(i,1);
			
			subplot(1,2,1)
			loglog(x,y/CVcon',col)
			axis([t_s t_e 0.01 10])
			hold on
			subplot(1,2,2)
			loglog(alarr(i),bearr(i),'*','MarkerEdgeColor',col)
			hold on
			axis([0.001 1000 0.001 1000])
			%w = waitforbuttonpress;
			%pause(1)
		else
			cnt2 = cnt2+1;
			FCk_m(cnt2,:) = FCk(i,:);
			CVexactk_m(cnt2,:) = CVexactk(i,:);
			alarr_m(cnt2,:) = alarr(i);
			bearr_m(cnt2,:) = bearr(i);			
		end
	end
end

%%%%%%%%%%%%%%Function to compute noise %%%%%%%%%%%
function [M0 Mss FC FCmax CVexact CVmax CVcon] = Find_CV(al, be, r1, kpon, kpoff, kTFon, kTFoff, ga)

%  Analytical solution
	V2 = r1/kpoff;
	R = kTFon/kTFoff;
	P = kpon/kpoff;
	V1 = kTFoff/kpoff;
	r2 = al*r1;
	
	N1 = be*R*P*( V1+V2+1+P+R*V1);
	N2 = (1+P+V2)*(1+R)*(1+V1+al*be*V2) + (1+P+V2)*P*(1+be*R) + R*V1*(1+R)*( 1 + P*be +al*be*V2 ); 
	P11 = N1/N2;
	P01 = R/(1+R) - (1+be*R)*P11/(be*(1+R));
	P10 = P/(1+P+V2) - (1+be*(al*V2+P))*P11/(be*(1+P+V2));

	M0 = r1*P/((1+P+V2)*ga); % Constitutive expression

	FC = 1 + (al*be*(1+P)-P*be-1)*P11/(P*be);
	Mss = M0 + (r1/ga)*P11*(al - (1+be*(P+al*V2))/(be*(1+P+V2))); % mrna expression

% Noise

	A1 = kTFoff*(kTFon + kpon + ga)/be + kTFoff*(kTFoff + kpoff)/be + kTFoff*(r2+ga);
	A2 = r1*P10*kpon;
	B1 = (kTFon + ga)*(kpon+kpoff+kTFon+r1+ga) + kTFon*kTFoff;
	C1 = kpon + kpoff + r1 + ga;
	C2 = kpon + kpoff/be + r2 + ga;
	D1 = kpon*Mss;
	M11e = (B1*D1 - A2*C1)/(B1*C2+A1*C1); % Exact
	M10e = (A1*M11e + A2 )/B1;

	Fanoexact = 1 + ( r1*M10e + r2*M11e )/(ga*Mss) - Mss;
	CVexact = (Fanoexact/Mss)^0.5;
	
	CVcon = ((ga/r1 + kpon/(kpon + r1 + kpoff + ga))*(kpon + r1 + kpoff)/kpon-1)^0.5;
	CVmax = ((ga/r2 + kpon/(kpon + r2 + kpoff/be + ga))*(kpon + r2 + kpoff/be)/kpon-1)^0.5;
  	FCmax = al*be*(1+P+V2)/(1+P*be+al*be*V2);
	FNcon = CVcon^2*M0;
	FNmax = CVmax^2*FCmax*M0;


end

