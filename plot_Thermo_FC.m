% Plot the Fold-change heatmap from thermodynamic model of simple
% regulation for finite R, and saturating R, fintie P and infinite P

R = 1;
P=0.1;

ala = logspace(log10(0.001),log10(1000),1000);
bea = logspace(log10(0.001),log10(1000),1000);
for j=1:length(ala)
for i=1:length(bea)
	FC1(i,j) = (1+P+ala(j)*bea(i)*R*(1+P))/(1+P+(1+P*bea(i))*R);
	FC2(i,j) = (1+ala(j)*bea(i)*R)/(1+bea(i)*R); %P=infinity, finite R

	FCM1(i,j) = ala(j)*bea(i)*(1+P)/(1+P*bea(i)); %R=infinity, finite P
	FCM2(i,j) = ala(j); % P and R=infinity
end
end
subplot(2,2,1)
imagesc(log10(FC1))
axis xy
title("P="+P+", R="+R)

subplot(2,2,2)
imagesc(log10(FC2))
axis xy
title("Infinite P and R="+R)

subplot(2,2,3)
imagesc(log10(FCM1))
axis xy
title("P="+P+", saturating R")

subplot(2,2,4)
imagesc(log10(FCM2))
axis xy
title("Both P and R saturating")

