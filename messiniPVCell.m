clear
clc
Voc = 46.22;
Isc = 9.06;
Ns = 72;
Np = 1;
A = 1.3;
Ki = 0.005258;
Gref = 1000;
q = 1.602*1e-19;
Rs = 0.008;
Eg = 1.1;
K = 1.3805e-23;
Trtc = 25; 
To = Trtc + 273.15;

V = linspace(0,100,1000); % i
V = V';
Trc = 10:15:70; % j
G = linspace(200,1000,5); % k
I = zeros(length(V),length(Trc),length(G)); 

for i = 1:length(V)
    for j = 1:length(Trc)
        for k = 1:length(G) 
            Tr(j) = Trc(j)+273.15;
            NskATo = Ns*K*A*To;
            Iph = (G(k)*(Isc+((To-Tr(j))*Ki)))/Gref;
            Irs = Isc/(exp(Voc*q/NskATo)-1);
            Is = Irs*((Tr(j)/To)^3)*exp((q*Eg/(K*A))*(1/Tr(j)-1/To));
            Err = 1;
            tolErr = 1e-5;
            while abs(Err)>tolErr
                Inew = Np*Iph-Np*Is*(exp(q*(V(i)-Rs*I(i,j,k))/NskATo)-1);
                Err = Inew-I(i,j,k);
                I(i,j,k) = Inew;
            end
        end
    end
end

figure(1)
P_Tconst = V.*ones(size(squeeze(I(:,2,:)))).*squeeze(I(:,2,:));
plot(V,squeeze(I(:,2,:)),'LineWidth',2),grid
ylim([0 inf])
xlabel('Voltage(V)')
ylabel('Current(A)')
legend('200 W/m^2','400 W/m^2','600 W/m^2','800 W/m^2','1000 W/m^2')
title('I–V characteristics, varying irradiance at constant temperature. = 25 °C')


figure(2)
P_Tconst = V.*ones(size(squeeze(I(:,2,:)))).*squeeze(I(:,2,:));
plot(V,P_Tconst,'LineWidth',2),grid
ylim([0 inf])
xlabel('Voltage(V)')
ylabel('Power(W)')
legend('200 W/m^2','400 W/m^2','600 W/m^2','800 W/m^2','1000 W/m^2')
title('P–V characteristics, varying solar radiation at constant temperature = 25 °C')

figure(3)
% P_Tconst = V.*ones(size(squeeze(I(:,:,end)))).*squeeze(I(:,:,end));
plot(V,squeeze(I(:,:,end)),'LineWidth',2),grid
ylim([0 inf])
xlabel('Voltage(V)')
ylabel('Current (A)')
legend('10 °C','25 °C','40 °C','55 °C','70 °C')
title('I–V characteristics, varying temperature at constant solar radiation = 1000 W/m^2')

figure(4)
P_Tconst = V.*ones(size(squeeze(I(:,:,end)))).*squeeze(I(:,:,end));
plot(V,P_Tconst,'LineWidth',2),grid
ylim([0 inf])
xlabel('Voltage(V)')
ylabel('Power(W)')
legend('10 °C','25 °C','40 °C','55 °C','70 °C')
title('I–V characteristics, varying temperature at constant solar radiation = 1000 W/m^2')


G_months = [99.52	160.92	203.17	236.31	214.99	250.65	203.28	193.55	215.01	181.58	122.02	87.22];
T_months = [11.45	17.09	20.22	26.85	32.68	30.08	29.96	28.92	28.79	25.99	20.97	14.84];

clear I Tr
I = zeros(length(V),length(G_months)); 
for i = 1:length(V)
    %for j = 1:length(T_months)
        for k = 1:length(G_months) 
            Tr(k) = T_months(k)+273.15;
            NskATo = Ns*K*A*To;
            Iph = (G_months(k)*(Isc+((To-Tr(k))*Ki)))/Gref;
            Irs = Isc/(exp(Voc*q/NskATo)-1);
            Is = Irs*((Tr(k)/To)^3)*exp((q*Eg/(K*A))*(1/Tr(k)-1/To));
            Err = 1;
            tolErr = 1e-5;
            while abs(Err)>tolErr
                Inew = Np*Iph-Np*Is*(exp(q*(V(i)-Rs*I(i,k))/NskATo)-1);
                Err = Inew-I(i,k);
                I(i,k) = Inew;
            end
        end
    %end
end

P_months = V.*ones(size(I)).*I;

figure(5)

subplot 321
plot(V,I(:,1:4),'LineWidth',2),grid
ylim([0 3])
xlabel('Voltage(V)')
ylabel('Current(A)')
legend('Jan','Feb','Mar','Apr')
title('I–V graph for the month of Jan-15 to April-15.')
%P_Tconst = V.*ones(size(squeeze(I(:,2,:)))).*squeeze(I(:,2,:));
subplot 322
plot(V,P_months(:,1:4),'LineWidth',2),grid
ylim([0 80])
xlabel('Voltage(V)')
ylabel('Power(W)')
legend('Jan','Feb','Mar','Apr')
title('P–V graph for the month of Jan-15 to April-15.')

subplot 323
plot(V,I(:,5:8),'LineWidth',2),grid
ylim([0 3])
xlabel('Voltage(V)')
ylabel('Current(A)')
legend('May','June','July','Aug')
title('I–V graph for the month of May-15 to August-15')
%P_Tconst = V.*ones(size(squeeze(I(:,2,:)))).*squeeze(I(:,2,:));
subplot 324
plot(V,P_months(:,5:8),'LineWidth',2),grid
ylim([0 80])
xlabel('Voltage(V)')
ylabel('Power(W)')
legend('May','June','July','Aug')
title('P–V graph for the month of May-15 to August-15')

subplot 325
plot(V,I(:,9:12),'LineWidth',2),grid
ylim([0 3])
xlabel('Voltage(V)')
ylabel('Current(A)')
legend('Sep','Oct','Nov','Dec')
title('I–V graph for the month of September-15 to December-15')
%P_Tconst = V.*ones(size(squeeze(I(:,2,:)))).*squeeze(I(:,2,:));
subplot 326
plot(V,P_months(:,9:12),'LineWidth',2),grid
ylim([0 80])
xlabel('Voltage(V)')
ylabel('Power(W)')
legend('May','June','July','Aug')
title('P–V graph for the month of September-15 to December-15')

% 
