

%% Method 1: Use the N2O Value From Park et al. 2004 (doi:10.1029/2003JD003731) and D17O from Boering 2004 (doi:10.1029/2003GL018451)
% 100 hPa = 248 ppb N2O = 1.57 per mil
% 380 K = 295 ppb N2O = 0.404 per mil 


% import netCDF of AgeTrop (mean age from tropopause to surface)
ncinfo("~/AgeTrop.hourly.nc")
ncname = "~/AgeTrop.hourly.nc";

% 33N = lat(32) 117W = lon(14)  [La Jolla]
age_array = ncread(ncname,"SpeciesConc_AgeTrop");
age_series = squeeze(age_array(14,32,1,:));

time = ncread(ncname,"time");

epoch = datetime(1980,01,01,00,00,00);
date = epoch + minutes(time);

figure
plot(date,age_series)
xlabel("Date")
ylabel("Mean Age (y)")
set(gca,"Box","On","FontSize",24)
title("Mean Age at 33\circN, 117\circW")


% isolate the 1990s

date_90s = date(date>=datetime(1990,01,01) & date < datetime(2000,01,01));
age_90s = age_series(date>=datetime(1990,01,01) & date < datetime(2000,01,01));

figure
plot(date,age_series)
xlabel("Date")
ylabel("Mean Age (y)")
set(gca,"Box","On","FontSize",24)
title("Mean Age at 33\circN, 117\circW")
hold on
plot(date_90s,age_90s)
hold off


% consolidate as monthly means

the90s = retime(timetable(date_90s,age_90s),'monthly','mean'); % monthly means
std90s = retime(timetable(date_90s,age_90s),'monthly',@std); % monthly standard deviations


the90s = join(the90s,std90s);
the90s.Properties.VariableNames = {'Mean', 'std'};

hold on
plot(the90s.date_90s,the90s.std)
hold off

figure
hold on
plot(date_90s,age_90s)
plot(the90s.date_90s,the90s.Mean)
hold off

% convert to days (from years)
the90s.Mean = the90s.Mean .* 365; % convert to daily
the90s.std = the90s.std .* 365;

figure
hold on
plot(date_90s,age_90s,"LineWidth",0.5,"Color",[200 200 200]./255)
plot(the90s.date_90s,the90s.Mean./365,"LineWidth",5)
hold off
set(gca,"FontSize",24,"Box","on")
xlabel("Date")
ylabel("Mean Age Strat -> Surface (y)")
title("GEOS-Chem Mean Age at 33\circN, 117\circW")
legend(["Hourly" "Monthly Average of Hourly"],"Location","southeast","box","off")

% calculate age spectra
t = 1:10000;
for i = 1:height(the90s)
    G_array(i,:) = G_diff(the90s.Mean(i),0.9 .* the90s.Mean(i),t); % factor of 0.9 suggested by fig 1 of Orbe et al. 2016 doi:10.1175/JAS-D-15-0289.1
end


figure
plot(t./365,G_array)
set(gca,"FontSize",24,"Box","on")
xlim([0 1])
xlabel("Transit Time (y)")
ylabel("Probability Density (d^{-1})")

% Calculate Surface Contributions
% name code mixsurf_FLUX_TAU, where A100 = 100 hPa and A380 = 380 K

for i = 1:height(G_array)
    mixsurf_H100_5(i) = sum(G_array(i,:) .* exp(-t./5) .* 1.57); % strat boundary condition = 1.57 permil
end
for i = 1:height(G_array)
    mixsurf_H100_7(i) = sum(G_array(i,:) .* exp(-t./7) .* 1.57); % strat boundary condition = 1.57 permil
end
for i = 1:height(G_array)
    mixsurf_H100_15(i) = sum(G_array(i,:) .* exp(-t./15) .* 1.57); % strat boundary condition = 1.57 permil
end
for i = 1:height(G_array)
    mixsurf_H100_30(i) = sum(G_array(i,:) .* exp(-t./30) .* 1.57); % strat boundary condition = 1.57 permil
end

for i = 1:height(G_array)
    mixsurf_A380_5(i) = sum(G_array(i,:) .* exp(-t./5) .* 0.404); % strat boundary condition = 0.404 permil
end
for i = 1:height(G_array)
    mixsurf_A380_7(i) = sum(G_array(i,:) .* exp(-t./7) .* 0.404); % strat boundary condition = 0.404 permil
end
for i = 1:height(G_array)
    mixsurf_A380_15(i) = sum(G_array(i,:) .* exp(-t./15) .* 0.404); % strat boundary condition = 0.404 permil
end
for i = 1:height(G_array)
    mixsurf_A380_30(i) = sum(G_array(i,:) .* exp(-t./30) .* 0.404); % strat boundary condition = 0.404 permil
end

figure
hold on
plot(the90s.date_90s,mixsurf_H100_5,"-","LineWidth",3,"Color",gca().ColorOrder(1,:))
plot(the90s.date_90s,mixsurf_H100_7,"-","LineWidth",3,"Color",gca().ColorOrder(2,:))
plot(the90s.date_90s,mixsurf_H100_15,"-","LineWidth",3,"Color",gca().ColorOrder(3,:))
plot(the90s.date_90s,mixsurf_H100_30,"-","LineWidth",3,"Color",gca().ColorOrder(4,:))
plot(the90s.date_90s,mixsurf_A380_5,"--","LineWidth",3,"Color",gca().ColorOrder(1,:))
plot(the90s.date_90s,mixsurf_A380_7,"--","LineWidth",3,"Color",gca().ColorOrder(2,:))
plot(the90s.date_90s,mixsurf_A380_15,"--","LineWidth",3,"Color",gca().ColorOrder(3,:))
plot(the90s.date_90s,mixsurf_A380_30,"--","LineWidth",3,"Color",gca().ColorOrder(4,:))
hold off
set(gca,"FontSize",24,"Box","on")
xlabel("Date")
ylabel("Stratospheric Contribution to Surface \Delta^{17}O-CO_2 (‰)")
legend(["\tau = 5 days" "\tau = 7 days" "\tau = 15 days" "\tau = 30 days"],"box","off","Location","southoutside","Orientation","horizontal")



%% Method 2: Apply a Seasonal Cycle to these Boundaries


%100 hPa
figure
plot(Hfluxes) %fluxes from Holton 1990 (doi:10.1175/1520-0469(1990)047<0392:OTGEOM>2.0.CO;2)

norm_Hfluxes = Hfluxes./mean(Hfluxes); %normalized by mean

figure
plot(norm_Hfluxes)

the90s.month = month(the90s.date_90s);

the90s_H100 = the90s;

for i = 1:height(the90s_H100)
    the90s_H100.Flux(i) = norm_Hfluxes(the90s_H100.month(i)) .* 1.57;
end

figure
plot(the90s_H100.date_90s,the90s_H100.Flux,"LineWidth",3)
set(gca,"FontSize",24,"Box","on")
xlabel("Date")
ylabel("Monthly 100 hPa \Delta^{17}O (‰)")

for i = 1:height(G_array)
    mixsurf2_H100_5(i) = sum(G_array(i,:) .* exp(-t./5) .* the90s_H100.Flux(i)); % strat boundary condition = 1.57 permil
end
for i = 1:height(G_array)
    mixsurf2_H100_7(i) = sum(G_array(i,:) .* exp(-t./7) .* the90s_H100.Flux(i)); % strat boundary condition = 1.57 permil
end
for i = 1:height(G_array)
    mixsurf2_H100_15(i) = sum(G_array(i,:) .* exp(-t./15) .* the90s_H100.Flux(i)); % strat boundary condition = 1.57 permil
end
for i = 1:height(G_array)
    mixsurf2_H100_30(i) = sum(G_array(i,:) .* exp(-t./30) .* the90s_H100.Flux(i)); % strat boundary condition = 1.57 permil
end

figure
hold on
plot(the90s.date_90s,mixsurf2_H100_5,"-","LineWidth",3,"Color",gca().ColorOrder(1,:))
plot(the90s.date_90s,mixsurf2_H100_7,"-","LineWidth",3,"Color",gca().ColorOrder(2,:))
plot(the90s.date_90s,mixsurf2_H100_15,"-","LineWidth",3,"Color",gca().ColorOrder(3,:))
plot(the90s.date_90s,mixsurf2_H100_30,"-","LineWidth",3,"Color",gca().ColorOrder(4,:))
hold off
set(gca,"FontSize",24,"Box","on")
xlabel("Date")
ylabel("Stratospheric Contribution to Surface \Delta^{17}O-CO_2 (‰)")
legend(["\tau = 5 days" "\tau = 7 days" "\tau = 15 days" "\tau = 30 days"],"box","off","Location","southoutside","Orientation","horizontal")

%380 K
figure
plot(appfluxes380) %fluxes from Appenzeller 1996 (doi:10.1029/96JD00821)

norm_appfluxes380 = appfluxes380./mean(appfluxes380); %normalized by mean

figure
plot(norm_appfluxes380)

the90s.month = month(the90s.date_90s);

the90s_A380 = the90s;

for i = 1:height(the90s_A380)
    the90s_A380.Flux(i) = norm_appfluxes380(the90s_A380.month(i)) .* 0.404;
end

figure
plot(the90s_A380.date_90s,the90s_A380.Flux,"LineWidth",3)
set(gca,"FontSize",24,"Box","on")
xlabel("Date")
ylabel("Monthly 380 K \Delta^{17}O (‰)")

for i = 1:height(G_array)
    mixsurf2_A380_5(i) = sum(G_array(i,:) .* exp(-t./5) .* the90s_A380.Flux(i)); % strat boundary condition = 1.57 permil
end
for i = 1:height(G_array)
    mixsurf2_A380_7(i) = sum(G_array(i,:) .* exp(-t./7) .* the90s_A380.Flux(i)); % strat boundary condition = 1.57 permil
end
for i = 1:height(G_array)
    mixsurf2_A380_15(i) = sum(G_array(i,:) .* exp(-t./15) .* the90s_A380.Flux(i)); % strat boundary condition = 1.57 permil
end
for i = 1:height(G_array)
    mixsurf2_A380_30(i) = sum(G_array(i,:) .* exp(-t./30) .* the90s_A380.Flux(i)); % strat boundary condition = 1.57 permil
end

figure
hold on
plot(the90s.date_90s,mixsurf2_A380_5,"-","LineWidth",3,"Color",gca().ColorOrder(1,:))
plot(the90s.date_90s,mixsurf2_A380_7,"-","LineWidth",3,"Color",gca().ColorOrder(2,:))
plot(the90s.date_90s,mixsurf2_A380_15,"-","LineWidth",3,"Color",gca().ColorOrder(3,:))
plot(the90s.date_90s,mixsurf2_A380_30,"-","LineWidth",3,"Color",gca().ColorOrder(4,:))
hold off
set(gca,"FontSize",24,"Box","on")
xlabel("Date")
ylabel("Stratospheric Contribution to Surface \Delta^{17}O-CO_2 (‰)")
legend(["\tau = 5 days" "\tau = 7 days" "\tau = 15 days" "\tau = 30 days"],"box","off","Location","southoutside","Orientation","horizontal")


save La_Jolla_Method_2_results mixsurf2_H100_5 mixsurf2_H100_7 mixsurf2_H100_15 mixsurf2_H100_30 mixsurf2_A380_5 mixsurf2_A380_7 mixsurf2_A380_15 mixsurf2_A380_30


% Combined results

figure
hold on
plot(the90s.date_90s,mixsurf2_H100_5,"-","LineWidth",3,"Color",gca().ColorOrder(1,:))
plot(the90s.date_90s,mixsurf2_H100_7,"-","LineWidth",3,"Color",gca().ColorOrder(2,:))
plot(the90s.date_90s,mixsurf2_H100_15,"-","LineWidth",3,"Color",gca().ColorOrder(3,:))
plot(the90s.date_90s,mixsurf2_H100_30,"-","LineWidth",3,"Color",gca().ColorOrder(4,:))
plot(the90s.date_90s,mixsurf2_A380_5,"--","LineWidth",3,"Color",gca().ColorOrder(1,:))
plot(the90s.date_90s,mixsurf2_A380_7,"--","LineWidth",3,"Color",gca().ColorOrder(2,:))
plot(the90s.date_90s,mixsurf2_A380_15,"--","LineWidth",3,"Color",gca().ColorOrder(3,:))
plot(the90s.date_90s,mixsurf2_A380_30,"--","LineWidth",3,"Color",gca().ColorOrder(4,:))
hold off
set(gca,"FontSize",24,"Box","on")
xlabel("Date")
ylabel("Stratospheric Contribution to Surface \Delta^{17}O-CO_2 (‰)")
legend(["\tau = 5 days" "\tau = 7 days" "\tau = 15 days" "\tau = 30 days"],"box","off","Location","southoutside","Orientation","horizontal")




%% Full Combined Results


figure
p = tiledlayout(2,1,"TileSpacing","compact");

nexttile
hold on
plot(the90s.date_90s,mixsurf_H100_5,"-","LineWidth",3,"Color",gca().ColorOrder(1,:))
plot(the90s.date_90s,mixsurf_H100_7,"-","LineWidth",3,"Color",gca().ColorOrder(2,:))
plot(the90s.date_90s,mixsurf_H100_15,"-","LineWidth",3,"Color",gca().ColorOrder(3,:))
plot(the90s.date_90s,mixsurf_H100_30,"-","LineWidth",3,"Color",gca().ColorOrder(4,:))
plot(the90s.date_90s,mixsurf_A380_5,"--","LineWidth",3,"Color",gca().ColorOrder(1,:))
plot(the90s.date_90s,mixsurf_A380_7,"--","LineWidth",3,"Color",gca().ColorOrder(2,:))
plot(the90s.date_90s,mixsurf_A380_15,"--","LineWidth",3,"Color",gca().ColorOrder(3,:))
plot(the90s.date_90s,mixsurf_A380_30,"--","LineWidth",3,"Color",gca().ColorOrder(4,:))
hold off
set(gca,"FontSize",24,"Box","on")
title("Method 1")

nexttile
hold on
plot(the90s.date_90s,mixsurf2_H100_5,"-","LineWidth",3,"Color",gca().ColorOrder(1,:))
plot(the90s.date_90s,mixsurf2_H100_7,"-","LineWidth",3,"Color",gca().ColorOrder(2,:))
plot(the90s.date_90s,mixsurf2_H100_15,"-","LineWidth",3,"Color",gca().ColorOrder(3,:))
plot(the90s.date_90s,mixsurf2_H100_30,"-","LineWidth",3,"Color",gca().ColorOrder(4,:))
plot(the90s.date_90s,mixsurf2_A380_5,"--","LineWidth",3,"Color",gca().ColorOrder(1,:))
plot(the90s.date_90s,mixsurf2_A380_7,"--","LineWidth",3,"Color",gca().ColorOrder(2,:))
plot(the90s.date_90s,mixsurf2_A380_15,"--","LineWidth",3,"Color",gca().ColorOrder(3,:))
plot(the90s.date_90s,mixsurf2_A380_30,"--","LineWidth",3,"Color",gca().ColorOrder(4,:))
hold off
set(gca,"FontSize",24,"Box","on")
title("Method 2")
legend(["\tau = 5 days" "\tau = 7 days" "\tau = 15 days" "\tau = 30 days"],"box","off","Location","southoutside","Orientation","horizontal")

xlabel(p,"Date","FontSize",36)
ylabel(p,"Stratospheric Contribution to Surface \Delta^{17}O-CO_2 (‰)","FontSize",36)


