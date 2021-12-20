function [rho_star,a_star] = cea(pcc,rof,d_t)

%% Eingansdaten
% Treibstoffe und ROF
fuel = 'H2';
oxid = 'O2(L)';
% rof = 6.0; % [-]
% Kontraktionsverhältnis und Düsenhals
% d_t = 28.0 * 0.001; % [m]
A_t = pi() / 4 * d_t^2; % [m^2]
acat = (pi() / 4 * 0.05^2) / A_t; % [-]
% Brennkammerdruck
pcc = pcc*1E-5; % [bar]

%% Erstellen des Inputfiles cea.inp
ceain(fuel,oxid,rof,pcc,acat);

%% Ausführen von cea2.exe (Gordon-McBride)
dos cea2.exe;
% ! cea2.exe;

%% Auslesen von cea.out
fid = fopen('cea.out');
C = textscan(fid,'%s');

% Drücke
% Infinite, End of Combustion, Hals
index_p_inf = strmatch('P,',C{1,1},'exact') + 2;
p_inf = str2double(C{1,1}(index_p_inf,1)) * 1E5; % [Pa]
index_p_EoC = index_p_inf + 1;
p_EoC = str2double(C{1,1}(index_p_EoC,1)) * 1E5; % [Pa]
index_p_star = index_p_EoC + 1;
p_star = str2double(C{1,1}(index_p_star,1)) * 1E5; % [Pa]

% Temperaturen
% Infinite, End of Combustion, Hals
index_T_inf = strmatch('T,',C{1,1},'exact') + 2;
T_inf = str2double(C{1,1}(index_T_inf,1)); % [K]
index_T_EoC = index_T_inf + 1;
T_EoC = str2double(C{1,1}(index_T_EoC,1)); % [K]
index_T_star = index_T_EoC + 1;
T_star = str2double(C{1,1}(index_T_star,1)); % [K]

% Dichten
% Infinite, End of Combustion, Hals
index_rho_inf = strmatch('RHO',C{1,1}) + 3;
rho_inf_char = char(C{1,1}(index_rho_inf));
if length(rho_inf_char) > 6
    rho_inf_exp = str2double(rho_inf_char(8));
    rho_inf_str = [rho_inf_char(1:5)];
    rho_inf = str2double(rho_inf_str) / 10^rho_inf_exp; % [kg/m^3]
    index_rho_EoC = index_rho_inf + 1;
else
    rho_inf = str2double(rho_inf_char); % [kg/m^3]
    index_rho_EoC = index_rho_inf + 2;
end
rho_EoC_char = char(C{1,1}(index_rho_EoC));
if length(rho_EoC_char) > 6
    rho_EoC_exp = str2double(rho_EoC_char(8));
    rho_EoC_str = [rho_EoC_char(1:5)];
    rho_EoC = str2double(rho_EoC_str) / 10^rho_EoC_exp; % [kg/m^3]
    index_rho_star = index_rho_EoC + 1;
else
    rho_EoC = str2double(rho_EoC_char); % [kg/m^3]
    index_rho_star = index_rho_EoC + 2;
end
rho_star_char = char(C{1,1}(index_rho_star));
if length(rho_star_char) > 6
    rho_star_exp = str2double(rho_star_char(8));
    rho_star_str = [rho_star_char(1:5)];
    rho_star = str2double(rho_star_str) / 10^rho_star_exp; % [kg/m^3]
else
    rho_star = str2double(rho_star_char); % [kg/m^3]
end

% Schallgeschwindigkeiten
% Infinite, End of Combustion, Hals
index_a_inf = strmatch('SON',C{1,1},'exact') + 2;
a_inf = str2double(C{1,1}(index_a_inf,1)); % [m/s]
index_a_EoC = index_a_inf + 1;
a_EoC = str2double(C{1,1}(index_a_EoC,1)); % [m/s]
index_a_star = index_a_EoC + 1;
a_star = str2double(C{1,1}(index_a_star,1)); % [Km/s]
