%-----Reading data from CAS, 2DS, PCASP, ORACLES data files-----%
%-----Written by Siddhant Gupta, Univ. of Oklahoma - 11/20/2019-----%

%-----Input files------%
% 1. CAS: Droplet size distribution (0.5 < D < 50 um), we advise using data from the 3-50 um size range
%       e.g. 16_09_06_06_40_17.Combined.conc.casn.sync.1Hz.nc
% 2. 2DS: Droplet size distribution (10 < D < 1280 um), we strongly advise using data for D > 50 um only
%       e.g. base160906064017.combined.H.2DS.cdf
% 3. PCASP: Aerosol size distribution (0.1 < D < 3 um), we advise excluding data from the first bin
%       e.g. 16_09_06_06_40_17.Combined.conc.spp200.sync.1Hz.nc
% 3. ORACLES: Bulk LWC from King hot-wire and basic variables (lat,lon,aircraft altitude, etc.)
%       e.g. 16_09_06_06_40_17.Combined.oracles.sync.nc

%-----Output parameters-----%
%       1. Droplet size distribution (Nd)
%       2. Droplet concentration (CAS: cas_Nt, 2DS: twods_Nt, total: Nt)
%       3. Aerosol concentration (PCASP: pcasp_Nt)
%       4. Effective radius (re) - cloud mode/precipitation mode
%       5. Liquid Water Content (lwc) - cloud mode/precipitation mode
%       6. Liquid Water Content (king_lwc) - bulk measurement from hot-wire
%       7. Rain rate (r)

%-----Identify/retrieve data files-----%
for j=1
k1='06';k2='10';k3='12';k4='14';k5='20';k6='25';
flightdate=['201609' eval(sprintf('k%d',j))];
ncid_c=dir(fullfile(pwd,['18_09_' eval(sprintf('k%d',j)) '*casn.*']));
ncid_s=dir(fullfile(pwd,['18_09_' eval(sprintf('k%d',j)) '*oracles.*']));
ncid_p=dir(fullfile(pwd,['18_09_' eval(sprintf('k%d',j)) '*.spp200*z.nc']));
ncid_t=dir(fullfile(pwd,['*1809' eval(sprintf('k%d',j)) '*H.2DS*']));
ncid=netcdf.open(ncid_c.name);
a=netcdf.open(ncid_p.name);
b=netcdf.open(ncid_s.name);
c=netcdf.open(ncid_t.name);

altitude=ncread(ncid_s.name,'Sys_GPSAlt'); % For 2018
king_lwc=ncread(ncid_s.name,'DLWCC'); % Bulk LWC from King hot-wire
end

%-----Read Unix time and convert to Time Vector for CAS file-----%
base_time_cas=netcdf.getVar(ncid,7);
time_offset_cas=netcdf.getVar(ncid,8);
offset_cas=[0:numel(time_offset_cas)-1]';
base_seconds_cas=offset_cas+double(base_time_cas);
timevec_cas=unix_to_timevec(base_seconds_cas);
%-----Read Unix time and convert to Time Vector for ORACLES file-----%
base_time=netcdf.getVar(b,7);
time_offset=netcdf.getVar(b,8);
offset=[0:numel(time_offset)-1]';
base_seconds=offset+double(base_time);
timevec=unix_to_timevec(base_seconds);
%-----Read Unix time and convert to Time Vector for PCASP file-----%
base_time_pcasp=netcdf.getVar(a,7);
time_offset_pcasp=netcdf.getVar(a,8);
offset_pcasp=[0:numel(time_offset_pcasp)-1]';
base_seconds_pcasp=offset_pcasp+double(base_time_pcasp);
timevec_pcasp=unix_to_timevec(base_seconds_pcasp);

%-----Read variables from 2DS file-----%
timevec_twods=netcdf.getVar(c,0);
bin_min_twods=netcdf.getVar(c,1)*1000; % Min D for each bin (converted from mm to um). 
bin_mid_twods=netcdf.getVar(c,3)*1000; % Max D for each bin (converted from mm to um). 
bin_dD_twods=netcdf.getVar(c,4)*1000;  % Bin width for each bin (converted from mm to um). 
twods_Nd=netcdf.getVar(c,7);           % 2DS size distribution
twods_Nt=sum(twods_Nd(:,5:128),2);     % 2DS droplet concentration (D > 50 um used). 

%-----Read CAS droplet concentration-----%
cas_Nt=netcdf.getVar(ncid,41); % 3 to 50 um

%-----Read PCASP aerosol concentration-----%
pcasp_Nt=netcdf.getVar(b,39); % excluding first bin

%-----Indices to plot data from different files together (files can have different variable lengths)-----%
x=min(length(cas_Nt),length(twods_Nt));
y=min(length(cas_Nt),length(king_lwc));

%-----Flag to remove noise (CAS N < 10 cm-3 or hot-wire LWC < 0.05 g m-3)-----%
king_lwc_limit=0.05; cas_limit=10;
index1=find(cas_Nt(1:y)<cas_limit | king_lwc(1:y)<king_lwc_limit);
cas_Nt(index1)=NaN; twods_Nt(index1)=NaN; king_lwc(index1)=NaN;

%-----Add CAS and 2DS droplet concentrations-----%
Nt=cat(2,cas_Nt(1:x),twods_Nt(1:x));
Nt=nansum(Nt,2);

%-----Flags for effective radius (Re) or liquid water content (LWC) calculation-----%
cas=1; % 0 to exclude CAS data to get precipitation mode Re/LWC
twods=1; % 0 to exclude 2DS data to get cloud mode Re/LWC

%-----Start radius for R-calculation-----% 
r_radius=17; % 1 for 3um,8 for 10um,12 for 25um,13 for 30um,15 for 40um,17 for 50um
%% %-----Calculate Effective Radius from CAS and 2DS size distribution-----%

%-----Retrieve CAS size distribution from CAS file-----%
channels=16;
cas_Nd=zeros(numel(timevec_cas),channels);
for k=23:38
    eval(sprintf('cas_Nd(:,%d-22)=netcdf.getVar(ncid,%d);',k,k));
end

binmid=[3.25;3.75;4.5;5.75;6.85;7.55;9.05;11.35;13.75;17.5;22.5;27.5;32.5;37.5;42.5;47.5];

%-----Remove CAS and 2DS size distribution for noise-----%
cas_Nd(index1,:)=NaN;
twods_Nd(index1,:)=NaN;
if cas==0
cas_Nd(:,:)=NaN;
else
end
if twods==0
twods_Nd(:,:)=NaN;
else
end

Nd=cat(2,cas_Nd(1:x,:),twods_Nd(1:x,5:128));
bin_mid=cat(1,binmid,bin_mid_twods(5:128));

original_binmid_all=repmat(bin_mid,1,x);

%-----Effective Radius calculation-----%
num_all=zeros(numel(bin_mid),x);
den_all=zeros(numel(bin_mid),x);

for k=1:x
for i=1:numel(bin_mid)
    num_all(i,k)=Nd(k,i).*((original_binmid_all(i,k)).^3);
    den_all(i,k)=Nd(k,i).*((original_binmid_all(i,k)).^2);
end
end

numsum_all=nansum(num_all);
densum_all=nansum(den_all);
frac_all=numsum_all./densum_all;
re=(frac_all/2)';
%% %-----Calculate LWC from CAS and 2DS size distributions-----%

%-----CAS LWC calculation-----%
cas_lwc_bin=zeros(numel(timevec_cas),numel(binmid));
binmid3=power(binmid,3);

for k=1:numel(timevec_cas)
    for i=1:numel(binmid)
    cas_lwc_bin(k,i)=pi/6.*cas_Nd(k,i).*binmid3(i);
    end
end
cas_lwc=nansum(cas_lwc_bin,2)./1000000; % Convert to g / m3

%-----2DS LWC calculation-----%
twods_lwc_bin=zeros(numel(timevec_twods),numel(bin_mid_twods(5:128)));
bin_mid3=power(bin_mid_twods(1:128),3);

for k=1:numel(timevec_twods)
    for i=5:numel(bin_mid_twods)
    twods_lwc_bin(k,i)=pi/6.*twods_Nd(k,i).*bin_mid3(i);
    end
end
twods_lwc=nansum(twods_lwc_bin,2)./1000000;  % Convert to g / m3

%-----Add CAS and 2DS liquid water content-----%
lwc=cas_lwc(1:x)+twods_lwc(1:x);
lwc(lwc==0)=NaN;
%% %-----Calculate rain Rate from CAS and 2DS size distributions-----%

original_binmid=power(original_binmid_all/10000,3); %D^3, cm3

%-----Fall Speed/Rain Rate calculation based on Droplet fall speed relationships based on droplet radii-----%
k1=1.19*power(10,6); k2=8*power(10,3); k3=2.01*power(10,3); % Rogers (1976)
u1=k1.*power(bin_mid(1:14)/10000,2); % D<40um, u1=k1r^2
u2=k2.*bin_mid(15:72)/10000; % 40<D<600um, u2=k2r
u3=k3.*power(bin_mid(73:140)/10000,0.5); % D>600um, u3=k3r^1/2
u=vertcat(u1,u2,u3);
u_all=repmat(u,1,x); % cm/s
r_all=zeros(numel(u),x);

for k=1:x
for i=r_radius:numel(bin_mid)
    r_all(i,k)=Nd(k,i).*u_all(i,k).*original_binmid(i,k);
end
end
r=nansum(r_all)*36000; %cm/s to mm/hr
%%

clearvars -except timevec* r re lwc *Nt Nd king_lwc bin_mid