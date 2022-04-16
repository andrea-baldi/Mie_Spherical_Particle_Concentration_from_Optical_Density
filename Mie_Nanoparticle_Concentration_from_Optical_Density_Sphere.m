% Baldi Lab 16/04/2022

% This script reads the relative permittivity data of a material and uses
% Mie theory to compute the extinction cross-section for a spherical particle
% of that material embedded in a lossless medium. The permittivity file has
% to be a tab-delimited text file with three columns: energy (in eV),
% epsilon1, epsilon2.
% 
% From the particle size, the script computes the particle concentration
% for a solution with a specified optical density at resonance.
clear all
close all

% MANUAL INPUT
prompt = {'Sphere radius (nm)','Refractive index of the surrounding medium (air = 1, water = 1.333)','Measured optical density at resonance','Optical path of the cuvette (cm)'};
dlg_title = 'Input parameters';
num_lines = 1;
def = {'','','',''};
answer = inputdlg(prompt,dlg_title,num_lines,def);
r=str2double(answer{1});
index=str2double(answer{2});
OD = str2double(answer{3}); % Optical density (= Absorbance = -Log(Transmittance) = Log(I0/I) ) at the wavelength of maximum extinction
optical_path = str2double(answer{4});

% DEFAULT SETTINGS
Esteps = 200;        % Number of energy steps
nmax = 10;           % Maximum order of the Bessel function

% CALCULATIONS
read=dlmread(uigetfile('*','Select a dielectric function file')); % Load a dielectric function file (energy in eV, eps1, eps2)
e = 1.60217646e-19;    % Elementary charge in SI units
h = 6.626068e-34;      % Planck's constant in SI units
hbar = 1.05457148e-34; % hbar in SI units
me = 9.10938215e-31;   % Electron rest mass in SI units
c = 2.99792458e8;      % Light speed in SI units

% Load Experimental data
eV_read=read(:,1);
e1_read=read(:,2);
e2_read=read(:,3);

% Interpolates the tabulated experimental data to make the data set smooth
Emin=read(1,1);
Emax=read(size(read,1),1);
energy = (Emin:(Emax-Emin)/Esteps:Emax)';
e1_read = interp1(eV_read, e1_read, energy, 'spline');
e2_read = interp1(eV_read, e2_read, energy, 'spline');

% Converts the frequency data to wavelength in meters
lambda = h*c./(e*energy);

% Creates the wavenumber k
k = 2*pi*index./lambda;

% Creates the total permittivity values of the particle and of the
% lossless medium
n_read = (((e1_read.^2 + e2_read.^2).^(1/2) + e1_read)./2).^(1/2);
k_read = (((e1_read.^2 + e2_read.^2).^(1/2) - e1_read)./2).^(1/2);
etot_read = n_read + 1i*k_read;
emed= (index^2)*ones(size(etot_read));
m=etot_read./index;

radius = r*1e-9; % Convertion to meters
x = k.*radius; % Size parameter
mx = m.*x;
scaele = 0; % initialise scattering matrix element
extele = 0; % initialise extinction matrix element
for n=1:nmax
    jnx = sqrt(pi./(2.*x)).*besselj(n+0.5,x);
    jnminx = sqrt(pi./(2.*x)).*besselj(n-0.5,x);
    jnmx = sqrt(pi./(2.*mx)).*besselj(n+0.5,mx);
    jnminmx = sqrt(pi./(2.*mx)).*besselj(n-0.5,mx);
    hnx = sqrt(pi./(2.*x)).*besselh(n+0.5,x);
    hnminx = sqrt(pi./(2.*x)).*besselh(n-0.5,x);
    xjnxdiff = x.*jnminx - n.*jnx;
    mxjnmxdiff = mx.*jnminmx - n.*jnmx;
    xhnxdiff = x.*hnminx - n.*hnx;
    an = (m.^2.*jnmx.*xjnxdiff-jnx.*mxjnmxdiff)./(m.^2.*jnmx.*xhnxdiff-hnx.*mxjnmxdiff);
    bn = (jnmx.*xjnxdiff-jnx.*mxjnmxdiff)./(jnmx.*xhnxdiff-hnx.*mxjnmxdiff);
    scaele = scaele + (2*n+1).*(an.*conj(an)+bn.*conj(bn));
    extele = extele + (2*n+1).*real(an+bn);
end

% Calculate scattering, extinciton, and absorption cross sections and efficiencies
% Csca = 2*pi./(k.^2).*scaele;
Cext = 2*pi./(k.^2).*extele;
% Cabs = Cext-Csca;
% Qsca = Csca./(pi*radius^2);
% Qext = Cext./(pi*radius^2);
% Qabs = Cabs./(pi*radius^2);

% Plot the extinction cross-section and choose the energy value corresponding to your optical density
screen_size = get(0,'ScreenSize');
f1 = figure(1);
plot(energy,Cext,'k','linewidth',2)
xlabel('Energy (eV)', 'FontSize', 10 );
ylabel(['Extinction cross-section of a ',num2str(radius*1E9),' nm radius sphere'],'FontSize',10);
title('Select the energy value corresponding to your optical density');
axis([Emin Emax 0 1.1*max(Cext)])
[x0,y0]=ginput(1);
[dummyvalue,indexplot]=min(abs(energy-x0));
psdensity=log(10)*OD/(Cext(indexplot)*optical_path/100); %Calculate the density in particles/m^3 using Beer-Lambert Law. The factor log(10) comes from the fact that most spectrometers plot the absorbance as the 10-base logarithm of I0/I, while the Beer-Lambert law expressed in terms of the extinction cross section, Cext, uses the natural logarithm.
psdensityml=psdensity./1E6 % density in particles/ml
txt = ['density = ',num2str(psdensityml,'%g'),' particles/ml'];
text(Emin*1.2,Cext(indexplot)*1.05,txt,'FontSize',12)
