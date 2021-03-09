% plots in 1D theoretical resolution vs aperture length
% plots in 1D nearest ambiguity location vs platform spacing
% plots in 2D min # of platforms required vs resolution and nearest ambiguity location
close all; clear; clc
%% parameters
c=3e8; % speed of light (m/s)
f=1.2e9; % center frequency (Hz)
wl=c/f; % wavelength (m)
H=700e3; % altitude (m)
theta=30; % look angle wrt nadir (deg)
alpha=30; % baseline tilt wrt horizontal (deg)
L=10e3:10:100e3; % aperture lengths (m)
dL=100:1:1000; % platform spacing (m)
res=2:0.01:30; % desired resolution (m)
d_amb=20:0.05:100; % desired nearest ambiguity location (m)
%% theoretical resolution vs aperture length and nearest ambiguity location vs platform spacing
p_res=[2 1 1.4]; %SAR SIMO MIMO constants for resolution
p_amb=[2 1 1]; %SAR SIMO MIMO constants for nearest ambiguity location
R=H/cosd(theta);
Ln=L*cosd(abs(theta-alpha));
figure(1);hold on
figure(2);hold on
for i=1:3
    dn=wl*R/p_res(i)./Ln; % 4dB resolution along perpendicular to look angle direction
    n_amb=wl*R./(p_amb(i)*dL*cosd(abs(theta-alpha))); % closest ambiguity location along perpendicular to look angle direction
    figure(1);plot(L/1e3,dn,'linewidth',2)
    if i==1;figure(2);plot(dL/1e3,n_amb/1e3,'linewidth',2);end
    if i==2;figure(2);plot(dL/1e3,n_amb/1e3,'linewidth',2);end
    if i==3;figure(2);plot(dL/1e3,n_amb/1e3,'k--','linewidth',2.5);end
end
figure(1)
title('Resolution vs Aperture Length')
xlabel('Aperture Length (along baseline) (km)')
ylabel('4 dB Resolution (along n) (m)')
legend('SAR','SIMO','MIMO')
set(gca,'fontsize',10)
xlim([min(L/1e3) max(L/1e3)])
figure(2)
title('Nearest Ambiguity Location vs Platform Spacing')
xlabel('Platform Spacing (along baseline) (km)')
ylabel('Nearest Ambiguity Location (along n) (km)')
legend('SAR','SIMO','MIMO')
set(gca,'fontsize',10)
xlim([min(dL/1e3) max(dL/1e3)])
%% min # of platform required for a specific resolution and nearest ambiguity location
n_platform=zeros(length(res),length(d_amb));
for mode=1:3 %1: SAR, 2:SIMO, 3:MIMO
    if mode==1;modd='SAR';
    elseif mode==2;modd='SIMO';
    elseif mode==3;modd='MIMO';end
    for i=1:length(res)
        L=wl*R/p_res(mode)./res(i)/cosd(abs(theta-alpha));
        dL=wl*H./(cosd(theta)*p_amb(mode)*d_amb*cosd(abs(theta-alpha)));
        n_platform(i,:)=ceil(L./dL);
    end
    max(max(n_platform))
    min(min(n_platform))
    figure;imagesc(d_amb,res,n_platform,[1 50]);colorbar;colormap jet
    title(['Minimum Number of Platforms Required (' modd ')'])
    xlabel('Nearest Ambiguity Location (along n) (m)')
    ylabel('4 dB Resolution (along n) (m)')
    set(gca,'fontsize',12);set(gca,'YDir','normal')
    xlim([min(d_amb) max(d_amb)])
    ylim([min(res) max(res)])
end
dock