%% generates raw data
close all; clear; clc
%% parameters
mode=1; %1: SAR, 2:SIMO, 3:MIMO
tx_element=1; % which element transmits for SIMO (max value N)
c=3e8; % speed of light (m/s)
f=1e9; % center frequency (Hz)
wl=c/f; % wavelength (m)
H=500e3; % altitude (m)
theta=60; % look angle wrt nadir (deg)
alpha=30; % baseline tilt wrt horizontal (deg)
B=80e6; % bandwidth (Hz) (used for matched filter output) 
tau=10e-6; % pulse width (s) (used for matched filter output)
dt=1e-9; % receive window resolution (s)
%% platforms and targets
% platform locations
s=-10e3:1000:10e3; % uniform platform spacing
% s=1e3*[-40 -37 -30 -17 4 21 30 35 38 39 40]; % non-uniform platform spacing
% target locations relative to reference point on ground (where perpendicular to baseline line intersects zero height ground)(assuming all targets are in same azimuth pixel) (xt and zt should be zero or positive)
xt_rel=[30 30 30 30 30 32 36 44 60 92]; % target horizontal location (relative to reference point on ground)
zt=    [10 12 16 24 40 40 40 40 40 40]; % target vertical location (relative to ground)
at=    [10 10 15 20 20 10 15 15 20 20]; % target reflectivities (relative)
N=length(s); % number of platforms
T=length(xt_rel); % number of targets
zs=H+s*sind(alpha);xs=s*cosd(alpha); % platform coordinates
xt=xt_rel+H*tand(theta); % target horizontal coordinate relative to nadir point of baseline center
% plot platform locations
xmin=min(min(xs/1e3),min(zs/1e3)-H/1e3);xmax=max(max(xs/1e3),max(zs/1e3)-H/1e3);
zmin=min(min(xs/1e3),min(zs/1e3)-H/1e3);zmax=max(max(xs/1e3),max(zs/1e3)-H/1e3);
figure;scatter(xs/1e3,(zs-H)/1e3,100,'filled');
if theta<45;x2=abs(zmin)*tand(theta);z2=zmin;elseif theta>=45;x2=xmax;z2=-xmax/tand(theta);end;line([0 x2],[0 z2],'Color','red','linewidth',2)
xlim([xmin xmax]);ylim([zmin zmax])
set(gca,'fontsize',12);title('Platform Positions and Look Angle');xlabel('horizontal distance - x (km)');ylabel('vertical distance - z (km)')
%% matched filter output and receive window
t_mf=-tau:dt:tau; % mf is zero outside +/- tau
mf=sinc((1-abs(t_mf)/tau).*t_mf*B).*(1-abs(t_mf)/tau); % equation for LFM chirp matched filter output
figure;plot(t_mf,20*log10(abs(mf)),'linewidth',2);ylim([-100 0]);
title('Matched Filter');ylabel('power (dB)');xlabel('time (s)');set(gca,'fontsize',12)
% receive window with matched filter centered on median range (Ravg)
Rmin=((max(xs)-min(xt))^2+(min(zs)-max(zt))^2)^0.5;
Rmax=((min(xs)-max(xt))^2+(max(zs)-min(zt))^2)^0.5;
Ravg=(Rmin+Rmax)/2;
Trx_min=2*(Rmin-Ravg)/c-tau;% receive window duration (s) 
Trx_max=2*(Rmax-Ravg)/c+tau;% receive window duration (s) 
t=Trx_min:dt:Trx_max;
Nt=length(t);
trx=zeros(1,Nt);
t_mf_ind=round(Nt/2-tau/dt):round(Nt/2+tau/dt);
trx(t_mf_ind)=mf; 
%% theoretical resolutions
if mode==1;p=2;end
if mode==2;p=1;end
if mode==3;p=1.4;end
dr=c/2/B % 4dB slant range resolution (m)
dr_x=dr/sind(theta) % horizontal projection of slant range resolution
dr_z=dr/cosd(theta) % vertical projection of slant range resolution
Lsa=(max(s)-min(s))*cosd(abs(theta-alpha));dn=wl*Ravg/p/Lsa % 4dB resolution along perpendicular baseline
dn_x=dn/cosd(theta) % horizontal projection of resolution along perpendicular baseline
dn_z=dn/sind(theta) % vertical projection of resolution along perpendicular baseline
%% random noise
na=0; % additive noise amplitude (set to 0 for no additive random noise)
sp_std=0;sp_mn=1; % multiplicative noise std.dev and mean (set std=0 and mean=1 for no multiplicative random noise)
nt=(na/2^0.5)*(randn(N,Nt)+1i*randn(N,Nt)); % additive random noise
st=(sp_std/2^0.5)*(randn(1,T)+1i*randn(1,T))+sp_mn; % multiplicative random noise
%% raw data simulation
mf_delayed=zeros(T,Nt);
w=waitbar(0);
if mode==1 % SAR (ping-pong)
    Pr=zeros(N,Nt);
    for i=1:N
        waitbar(i/N);
        Ri=((xt-xs(i)).^2+(zt-zs(i)).^2).^0.5;
        delays=2*(Ri-Ravg)/c;delays_ind=round(delays/dt); % delays for each target for a specific platform (s)
        for j=1:T
            if delays_ind(j)>=0
                shifted_trx=[zeros(1,delays_ind(j)) trx(1:end-delays_ind(j))];
            elseif delays_ind(j)<0
                shifted_trx=[trx(1-delays_ind(j):end) zeros(1,-delays_ind(j))];
            end
            mf_delayed(j,:)=at(j).*st(j).*exp(-1i*4*pi/wl*Ri(j)).*shifted_trx;
        end
        Pr(i,:)=sum(mf_delayed+nt(i,:),1);
    end
elseif mode==2 % SIMO (semi-active)
    Pr=zeros(N,Nt);
    for i=1:N
        waitbar(i/N);
        R1=((xt-xs(tx_element)).^2+(zt-zs(tx_element)).^2).^0.5;
        R2=((xt-xs(i)).^2+(zt-zs(i)).^2).^0.5;
        delays=(R1+R2-2*Ravg)/c;delays_ind=round(delays/dt); % delays for each target for a specific RX platform (s)
        for j=1:T
            if delays_ind(j)>=0
                shifted_trx=[zeros(1,delays_ind(j)) trx(1:end-delays_ind(j))];
            elseif delays_ind(j)<0
                shifted_trx=[trx(1-delays_ind(j):end) zeros(1,-delays_ind(j))];
            end
            mf_delayed(j,:)=at(j).*st(j).*exp(-1i*2*pi/wl*(R1(j)+R2(j))).*shifted_trx;
        end
        Pr(i,:)=sum(mf_delayed+nt(i,:),1);
    end
elseif mode==3 % MIMO (active)
    Pr=zeros(N,N,Nt);
    for i=1:N
        for j=1:N
            waitbar(((i-1)*N+j)/(N*N));
            R1=((xt-xs(i)).^2+(zt-zs(i)).^2).^0.5;
            R2=((xt-xs(j)).^2+(zt-zs(j)).^2).^0.5;
            delays=(R1+R2-2*Ravg)/c;delays_ind=round(delays/dt); % delays for each target for a specific TX-RX platform pair (s)
            for k=1:T
                if delays_ind(k)>=0
                    shifted_trx=[zeros(1,delays_ind(k)) trx(1:end-delays_ind(k))];
                elseif delays_ind(k)<0
                    shifted_trx=[trx(1-delays_ind(k):end) zeros(1,-delays_ind(k))];
                end
                mf_delayed(k,:)=at(k).*st(k).*exp(-1i*2*pi/wl*(R1(k)+R2(k))).*shifted_trx;
            end
            Pr(i,j,:)=sum(mf_delayed+nt(i,:),1);
        end
    end
end
close(w)    
%% save data and plot figures
if mode==1 || mode==2
    if mode==1;ttl='SAR';end;if mode==2;ttl='SIMO';end
    figure;hold on;imagesc(t,1:N,20*log10(abs(Pr)));colormap jet;colorbar;title([ttl ' Raw Data Amplitude (dB)']);ylabel('platform');xlabel('receive window (s)');set(gca,'fontsize',12);ylim([1 N])
    figure;hold on;imagesc(t,1:N,angle(Pr)*180/pi);colormap parula;colorbar;title([ttl ' Raw Data Phase (deg)']);ylabel('platform');xlabel('receive window (s)');set(gca,'fontsize',12);ylim([1 N])
    save('raw.mat','mode','Pr','H','wl','xs','zs','xt','xt_rel','zt','at','tx_element','theta','alpha','t','tau','B','Ravg')
elseif mode==3
    Pr=reshape(Pr,N*N,Nt);
    figure;imagesc(t,1:N*N,20*log10(abs(Pr)));colormap jet;colorbar;title('MIMO Raw Data Amplitude (dB)');xlabel('receive window (s)');ylabel('transmit/receive platform pair');set(gca,'fontsize',12)
    figure;imagesc(t,1:N*N,angle(Pr)*180/pi);colormap parula;colorbar;hcb=colorbar;title(hcb,'deg');title('MIMO Raw Data Phase (deg)');xlabel('receive window (s)');ylabel('transmit/receive platform pair');set(gca,'fontsize',12)
    save('raw.mat','mode','Pr','H','wl','xs','zs','xt','xt_rel','zt','at','tx_element','theta','alpha','t','tau','B','Ravg', '-v7.3') %-v7.3 allows saving large file
end