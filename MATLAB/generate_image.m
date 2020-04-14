%% generates image from raw data using time-domain back-projection
close all; clear; clc
%% scene extend and resolution (relative to reference point on ground, where perpendicular to baseline line intersects zero height ground)
xsc=20:0.5:100;dxsc=xsc(2)-xsc(1);xsc_min=min(xsc);
zsc=5:0.5:45;dzsc=zsc(2)-zsc(1);zsc_min=min(zsc);
K=length(xsc); % number of scence pixels in horizontal direction
M=length(zsc); % number of scence pixels in vertical direction
I=zeros(M,K); % image 2D array
%% load raw data
load('raw.mat'); %'mode','Pr','H','wl','xs','zs','xt','xt_rel','zt','at','tx_element','theta','alpha','t','tau','B','Ravg'
N=length(xs); % number of platforms
T=length(xt); % number of targets
Nt=length(t); % number of time samples
dt=t(2)-t(1); % time resolution
%% focus using time-domain back-projection
c=3e8; % speed of light (m/s)
w=waitbar(0);
if mode==1 % SAR
    for i=1:K
        x=xsc(i)+H*tand(theta);
        for j=1:M
            waitbar(((i-1)*M+(j-1))/(K*M));
            z=zsc(j);
            Ri=((x-xs).^2+(z-zs).^2).^0.5;
            delays=2*(Ri-Ravg)/c;delays_ind=round(Nt/2+delays/dt); % delays for each platform for a specific scene pixel (s)              
            delays_ind2=(1:N)+(delays_ind-1)*N;
            I(j,i)=abs(sum(Pr(delays_ind2).*exp(1i*4*pi/wl*Ri)));
        end
    end
elseif mode==2% SIMO
    for i=1:K
        x=xsc(i)+H*tand(theta);
        for j=1:M
            waitbar(((i-1)*M+(j-1))/(K*M));
            z=zsc(j);
            R1=((x-xs(tx_element)).^2+(z-zs(tx_element)).^2).^0.5;
            R2=((x-xs).^2+(z-zs).^2).^0.5;
            delays=(R1+R2-2*Ravg)/c;delays_ind=round(Nt/2+delays/dt); % delays for each platform for a specific scene pixel (s)              
            delays_ind2=(1:N)+(delays_ind-1)*N;
            I(j,i)=abs(sum(Pr(delays_ind2).*exp(1i*2*pi/wl*(R1+R2))));
        end
    end  
elseif mode==3 % MIMO
    xs_rx=repmat(xs,1,N);
    zs_rx=repmat(zs,1,N);
    xs_tx=kron(xs,ones(1,N));
    zs_tx=kron(zs,ones(1,N));
    for i=1:K
        x=xsc(i)+H*tand(theta);
        for j=1:M
            waitbar(((i-1)*M+(j-1))/(K*M));
            z=zsc(j);
            R1=((x-xs_tx).^2+(z-zs_tx).^2).^0.5;
            R2=((x-xs_rx).^2+(z-zs_rx).^2).^0.5;
            delays=(R1+R2-2*Ravg)/c;delays_ind=round(Nt/2+delays/dt); % delays for each platform for a specific scene pixel (s)
            delays_ind2=(1:N*N)+(delays_ind-1)*N*N;
            I(j,i)=abs(sum(Pr(delays_ind2).*exp(1i*2*pi/wl*(R1+R2))));
        end
    end
end
close(w)
%% PLOTS
% plot input scence
D=zeros(M,K);
xt_ind=round((xt_rel-xsc_min)/dxsc);zt_ind=round((zt-zsc_min)/dzsc);
for i=1:T;D(zt_ind(i)+1,xt_ind(i)+1)=at(i);end
figure;hold on;imagesc(xsc,zsc,D);xlabel('horizontal distance (m)');ylabel('vertical distance (m)');title('Input Scene (linear amplitude)')
colormap gray;xlim([min(xsc) max(xsc)]);ylim([min(zsc) max(zsc)]);colorbar
% plot image
tx_el='';if mode==1;ttl='(SAR)';end;if mode==2;ttl='(SIMO)';tx_el=[' TX element:' num2str(tx_element)];end;if mode==3;ttl='(MIMO)';end
figure;hold on;xlabel('horizontal distance (m)');ylabel('vertical distance (m)');title(['Generated Scene (linear amplitude) ' ttl tx_el])
imagesc(xsc,zsc,abs(I));colormap jet;xlim([min(xsc) max(xsc)]);ylim([min(zsc) max(zsc)]);colorbar
% plot(xt_rel,zt,'ko','MarkerSize',5,'linewidth',2) % to add marker for each target
saveas(gcf,'image.jpg') % save image as jpg