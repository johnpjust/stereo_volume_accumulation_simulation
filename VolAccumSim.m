function [volsumGT, volsumMeas, slatvecindcnt]=VolAccumSim(N, shape, xfill, hfill, elevSpd, slatskip, slat_pat)
%% set Geometry and Sampling-Rate parameters %%
%quantized/smallest x-increment 
deltax = 0.1; %inches ***DO NOT CHANGE***

%slat spacing layout options for one rotation
link_size = 2; %inches per link

%slat pattern
% slat_pat = 1; %0 = typical 23 slats;; 1=options 1: 24 slats;; 2=option 2:  24 slats;; 3=options 3: 26 slats
typ_slat_spa = 20; %20" used for fill equation

%sampling frequency
fs = 7.5;

%slat/elev dimensions (height, width, depth)
Hs=7;
Welv = 34.3; %Welv=34.25;Ws=33.25;
Ds=1.5;
Helev = 16.3; %16.25

%FOV = 51 deg, Angle relative to elev = 37.3, Hc=.987 cm,
%offsetFromTop of FOV=3cm??
%camera mount/view (90-deg height, forward/backward angles, wind len)
Hc = 38.7; %820.3644mm default 98cm
theta_c_f = 37.3-51/2; %from perp -- must be > 0
theta_c_b = 37.3+51/2; %from perp -- must be > theta_c_f and << 90
x0 = 0; %reference point on x-axis

%% Initialize Run/Fill Parameters -- User Entered %%

%slat layout AUTO-CALC -- do not change this part
switch slat_pat
    case 1
        slat_pattern = link_size*[10*ones(9,1);12;10*ones(14,1)]; %options 1: 24 slats
    case 2
        slat_pattern = link_size*[10*ones(9,1);12;10*ones(13,1);12]; %option 2:  24 slats
    case 3
        slat_pattern = link_size*[10*ones(9,1);12;10*ones(15,1);12]; %options 3: 26 slats
    otherwise
        slat_pattern = link_size*10*ones(23,1); %typical 23 slats
end

% % % %change these parameters
% % % %number of rotations
% % % N = 4;
% % % %specify volume stacking shape
% % % shape = 2; %1=rectangular, 2=triangular, 3=sqrt, 4=quadratic, 5=random
% % %             %shape = unidrnd(4,[length(slat_pattern)*N,1])
% % % %specify % in decimal x-fill for each slat (% of x-space between slats)
% % % xfill = .35; %if >= 1 then random
% % % %specify % in decimal H-fill for each slat (% of sidewall height)
% % % hfill = .35; %if >= 1 then random
% % % elevSpd = 1.9; %m/s 5 == random
% % % slatskip = 2; %1=sequential, 2=alternating, 3=every third is zero

%% All calcs below are automatic from set parameters %%

%calc parameters according to user-entered values
xc = round(Hc*tand(theta_c_f),1);%158.5713 #camera loc on x-axis
xL = min(xc+23.6,Hc*tand(theta_c_b)); %end of FOV on x-axis

%camera window vectors (lower heights and upper heights at each position)
camWindLowlim = [Hc*(1-(x0:deltax:xc)/xc),zeros(1,uint32((xL-xc)/deltax))]';
camWindHighlim = [Hc*(1-(x0:deltax:xL)/(Hc*tand(theta_c_b)))]';
%% "Ground Truth" vol vector creation %%
%volume vector 
    %slat vector (number of total slats to run, positions, full/empty)
    slatlen = uint32(Ds/deltax);
    slatind = uint32(cumsum(repmat(slat_pattern,[N,1]))/deltax);
    
    slatvec = zeros(slatind(end),1);

    %xfill, hfill vectors
    if xfill < 1
        xfill = typ_slat_spa*xfill.*ones(length(slat_pattern)*N,1);
    else
        xfill = typ_slat_spa*unifrnd(0,.5,[length(slat_pattern)*N,1]);
    end
    if hfill < 1
        hfill = Helev*hfill.*ones(length(slat_pattern)*N,1);
    else
        hfill = Helev*unifrnd(0,.5,[length(slat_pattern)*N,1]);
    end

    if slatskip == 2
        hfill(2:2:end) = 0;
    elseif slatskip == 3
        hfill(setdiff(1:length(hfill),1:3:length(hfill)))=0;
    end
    
    % fill with slats for true volume if desirable
%     emptyslatvec = zeros(slatind(end),1);
%     for n=1:length(slatind)
%         emptyslatvec(slatind(n):(slatind(n)+slatlen)) = Hs;
%     end

    %fill slat vector with volume based on shape, xfill, hfill
    if shape > 4
        shape = unidrnd(4,[length(slat_pattern)*N,1]);
    else
        shape = shape*ones(length(slat_pattern)*N,1);
    end
    for n=1:length(slatind)
        volind = double(slatind(n)-round(xfill(n)/deltax):slatind(n));
        xtemp = deltax*(volind-min(volind));
        switch shape(n)
            case 2
                slatvec(volind) = xtemp*hfill(n)/xfill(n); %triangular
            case 3
                slatvec(volind) = hfill(n)*sqrt((xtemp/xfill(n))); %quadratic_out
            case 4
                slatvec(volind) = hfill(n)*(xtemp/xfill(n)).^2; %quadratic_in
            otherwise
                slatvec(volind) = hfill(n); %quare
        end
    end

    %randomize starting position with pre-appended zeros
    randstartpos = zeros(uint32(xL/deltax) + unidrnd(500),1)/deltax;
    slatvec = [randstartpos;slatvec];
    slatind = slatind + length(randstartpos);
    
    %plot vol vector to check
%     plot((0:(length(slatvec)-1))*deltax,slatvec);
% elev speed vector (Gaussian, step changes, etc) - make at least as 
% large as slat vector, maybe more

% plot volume vector to visual it looks appropriate (highlight slats)


%% accumulate volume (shift, multiply, add)
elevSpd = elevSpd*39.3701.*ones(round(length(slatvec)),1);
slatvecindcnt = zeros(N*length(slatvec),1);
repeatcnt = 1;
% accum loop
n=1;
volsumMeas = 0;
while n < length(slatvec)
    %calculate bool window vector
    if length(slatvec)-n >= length(camWindHighlim)
        vol = slatvec(n:(n+length(camWindHighlim)-1));
    else
        vol = slatvec(n:end);
    end
    
    slatindsub = double(slatind(slatind >= n & slatind <= (n+length(vol)-round(Ds/deltax))));
    seq = [n:(n+length(camWindHighlim)-1)]';
     Hshad = zeros(2*length(vol),1);
     if ~isempty(slatindsub)
        slatindsub = slatindsub - n+1;
        for n_in=1:length(slatindsub)
            tan_theta_slat = (Hc-Hs)/(deltax*slatindsub(n_in)+Ds);
            tan_theta_vol = (Hc-vol(slatindsub(n_in)))/(deltax*slatindsub(n_in));
            xshad_slat = Hs/tan_theta_slat;
            xshad_vol = vol(slatindsub(n_in))/tan_theta_vol;
            Hshad_slat = [Hs*ones(1,round(Ds/deltax)-1),Hs-tan_theta_slat*(0:deltax:xshad_slat)]';
            Hshad_vol = [vol(slatindsub(n_in))-tan_theta_vol*(0:deltax:xshad_vol)]';
            shadlendiff = length(Hshad_slat) - length(Hshad_vol);
            if shadlendiff > 0
                Hshad_vol = [Hshad_vol;zeros(shadlendiff,1)];
            elseif shadlendiff < 0
                Hshad_slat = [Hshad_slat;zeros(abs(shadlendiff),1)]; 
            end
            Hshad(slatindsub(n_in):(slatindsub(n_in)+length(Hshad_vol)-1)) = max(Hshad_vol,Hshad_slat);
        end
        Hshad = Hshad(1:(length(vol)));
        camWindLowlim_adj = max(camWindLowlim(1:length(vol)),Hshad);
        camWindHighlim_adj = camWindHighlim(1:length(vol));

        vol(vol<camWindLowlim_adj | vol>camWindHighlim_adj) = 0;
        if sum(vol) > 0
            slatvecindcnt(repeatcnt:repeatcnt+length(seq((vol>camWindLowlim_adj & vol<camWindHighlim_adj)))-1) = seq((vol>camWindLowlim_adj & vol<camWindHighlim_adj));
            repeatcnt = repeatcnt + length(seq((vol>camWindLowlim_adj & vol<camWindHighlim_adj)));
        end
%         figure(2)
%         plot(slatvecindcnt(slatvecindcnt ~=0))
%         
%         figure(1)
%         plot(vol,'b.','Linewidth',3);hold on;plot(camWindHighlim_adj,'r');
%         plot(camWindLowlim_adj,'r');plot(Hshad,'k','Linewidth',3);
%         plot(vol,'g')
%         clf
    end
    %multiply bool window and vol vector over correct indecies (ignore slat
    %volumes)
    
    %update positional index of camera window w.r.t vol vector using elevspd
    volsumMeas = volsumMeas + sum(vol);
    n = n+round(elevSpd(n)/(fs*deltax));
end
volsumGT=sum(slatvec);
slatvecindcnt(slatvecindcnt==0) = [];

