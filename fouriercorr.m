% This function references from the fourier_shell_correlation within PyNX,
% which is created by Julio Cesar da Silva (jdasilva@esrf.fr) and modified
% by Vincent Favre-Nicolin (favre@esrf.fr) at ESRF.

function fouriercorr(ptycho_package1,ptycho_package2,row_ROI,col_ROI)
    % ROI in pixel [row_range,col_range]
    pixel_size = ptycho_package1.object_info.x_res;
    img1 = ptycho_package1.object_info.real_space(row_ROI,col_ROI);
    img2 = ptycho_package2.object_info.real_space(row_ROI,col_ROI);
    % ---Compute FSC and threshold---
    % default value
    snrt = [0.2071,0.5]; % 0.2071 for 1/2 bit and 0.5 for 1 bit.
    ring_thick = 3;
    rad_apod = 60;
    axial_apod = 20;
    % Apodization
    n = size(img1);
    circular_region = circle(img1,rad_apod);
    if ndims(img1) == 2
        window = circular_region;
    elseif ndims(img1) == 3
        fprintf('3D not yet.')
    end
    
    % FSC comoutation
    F1 = fft_custom(img1.*window);
    F2 = fft_custom(img2.*window);
    index = ringthickness(img1);
    [f,fnyquist] = nyquist(img1);
    
    C  = zeros(1,length(f));
    C1 = zeros(1,length(f));
    C2 = zeros(1,length(f));
    npts = zeros(1,length(f));
    
    for ii = 1:length(f)
        freq = f(ii);
        if ring_thick == 0
            mask = (index == freq);
        else
            mask = (index >= (freq - ring_thick /2)) & (index <= (freq + ring_thick / 2));
        end
        
        auxF1 = F1(mask);
        auxF2 = F2(mask);
        
        C(ii) = sum(auxF1 .* conj(auxF2));
        C1(ii) = sum(auxF1 .* conj(auxF1));
        C2(ii) = sum(auxF2 .* conj(auxF2));
        npts(ii) = sum(mask(:));
    end
    FSC = abs(C)./sqrt(C1.*C2);
    
    T = zeros(length(snrt),length(npts)); % all consider as multi threshold case
    for i = 1:length(snrt)
        snrt_i = snrt(i);
        Tnum = (snrt_i + (2 * sqrt(snrt_i) ./ sqrt(npts)) + 1 ./ sqrt(npts));
        Tden = (snrt_i + (2 * sqrt(snrt_i) ./ sqrt(npts)) + 1);
        T(i,:) = Tnum ./ Tden;
    end

    % determine resolution
    difference = FSC-T;
    s = log10(pixel_size);
    if     s < -6
        unit_name = 'nm';
        s = 1e9;
    elseif s < -3
        unit_name = 'um';
        s = 1e6;
    elseif s < 0
        unit_name = 'mm';
        s = 1e3;
    else
        unit_name = 'm';
        s = 1;
    end
   for ii = 1:length(snrt)
       cut_index = find(difference(ii,:)<0,1);
       resolution = pixel_size * s / (f(cut_index)/fnyquist);
       fprintf('SNRT = %f, resolution = %.2f %s\n',snrt(ii),resolution,unit_name);
   end
    
    FSCPlot(f,fnyquist,FSC,T,pixel_size)
    
    
function FSCPlot(f,fnyquist,FSC,T,pixel_size)
    s = log10(pixel_size);
    if     s < -6
        unit_name = 'nm';
        s = 1e9;
    elseif s < -3
        unit_name = 'um';
        s = 1e6;
    elseif s < 0
        unit_name = 'mm';
        s = 1e3;
    else
        unit_name = 'm';
        s = 1;
    end

    x1 = f/fnyquist; % spatial frequency/Nyquist
    x2 = pixel_size * s ./x1; % resolution
    fig = figure;
    fig.Position = [10 10 800 700];
    plot(f/fnyquist,FSC,'linewidth',4)
    for ii = 1:min(size(T))
        hold on
        plot(f/fnyquist,T(ii,:),'--','linewidth',3)
        hold off
    end
    legend('FSC','1/2 bit','1 bit','Box','off')
    
    % spatial frequency/Nyquist part
    ax1 = gca;
    ax1.FontSize = 14;
    box off;
    xlabel(ax1,'Spatial frequency/Nyquist','FontSize',14);
    ylabel(ax1,'Magnitude','FontSize',14);
    set(ax1,'XAxisLocation','bottom','LineWidth',2);
    set(ax1,'YTick',0:0.2:1);
    %ax1.XAxis.TickLength = [0 0];
    % resolution part
    %first_x1_ticklabel = str2double(ax1.XAxis.TickLabels{2});
    %first_resolution_ticklabel = pixel_size*s/first_x1_ticklabel;
    %resolution_ticks = [round(first_resolution_ticklabel):-5:10.1,10:-1:1]
    resolution_ticks = fliplr([1:1:10,20,40]);
    spatial_freq_ticks =  pixel_size*s./resolution_ticks;
    ax2 = axes('Position',get(ax1,'Position'),'Color','none');
    set(ax2, 'XAxisLocation', 'top','YAxisLocation','Right','Linewidth',2);
    set(ax2,'XLim',get(ax1,'XLim'),'YLim',get(ax1,'YLim'));
    set(ax2,'XTick',spatial_freq_ticks,'YTick',get(ax1,'YTick'));
    xticklabels(ax2,arrayfun(@(x) sprintf('%d',x),resolution_ticks,'UniformOutput',false));
    ax2.FontSize = 14;
    xlabel(ax2,sprintf('Resolution (%s)',unit_name),'FontSize',14);
    grid on

    
    
    % basic plot
    %{
    ax1 = gca;
    xlabel(ax1,'spatial frequency/Nyquist','FontSize',24);
    ylabel(ax1,'Magnitude','FontSize',24);
    ax2 = axes('Position',get(ax1,'Position'),'Color','none');
    set(ax2, 'XAxisLocation', 'top','YAxisLocation','Right');
    set(ax2,'XLim',get(ax1,'XLim'),'YLim',get(ax1,'YLim'));
    set(ax2,'XTick',get(ax1,'XTick'),'YTick',get(ax1,'YTick'));
    ax2xTickLabels=cellfun(@(X)num2str((pixel_size*s/str2double(X))),get(ax1,'xticklabels'),'UniformOutput',false);
    ax2xTickLabels{1} = '';
    get(ax1,'xticklabels');
    set(ax2,'XTickLabel',ax2xTickLabels,'YTickLabel',{});
    xlabel(ax2,unit_name,'FontSize',24);
    %}
    
function index = ringthickness(img)
    n = size(img);
    nmax = max(n);
    x = (fix( -n(2)/2 ) : ceil( n(2)/2 -1 )) * floor(nmax/2) / floor(n(2)/2);
    y = (fix( -n(1)/2 ) : ceil( n(1)/2 -1 )) * floor(nmax/2) / floor(n(1)/2);
    if ndims(img) == 3
        z = (fix( -n(3)/2 ) : ceil( n(3)/2 -1 )) * floor(nmax/2) / floor(n(3)/2);
        [X,Y,Z] = meshgrid(x,y,z);
        sumsquares = X.^2 + Y.^2 + Z.^2;
    elseif ndims(img) == 2
        [X,Y] = meshgrid(x,y);
        sumsquares = X.^2 + Y.^2;
    end
    index = round(sqrt(sumsquares));

function fft_img = fft_custom(img)
    fft_img = ifftshift(fftn(fftshift(img)));
    


function [f,fnyquist] = nyquist(img)
    nmax = max(size(img));
    fnyquist = floor(nmax / 2);
    f = 0:fnyquist;


function window = apodization(img)
    n = size(img);
    ndim = ndims(img);
    if ndim == 2
        window = hanning(n(1))*hanning(n(2)); % equal to np.outer(hanning(n[0]),hanning(n[1])) in the original code.
    elseif ndim == 3
        error('3D not yet.')
    end
    
function circle_region = circle(img,rad_apod)
    n = size(img);
    ndim = ndims(img);
    if ndim == 2
        shape_x = n(2);
        shape_y = n(1);
    elseif ndim == 3
        error('3D not yet.');
    end
    x_array = 1:shape_x;
    y_array = 1:shape_y;
    x_cen = round(shape_x/2);
    y_cen = round(shape_y/2);
    [X,Y] = meshgrid(x_array-x_cen,y_array-y_cen);
    circle_region = 1 - radtap(X,Y,rad_apod,x_cen-rad_apod);
    
function taperfunc = radtap(X,Y,tappix,zerorad)
    % Creates a central cosine tapering.
    % It receives the X and Y coordinates, tappix is the extent of
    % tappering, zerorad is the radius with no data (zeros).
    
    tau = 2*tappix; % period of cosine funcition (only half a period is used)
    
    R = sqrt(X.^2 + Y .^2);
    taperfunc = 0.5 * ( 1 + cos(2*pi*(R - zerorad - tau/2)/tau) );
    taperfunc = (R > zerorad + tau / 2) * 1.0 + taperfunc .* (R <= zerorad + tau / 2);
    taperfunc = taperfunc .* (R >= zerorad);
 
        
function w = hanning(N)
    % hann_custom generates a Hann window
    % N: the length of the window
    % Returns:
    % w: Hann window of length N
    
    if N <= 0
        error('Window length N must be positive.');
    end
    n = (0:N-1)';  % Sample indices
    w = 0.5 * (1 - cos(2 * pi * n / (N - 1)));
