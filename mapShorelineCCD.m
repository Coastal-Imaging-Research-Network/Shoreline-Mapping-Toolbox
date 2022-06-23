function sl = mapShorelineCCD(xgrid,ygrid,Iplan,transects,plotoption,editoption)
%
%function sl = mapShorelineCCD(xgrid,ygrid,Iplan,transects,editoption)
%
%function that maps a shoreline from a rectified planview image using the
%Colour Channel Divergence (CCD) technique (Red minus the Blue channel - can be adapted).
%
%xgrid -- a 1xM array of local UTM eastings (note: this is different from
%               Argus cooedinates)
%ygrid -- a 1xN array of local UTM northings
%Iplan -- a MxNx3 uint8 image of the planview
%transects -- a structure that seeds the shoreline search algorithm
%using a series of cross-shore transects spaced at ~5m apart
% transects.x -- the start and end points of each transect
% x-coordinates (2 x M matrix, where 1st row = start, 2nd row = end)
% transects.y -- the start and end points of the transects
% y-coordinates (2 x M matrix, where 1st row = start, 2nd row = end)
%editoption == 1 if you want to manually edit the shoreline (using
%               interactive point editor)
%
%Created by Mitch Harley
%June 2018


if nargin==4
    plotoption=1;
    editoption=1;
elseif nargin==5
    editoption=0;
end

if plotoption==1
    f1=figure;
    image(xgrid,ygrid,Iplan)
    axis image; axis xy
end
[X,Y] = meshgrid(xgrid,ygrid);

%Find threshold
P = improfile(xgrid,ygrid,Iplan,transects.x,transects.y); %Sample pixels at transects to determine threshold
Iblack = find((P(:,:,1)-P(:,:,3))==0); %Remove black pixels from interfering with method (at edge of image)
P(Iblack,:,:) = [];
if isempty(P)
    disp('Image appears to be a greyscale image so shoreline could not be mapped using CCD approach. Returning empty values')
    sl.x =[];
    sl.y = [];
    sl.method = 'CCD';
    sl.threshold = [];
else
    [pdf_values,pdf_locs] = ksdensity(P(:,:,1)-P(:,:,3)); %find smooth pdf of CCD (Red minus Blue) channel
    xlabel_type = 'Red minus blue';
    thresh_weightings = [1/3 2/3]; %This weights the authomatic thresholding towards the sand peak (works well in SE Australia)
    %[peak_values,peak_locations]=findpeaks(pdf_values,pdf_locs); %Find peaks
    thresh_otsu = multithresh(P(:,:,1)-P(:,:,3)); %Threshold using Otsu's method
    I1 = find(pdf_locs<thresh_otsu);
    [~,J1] = max(pdf_values(I1));
    I2 = find(pdf_locs>thresh_otsu);
    [~,J2] = max(pdf_values(I2));
    thresh = thresh_weightings(1)*pdf_locs(I1(J1)) + thresh_weightings(2)*pdf_locs(I2(J2)); %Skew average towards the positive (i.e. sand pixels)
    disp(['Sand/water threshold determined to be ' num2str(thresh, '%0.0f')])
    
    if plotoption==1
        f2 = figure;
        plot(pdf_locs,pdf_values)
        hold on
        plot(pdf_locs([I1(J1) I2(J2)]),pdf_values([I1(J1) I2(J2)]),'ro')
        YL = ylim;
        plot([thresh thresh], YL,'r:','linewidth',2)
        %plot([thresh_otsu thresh_otsu], YL,'g:','linewidth',2)
        xlabel(xlabel_type,'fontsize',10)
        ylabel('Counts','fontsize',10)
    end
    
    %Extract contour
    RminusBdouble = double(Iplan(:,:,1))- double(Iplan(:,:,3));
    
    %Mask out data outside of ROI (ROI determined from transects)
    ROIx = [transects.x(1,:) fliplr(transects.x(2,:))]; %Transects file determines the ROI
    ROIy = [transects.y(1,:) fliplr(transects.y(2,:))]; %Transects file determines the ROI
    Imask = ~inpoly([X(:) Y(:)],[ROIx',ROIy']); %use the function inpoly instead of inpolygon as it is much faster
    RminusBdouble(Imask) = NaN; %Mask data
    %Imask = find(RminusBdouble==0); %Also remove regions of black colour
    %RminusBdouble(Imask) = NaN; %Mask data
    c = contours(X,Y,RminusBdouble,[thresh thresh]);
    
    
    %Now look at contours to only find the longest contour (assumed to be the
    %shoreline)
    II = find(c(1,:)==thresh);
    if II==1 %If only one line
        startI = 2;
        endI = size(c,2);
    else
        D = diff(II);
        [~,J] = max(D); %Select contour that is the longest continuous contour
        if J == 1
            startI = 2;
        else
            startI = 1+J+sum(D(1:J-1));
        end
        endI = startI+D(J)-J-1;
    end
    xyz.y = c(1,startI:endI)';
    xyz.x = c(2,startI:endI)';
    points = [xyz.y xyz.x];
    
    %Now loop through transects to extract shorelines only at the transects
    sl.x = NaN(1,length(transects.x));
    sl.y = NaN(1,length(transects.y));
    warning off
    for i = 1:length(transects.x)
        angle = atan(diff(transects.y(:,i))/diff(transects.x(:,i)));
        points_rot = rotatePoints(points,angle,[transects.x(1,i) transects.y(1,i)],'rads');
        max_distance = sqrt(diff(transects.y(:,i))^2+ diff(transects.x(:,i))^2);
        I = find(points_rot(:,2)>-1&points_rot(:,2)<1&points_rot(:,1)>0&points_rot(:,1)<max_distance); %Only find points greater than zero so that they are in the roi
        if ~isempty(I)
            [~,Imin] = min(points_rot(I,1));
            %[~,Imin] = max(points_rot(I,1));
            
            %Do test to see if point falls in a black pixel - remove if yes
            blacktest_tol = 2; %Tolerance in metres
            points_blacktest = [points_rot(I(Imin),1)+blacktest_tol,points_rot(I(Imin),2)]; %Select a point that is a certain tolerance seaward of the transect to see if it is black (default 2m)
            points_blacktest_unrot = unrotatePoints(points_blacktest,angle,[transects.x(1,i) transects.y(1,i)],'rads');
            PP = interp2(X,Y,RminusBdouble,points_blacktest_unrot(1),points_blacktest_unrot(2)); 
            if PP~=0 %If the "black test" pixel is not black
                sl.x(i) = points(I(Imin),1);
                sl.y(i) = points(I(Imin),2);
            end         
        end
        
    end
    
    
    
    if editoption ==1 && plotoption == 1
        figure(f1)
        h = impoly(gca,[sl.x' sl.y'],'closed',0);
        disp('Please edit your shoreline now or press any key to continue')
        pause
        newpos = getPosition(h);
        sl.x = newpos(:,1);
        sl.y = newpos(:,2);
    elseif editoption ==0
        sl.x = sl.x';
        sl.y = sl.y';
    end
    
    %Remove NaNs
    I = find(isnan(sl.x));
    sl.x(I) = [];
    sl.y(I) = [];
    sl.method = 'CCD';
    sl.threshold = thresh;
end
