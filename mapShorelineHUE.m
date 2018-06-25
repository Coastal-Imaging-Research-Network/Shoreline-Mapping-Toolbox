function sl = mapShorelineHUE(xgrid,ygrid,Iplan,transects,editoption)
%
%function sl = mapShorelineHUE(xgrid,ygrid,Iplan,transects,editoption)
%
%function that maps a shoreline from a rectified planview image using the
%Hue channel.
%
%xgrid -- a 1xM array of local UTM eastings
%ygrid -- a 1xN array of local UTM northings
%Iplan -- a MxNx3 uint8 image of the planview
%transects -- a structure that seeds the shoreline search algorithm
%using a series of cross-shore transects spaced at ~5m apart
   % transects.x -- the start and end points of each transect
   % x-coordinates (2 x M matrix, where 1st row = start, 2nd row = end)
   % transects.y -- the start and end points of the transects
   % y-coordinates (2 x M matrix, where 1st row = start, 2nd row = end)
%editoption == 1 if you want to manually edit the shoreline
%
%Created by Mitch Harley
%June 2018

   
if nargin==4
    editoption=1;
end

f1=figure;
image(xgrid,ygrid,Iplan)
axis image; axis xy
[X,Y] = meshgrid(xgrid,ygrid);

%Find threshold
HSV = rgb2hsv(Iplan);
P = improfile(xgrid,ygrid,HSV,transects.x,transects.y); %Sample pixels at transects to determine threshold
[pdf_values,pdf_locs] = ksdensity(P(:,:,1)); %find smooth pdf of Hue channel
xlabel_type = 'Hue Channel';
thresh_weightings = [1/3 2/3]; %This weights the authomatic thresholding towards the sand peak (works well in SE Australia)
[peak_values,peak_locations]=findpeaks(pdf_values,pdf_locs); %Find peaks
thresh_otsu = multithresh(P(:,:,1)); %Threshold using Otsu's method
f2 = figure;
plot(pdf_locs,pdf_values)
hold on
I1 = find(peak_locations<thresh_otsu);
[~,J1] = max(peak_values(I1));
I2 = find(peak_locations>thresh_otsu);
[~,J2] = max(peak_values(I2));
plot(peak_locations([I1(J1) I2(J2)]),peak_values([I1(J1) I2(J2)]),'ro')
%thresh = mean(peak_locations([I1(J1) I2(J2)])); %only find the last two peaks
thresh = thresh_weightings(1)*peak_locations(I1(J1)) + thresh_weightings(2)*peak_locations(I2(J2)); %Skew average towards the positive (i.e. sand pixels)

YL = ylim;
plot([thresh thresh], YL,'r:','linewidth',2)
%plot([thresh_otsu thresh_otsu], YL,'g:','linewidth',2)
xlabel(xlabel_type,'fontsize',10)
ylabel('Counts','fontsize',10)
disp(['Sand/water threshold determined to be ' num2str(thresh, '%0.0f')])

%Extract contour
HUEdouble = double(HSV(:,:,1));

%Mask out data outside of ROI (ROI determined from transects)
ROIx = [transects.x(1,:) fliplr(transects.x(2,:))]; %Transects file determines the ROI
ROIy = [transects.y(1,:) fliplr(transects.y(2,:))]; %Transects file determines the ROI
Imask = ~inpoly([X(:) Y(:)],[ROIx',ROIy']); %use the function inpoly instead of inpolygon as it is much faster
HUEdouble(Imask) = NaN; %Mask data
%Imask = find(HUEdouble==0); %Also remove regions of black colour
%HUEdouble(Imask) = NaN; %Mask data
c = contours(X,Y,HUEdouble,[thresh thresh]);
figure(f1)

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
            sl.x(i) = points(I(Imin),1);
            sl.y(i) = points(I(Imin),2);
        end
end

if editoption ==1
    h = impoly(gca,[sl.x' sl.y'],'closed',0);
    disp('Please edit your shoreline now or press any key to continue')
    pause
    newpos = getPosition(h);
    sl.x = newpos(:,1);
    sl.y = newpos(:,2);
end

%Remove NaNs
I = find(isnan(sl.x));
sl.x(I) = [];
sl.y(I) = [];
sl.method = 'HUE';
sl.threshold = thresh;