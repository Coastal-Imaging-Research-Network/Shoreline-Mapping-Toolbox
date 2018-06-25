function points_rotated = rotatePoints(points,angle,rotation_centre,option)
%
%function points_rotated = rotate_points(points,angle,rotation_centre,option)
%
%Function to rotate points using rotation matrix
%Created by Mitch Harley


%Set defaults
if nargin==3
    option = 'rads'; %default option is angle in radians
end

%First subtract rotation_center
points_new = points(:,1:2) - repmat(rotation_centre,size(points,1),1);

%Now convert angle to radians
switch option
    case 'rads'
        angle = angle;
    case 'degs'
        angle = angle*pi/180;
    case 'polyfit'
        angle = atan(angle(1));
end

%Now rotate the points
points_rotated = points_new*[cos(angle) -sin(angle); sin(angle) cos(angle)];

%Add elevation column if necessary
if size(points,2)==3
    points_rotated = [points_rotated points(:,3)];
end