function points = unrotate_points(points_rotated,angle,rotation_centre,option)

%Set defaults
if nargin==3
    option = 'rads'; %default option is angle in radians
end

points_rotated_new = points_rotated(:,1:2);

%Convert angle to radians
switch option
    case 'rads'
        angle = angle;
    case 'degs'
        angle = angle*pi/180;
    case 'polyfit'
        angle = atan(angle(1));
end

%Now rotate the points back
points_new = points_rotated_new*[cos(-angle) -sin(-angle); sin(-angle) cos(-angle)];

%Add rotation_center back to points
points = points_new + repmat(rotation_centre,size(points_new,1),1);

%Add elevation column if necessary
if size(points_rotated,2)==3
    points = [points points_rotated(:,3)];
end