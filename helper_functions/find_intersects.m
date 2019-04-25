function [mStructs, numM, shitlist,ylistr] = find_intersects(lines, xb,yb, raylength,maxsampleintersects)
% function [mStructs, numM, shitlist,ylistr] = find_intersects2(lines, xb,yb, raylength,maxsampleintersects)
% OR:
% [mStructs, numM, shitlist,ylistr] = find_intersects2(lines, xb,yb)
% Johannes Hendriks, 2018
% The University of Newcastle, Australia
%
% Inputs:
%   xb: the x positions of the boundary, use NaNs to represent non
%   continuous
%   yb: the y positions of the boundary, use NaNs to represent non
%   continuous
%   lines: The ray information
%   raylength (Optional): long enough to go in and out of the sample
%   maxsampleintersects (optional): to check that not too many intersects
%   were found
% 
% Outputs:
%   mStructs:
%   numM: number of successful measurements
%   shitlist: the list of rays that didnt intersect successfully
%   ylistr: the list of rays that did intersect successfully
%
% This function takes in ray and sample data and determines intersect
% locations. Stores data into an mStructs structure. 
%
% This function utilises the inbuilt function polyxpoly in order to
% vectorise the finding of intersects
%
% This code can handle non-convex samples with multiple in/outs per ray.
%
% Version history
%   parFindIntersects (Johannes Hendriks, 2015).
%   as_parFindIntersects (Alex Gregg, 2016).
%   vecFindItersects (Johannes Hendriks, 2015).
%   vecFindItersects_multi (Alex Gregg, 2016).
%   find_intersects (Johannes Hendriks, 2018)

%% Intial data re-organisation for intersect finding
% Determine how many rays (successful or unsuccessful) we have.
[~,~,numRays] = size(lines);
ylist = 1:numRays;


% To detect intersections, we must extend the ray a distance out from the
% pixel. If this length is not specified, calculate it.
if nargin < 4      % required ray length not provided
    d = max(sqrt(x_verts.^2 + y_verts.^2)); % max distance from origin ('defined center of pixel grid') to a vertex
    raylength = 2*d*1.5;                   % twice radius like distance and then factor of safety
end

% For validation, we specify a maximum number of intersects, if not, assume
% 4
if nargin < 5
    warning('No maximum number of intersects specified... using 4.')
    maxsampleintersects = 4;
    % Update this if your sample has the possibility of more 
end


% Build vectors X and Y for the rays in a form suitable for polyxpoly;
% Begin by calculating nHat for every ray
nhats = squeeze(lines(:,2,:) - lines(:,1,:));
nhats = nhats./repmat(sqrt(sum(nhats.^2)),2,1);
% build X and Y vectors for the rays
X = [squeeze(lines(1,1,:))' - nhats(1,:)*raylength/2;                       % x coord of start of segment
    squeeze(lines(1,1,:))' + nhats(1,:)*raylength/2;                        % x coord of end of segment
    NaN(1,numRays)];                                                        % breaks segments up
Y = [squeeze(lines(2,1,:))'- nhats(2,:)*raylength/2;                        % y coord of start of segment
    squeeze(lines(2,1,:))'+ nhats(2,:)*raylength/2;                         % y coord of end of segment
    NaN(1,numRays)];                                                        % breaks segments up
X = X(:);                                                                   % reorders into single line where the NaNs represent breaks
Y = Y(:);                       

% Detect Intersections with Polyxpoly and remove NaN segments from index
% numbering. xi is the unordered x co-ordinates of each intersect, yi is
% the corresponding y co-ordinates of each intersect, ii is the face-ray
% data for each intersect (col 1 is ray no, col 2 is face no), each row is
% an intersect).
[xi,yi,ii] = polyxpoly(X(:),Y(:),xb,yb);

ii(:,1) = (ii(:,1)+2)/3;        % remove the NaN segments from the index numbering


% Re-order the nhats and pixels so that they are in the same order as
% intersects.

nhati = nhats(:,ii(:,1))';       % direction of ray that created each intersect
pixeli = squeeze(lines(:,1,ii(:,1)));   % pixel location that created each intersect
pixeli = transpose(pixeli);

%% Sort by Ray
% Sort ii by first column,  extract indices (sort by ray), extract unsorted ray inds.
rayinds = ii(:,1);
[rayinds,Is] = sort(rayinds,'ascend');
ii = ii(Is,:);
faceinds = ii(:,2);
% Sort x,y,nhat,pixel the same way.
xi = xi(Is);
yi = yi(Is);
nhati = nhati(Is,:);
pixeli = pixeli(Is,:);

%% Calculate maximum number of intersects per ray

% Determine num of intersects per ray. Each row corresponds to a ray.
[urayinds,~,IC] = unique(rayinds);

% determine the number of intersects per ray
counts = hist(rayinds,[1:numel(lines(1,1,:))]);

if any(counts==0)
%     warning('Some lines did not give rise to intersects by polyxpoly, indices stored.')
    % Find the indices of these unsucessful lines
    shitlist = find(~counts);
%     display(['Bad indices: ',num2str(shitlist)])
else
    shitlist = [];
end

% Delete shit ones from shitlist
ylistr = ylist;
ylistr(shitlist) = [];

% [intsperray, ~] = hist(rayinds,urayinds);

% counts per successful ray
intsperray = counts(urayinds);
intsperray = intsperray';

% Delete any 0 intersect rays and delete any >max intersect rays
% Cond 1: If polyxpoly fucks up and lets through an unsuccesful ray
% Cond 2: If we dont separate the rays correctly
% Needs to be tested (case 2)....


del = or(or(intsperray==0,intsperray>maxsampleintersects), mod(intsperray,2));
if any(del)
    delrays = del(IC,:);
%     warning('0 or >max intersect rays detected. Deleting...')
    xi(delrays,:) = [];
    yi(delrays,:) = [];
    nhati(delrays,:) = [];
    pixeli(delrays,:) = [];
    intsperray(del,:) = [];     % different because already on a ray basis
    ylistr(del) = [];
    rayinds(delrays,:) = [];
end

%remove more


%% Sort within rays by intersect order
% Create a matrix of x,y co-ords
Xis = [xi, yi];

% calculate the distance of each intersect from centre c
c = dot(Xis,nhati,2);

% distribute the c, x and y into temp matrices
% each row is a ray, each column is an intersect
% size is numrays x maximum num intersects
% where a ray has less than the max number of intersects, NaNs are used.
% Also get nhat on a per ray basis
% First pre-allocate
tempc = NaN(numel(intsperray),max(intsperray));
tempx = NaN(numel(intsperray),max(intsperray));
tempy = NaN(numel(intsperray),max(intsperray));
tempfaces = NaN(numel(intsperray),max(intsperray));

% Get only one nhat for each measurement -

[check,Irayinds,~] = unique(rayinds,'rows','stable');
    
nhatr = nhati(Irayinds,:);
pixelr = pixeli(Irayinds,:);

% counter for distributing values...
count = 0;
for i = 1:numel(intsperray)
    tempc(i,1:intsperray(i)) = c((1+count):intsperray(i)+count);
    tempx(i,1:intsperray(i)) = Xis((1+count):intsperray(i)+count,1);
    tempy(i,1:intsperray(i)) = Xis((1+count):intsperray(i)+count,2);
    tempfaces(i,1:intsperray(i)) = faceinds((1+count):intsperray(i)+count);
    count = count + intsperray(i);
end


% Sort in terms of c value distance from center...
[tempcs, Ics] = sort(tempc,2,'ascend');

% Play with the indices to re-order the other temp matrices the same way
% ix2 is the index along dimension 2, and we want dimension 1 to remain unchanged
Ics1 = repmat([1:size(tempc,1)]', [1 size(tempc,2)]); %//'
% Convert to linear index equivalent of the reordering of the sort() call
Ics2 = sub2ind(size(tempc), Ics1, Ics);

% Sort the other matrices the same way (re-order columns so in intersect
% order)
tempx = tempx(Ics2);
tempy = tempy(Ics2);




% Determine how many successful measurements were made in total.
[numM, maxintersections] = size(tempx);


%% Calculate remaining information
% Ray s values (shift distance from center by first value)
temps = tempcs - repmat(tempcs(:,1),1,maxintersections);

% Entry and exit intersect co-ordinates
xitemp = NaN(2*numM,maxintersections/2);
xjtemp = NaN(2*numM,maxintersections/2);

% Entry
xitemp(1:2:end,:) = tempx(:,1:2:end);
xitemp(2:2:end,:) = tempy(:,1:2:end);
% Exit
xjtemp(1:2:end,:) = tempx(:,2:2:end);
xjtemp(2:2:end,:) = tempy(:,2:2:end);

% Irradiated Length
tempst = temps;
tempst(isnan(temps))=0;

L = (sum(tempst(:,2:2:end),2)-sum(tempst(:,1:2:end),2));
%% Data Storage - Initialise, replicate the structure then populate.

% Initialise the measurement structure
mStruct.exit = NaN(2,maxintersections);   % Exit intersections. Row 1 x, Row 2 y, each column a segment.
mStruct.entry = NaN(2,maxintersections);   % Entry intersections. Row 1 x, Row 2 y, each column a segment.
mStruct.nHat = NaN(2,1);                  % unit normal of measurement ray    
mStruct.pixel = NaN(2,1);               % position of measurement pixel
mStruct.L = NaN;                        % Total irradiated ray length
mStruct.y = NaN;                        % measurement value (NOT CALCULATED HERE)
mStruct.nsegs = NaN;
mStructs = repmat(mStruct,numM,1);          % replicate the structure

% Store xiDef
xifDefcell = mat2cell(xitemp,[2*ones(1,numM)],[maxintersections/2]);
[mStructs.entry] = xifDefcell{:};

% Store xjDef
xjDefcell = mat2cell(xjtemp,[2*ones(1,numM)],[maxintersections/2]);
[mStructs.exit] = xjDefcell{:};

% Store nHat
temp = nhatr';
temp = temp(:);
nHatCells = mat2cell(temp,2*[ones(1,numM)],[1]);
[mStructs.nHat] = nHatCells{:};

% Store Pixel
temp = pixelr';
temp = temp(:);
pixelCells = mat2cell(temp(:),[2*ones(1,numM)],[1]);
[mStructs.pixel] = pixelCells{:};

% Store L
% LCell = mat2cell(L,[ones(1,numM)],[1]);
LCell = num2cell(L);
[mStructs.L] = LCell{:};

% Store nsegs
% nsegsCells = mat2cell(intsperray/2,[ones(1,numM)],[1]);
nsegsCells = num2cell(intsperray/2);
[mStructs.nsegs] = nsegsCells{:};

% store the rayInd of the successful measurement
% riCell = mat2cell(ylistr,[ones(1,numM)],[1]);
riCell = num2cell(ylistr);
[mStructs.rayInd] = riCell{:};

end