function [ embryo ] = csimmath_placeCrystal( simParams, varargin )
%%% Written by Jonathan K. Meyers (ORCID 0000-0002-6698-3420)
%%% Takes the simulation parameters
%%% Creates an embryo and puts it somewhere in the simulation

%%%Optional arguments can override position/location
% 'Cx', num in the coordinates of the simulation
% 'Cy', num in the coordinates of the simulation
% 'rot', num between (0,90)

embryo = struct; %initialize
embryo.id = simParams.crystalCount + 1; %get new ID
embryo.timeStep = simParams.timeStep; %this has been incremented already


%get random position and rotation
embryo.Cx = rand() * ...
    simParams.simZoneSize_width - simParams.simZoneSize_width/2;
embryo.Cy = rand() * ...
    simParams.simZoneSize_width - simParams.simZoneSize_width/2;
embryo.rot = rand() * 90; %excludes 0 and 90
% (including 0 or 90 would make it a little more complicated)

%mostly for testing purposes, you can override the position/rotation
for vv = 1:length(varargin)
    if strcmpi(varargin{vv},'Cx')
        embryo.Cx = varargin{vv+1};
    elseif strcmpi(varargin{vv},'Cy')
        embryo.Cy = varargin{vv+1};
    elseif strcmpi(varargin{vv},'rot')
        embryo.rot = varargin{vv+1};
    end
end


%size dictated by simulation embryo size 
embryo.L = simParams.embryoSize;

%calculate the slopes based on the angle provided
embryo.m1 = tand(embryo.rot);
%perpendicular line has slope of negative reciprocal of the first
embryo.m2 = -1/embryo.m1;

%find intercepts of all lines using some geometry and trigonometry
embryo.b11 = embryo.Cy - embryo.m1 * embryo.Cx - embryo.L/2 * ...
    ( embryo.m1 * sind(embryo.rot) + cosd(embryo.rot) );
embryo.b12 = embryo.Cy - embryo.m1 * embryo.Cx + embryo.L/2 * ...
    ( embryo.m1 * sind(embryo.rot) + cosd(embryo.rot) );
embryo.b21 = embryo.Cy - embryo.m2 * embryo.Cx - embryo.L/2 * ...
    ( embryo.m2 * cosd(embryo.rot) - sind(embryo.rot) );
embryo.b22 = embryo.Cy - embryo.m2 * embryo.Cx + embryo.L/2 * ...
    ( embryo.m2 * cosd(embryo.rot) - sind(embryo.rot) );


%find the corners of the square and assign the X limits of each edge
slopeDif = embryo.m1 - embryo.m2;
embryo.x11 = sort([...
    (embryo.b21 - embryo.b11), (embryo.b22 - embryo.b11)] ./ slopeDif);
embryo.x12 = sort([...
    (embryo.b21 - embryo.b12), (embryo.b22 - embryo.b12)] ./ slopeDif);
embryo.x21 = sort([...
    (embryo.b21 - embryo.b11), (embryo.b21 - embryo.b12)] ./ slopeDif);
embryo.x22 = sort([...
    (embryo.b22 - embryo.b11), (embryo.b22 - embryo.b12)] ./ slopeDif);


%record the vertex positions
embryo.xcoords = [...
    (embryo.b21 - embryo.b11),...
    (embryo.b22 - embryo.b11),...
    (embryo.b22 - embryo.b12),...
    (embryo.b21 - embryo.b12)...
    ] ./ slopeDif;
embryo.ycoords = embryo.m1 .* embryo.xcoords + [...
    embryo.b11,...
    embryo.b11,...
    embryo.b12,...
    embryo.b12...
    ];


%initialize these in the struct. Don't know them.
embryo.placeAttempts = [];
embryo.active = [];
embryo.L_final = [];
embryo.xcoords_final = [];
embryo.ycoords_final = [];


end

