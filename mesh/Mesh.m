%% Panair mesh creation
%
% Adrien Crovato
% ULg 2015-2016
%
% This scripts reads airfoil data for panair/aero wing mesh creation

%% Initilization

clear all; 
close all;

% Folders and paths
nAirf = 3; % number of airfoil used to generate the wing
fname = cell(1,nAirf); % name of the airfoils data file
fpath = ''; % folder where airfoils data are
fname{1,1}(1,:) = 'nasaSC20714';
fname{1,2}(1,:) = 'nasaSC20712';
fname{1,3}(1,:) = 'nasaSC20712';
outFname = 'mesh.dat'; % output file name

file = cell(1,nAirf); % path to file
for i=1:nAirf
    file{1,i}(1,:) = [fpath fname{1,i}(1,:) '.dat'];
end
meshGen = 2; % Output type (0=none, 1=Panair, 2=Aero)

% Wing geometry (multi planform defintion for HALF wing)
meanC4AirfSweep = 0*pi/180; % mean quarter-chord sweep (0 if airfoil coordinates are provided STREAMWISE, mean sweep angle if CHORDWISE)
span = [2; 4]; % panform span
taper = [0.6; 0.4]; % planform taper
sweepLE = [25; 25]*pi/180; % leading edge planform sweep
% sweepLE = atan(tan(sweepC4)+(1-taper)./((AR)*(1+taper))); %if quarter chord sweep provided instead
dihedral = [5; 2]*pi/180; % planform dihedral
twist = [0.; 0.; -1.]*pi/180; % planform twist
rootChord = 2; % root chord length

chord = zeros(nAirf,1);
area = zeros(nAirf-1,1);
chord(1) = rootChord;
for i =2:nAirf
    chord(i) = chord(i-1)*taper(i-1);
    area(i-1) = (chord(i-1)+chord(i))*span(i-1)/2;
end
spanPos = zeros(nAirf,1);
for i=2:nAirf
    spanPos(i) = sum(span(1:i-1)); % spanwise airfoil stations
end
S = sum(area); % half wing area
b = sum(span); % half wing span
AR = 2*b^2/S; % wing aspect ratio

% Mesh parameters
meshType = 1; % mesh distribution type (1=half-cosine, 2=linear, 3=geometric)
prog = 1.1; % if mesh type is geometric, defines progression
nSpanStat = [5; 11]; % number of spanwise stations (spanwise panels - 1) for EACH planform
nPoints = 101; % number of chordwise points (chordwise panels - 1), MUST BE ODD AND DIVISIBLE BY 4 when 1 is substracted

%% Scan

data = cell(1,nAirf);
dataBot = cell(1,nAirf);
dataTop = cell(1,nAirf); % airfoils coordinates container

for j=1:nAirf
    fid = fopen(file{1,j}(1,:)); % open file
    header = fgetl(fid); % read header
    dataTmp = fscanf(fid,'%g %g', [2 inf]); % scan data
    fid =fclose(fid); % close file
    dataTmp = dataTmp'; % transpose to correct format
    nBasePoint = size(dataTmp,1); % total number of points

    % Store coordinates
    data{1,j} = [dataTmp(:,1) spanPos(j)*ones(size(dataTmp,1),1) dataTmp(:,2)];
    dataBot{1,j} = [dataTmp((nBasePoint+3)/2:end,1) spanPos(j)*ones((size(dataTmp,1)-1)/2,1) dataTmp((nBasePoint+3)/2:end,2)];
    dataTop{1,j} = [dataTmp(1:(nBasePoint-1)/2,1) spanPos(j)*ones((size(dataTmp,1)-1)/2,1) dataTmp(1:(nBasePoint-1)/2,2)];
    clear dataTmp;
end

%% Wing creation

% Airfoil stretching (if LE perpedicular airfoil provided)
if (meanC4AirfSweep ~= 0)
    for i = 1:nAirf
        data{1,i}(:,3) = data{1,i}(:,3)*cos(meanC4AirfSweep);
        dataBot{1,i}(:,3) = dataBot{1,i}(:,3)*cos(meanC4AirfSweep);
        dataTop{1,i}(:,3) = dataTop{1,i}(:,3)*cos(meanC4AirfSweep);
    end
end

% Chord scaling
for i = 1:nAirf
    data{1,i}(:,[1 3]) = data{1,i}(:,[1 3])*chord(i);
    dataBot{1,i}(:,[1 3]) = dataBot{1,i}(:,[1 3])*chord(i);
    dataTop{1,i}(:,[1 3]) = dataTop{1,i}(:,[1 3])*chord(i);
end

% Twist rotation
for i = 1:nAirf
    data{1,i}(:,[1 3]) = data{1,i}(:,[1 3])*[cos(twist(i)) -sin(twist(i)); sin(twist(i)) cos(twist(i))];
    dataBot{1,i}(:,[1 3]) = dataBot{1,i}(:,[1 3])*[cos(twist(i)) -sin(twist(i)); sin(twist(i)) cos(twist(i))];
    dataTop{1,i}(:,[1 3]) = dataTop{1,i}(:,[1 3])*[cos(twist(i)) -sin(twist(i)); sin(twist(i)) cos(twist(i))];
end

% LE sweep translation
for i = 2:nAirf
    data{1,i}(:,1) = data{1,i}(:,1) + min(data{1,i-1}(:,1)) + tan(sweepLE(i-1))*span(i-1);
    dataBot{1,i}(:,1) = dataBot{1,i}(:,1) + min(data{1,i-1}(:,1)) + tan(sweepLE(i-1))*span(i-1);
    dataTop{1,i}(:,1) = dataTop{1,i}(:,1) + min(data{1,i-1}(:,1)) + tan(sweepLE(i-1))*span(i-1);
end

% Dihedral transaltion
for i = 2:nAirf
    data{1,i}(:,3) = data{1,i}(:,3) + sum(tan(dihedral(1:i-1)).*(span(1:i-1)));
    dataBot{1,i}(:,3) = dataBot{1,i}(:,3) + sum(tan(dihedral(1:i-1)).*(span(1:i-1)));
    dataTop{1,i}(:,3) = dataTop{1,i}(:,3) + sum(tan(dihedral(1:i-1)).*(span(1:i-1)));
end

% Display planform
figure
hold on
set(gca,'YDir','reverse');
for i = 1:nAirf
    plot3([dataBot{1,i}(1,2) dataBot{1,i}(end,2)],[dataBot{1,i}(1,1) dataBot{1,i}(end,1)],[dataBot{1,i}(1,3) dataBot{1,i}(end,3)],'Linewidth', 2, 'Color', 'b');
end
for i = 1:nAirf-1
    plot3([dataBot{1,i}(1,2) dataBot{1,i+1}(1,2)],[dataBot{1,i}(1,1) dataBot{1,i+1}(1,1)],[dataBot{1,i}(1,3) dataBot{1,i+1}(1,3)],'Linewidth', 2, 'Color', 'b');
    plot3([dataBot{1,i}(end,2) dataBot{1,i+1}(end,2)],[dataBot{1,i}(end,1) dataBot{1,i+1}(end,1)],[dataBot{1,i}(end,3) dataBot{1,i+1}(end,3)],'Linewidth', 2, 'Color', 'b');
end
xlabel('$y$', 'Fontsize', 16, 'Interpreter', 'Latex');
ylabel('$x$', 'Fontsize', 16, 'Interpreter', 'Latex');
zlabel('$z$', 'Fontsize', 16, 'Interpreter', 'Latex');
title('Wing panform', 'Fontsize', 20, 'Interpreter', 'Latex');

% Display airfoils
figure
hold on
view([0,1,0])
for i=1:nAirf
    plot3(data{1,i}(:,1),data{1,i}(:,2),data{1,i}(:,3),'Linewidth', 2, 'Marker', 'o');
end
xlabel('$x$', 'Fontsize', 16, 'Interpreter', 'Latex');
ylabel('$y$', 'Fontsize', 16, 'Interpreter', 'Latex');
zlabel('$z$', 'Fontsize', 16, 'Interpreter', 'Latex');
title('Wing airfoils', 'Fontsize', 20, 'Interpreter', 'Latex');

%% Mesh creation

interBot = zeros((nPoints-1)/2,3,sum(nSpanStat));
interTop = zeros((nPoints-1)/2,3,sum(nSpanStat));
inter = zeros(nPoints,3,sum(nSpanStat));

interSpanPos = zeros(sum(nSpanStat),1);
index = 1;
for k=1:length(nSpanStat)
    % Min/max value
    minV = zeros(1,2);
    minI = zeros(1,2);
    [minV(1),minI(1)] = min(data{1,k}(:,1));
    [minV(2),minI(2)] = min(data{1,k+1}(:,1));
    % Spanwise location
    interSpanPos(index:index+nSpanStat(k)-1,1) = linspace(spanPos(k),spanPos(k+1),nSpanStat(k))';
    % Chordwise interpolation between airfoils k & k+1
    if meshType == 1
        beta = linspace(0.0175,pi,(nPoints-1)/2)';
        interXTmp = 0.5*(1-cos(beta)); % Half cosine based spacing
    elseif meshType == 3
        a = 0.5*(1-prog)/(1-prog^((nPoints-1)/4));
        dummy = zeros((nPoints-1)/4,1);
        dummy(1) = a;
        for i = 2:(nPoints-1)/4
            dummy(i) = dummy(i-1)+a*prog^(i-1);
        end
        interXTmp = [dummy;[1-flipud(dummy(1:end-1));1]]; % geometric spacing
    else
        interXTmp = linspace(0.001,1,(nPoints-1)/2)'; % linear spacing
    end
    interX1Tmp = interXTmp*(max(data{1,k}(:,1))-minV(1)) + minV(1);
    interX2Tmp = interXTmp*(max(data{1,k+1}(:,1))-minV(2)) + minV(2);    
    interBotY1Tmp = interp1(dataBot{1,k}(:,1),dataBot{1,k}(:,3),interX1Tmp,'spline','extrap');
    interTopY1Tmp = interp1(dataTop{1,k}(:,1),dataTop{1,k}(:,3),interX1Tmp,'spline','extrap');
    interBotY2Tmp = interp1(dataBot{1,k+1}(:,1),dataBot{1,k+1}(:,3),interX2Tmp,'spline','extrap');
    interTopY2Tmp = interp1(dataTop{1,k+1}(:,1),dataTop{1,k+1}(:,3),interX2Tmp,'spline','extrap');
    
    for i = 1:nSpanStat(k)
        coeff = [1-((interSpanPos(index+i-1)-spanPos(k))/(spanPos(k+1)-spanPos(k))) (interSpanPos(index+i-1)-spanPos(k))/(spanPos(k+1)-spanPos(k))];
        interBot(:,1,index+i-1) = coeff(1)*interX1Tmp + coeff(2)*interX2Tmp;
        interBot(:,2,index+i-1) = interSpanPos(index+i-1)*ones((nPoints-1)/2,1);
        interBot(:,3,index+i-1) = coeff(1)*interBotY1Tmp + coeff(2)*interBotY2Tmp;
        interTop(:,1,index+i-1) = coeff(1)*interX1Tmp + coeff(2)*interX2Tmp;
        interTop(:,2,index+i-1) = interSpanPos(index+i-1)*ones((nPoints-1)/2,1);
        interTop(:,3,index+i-1) = coeff(1)*interTopY1Tmp + coeff(2)*interTopY2Tmp;
        
        inter(:,:,index+i-1) = [flipud(interTop(:,:,index+i-1));[coeff(1)*min(data{1,k}(:,1))+coeff(2)*min(data{1,k+1}(:,1)) interSpanPos(index+i-1) coeff(1)*data{1,k}(minI(1),3)+coeff(2)*data{1,k+1}(minI(2),3)];interBot(:,:,index+i-1)];
    end
    index = sum(nSpanStat(1:k))+1;
end

% Delete intermediate (doublons) stations
for k=2:length(nSpanStat)
    rem = sum(nSpanStat(1:k-1))+3-k;
    inter(:,:,rem) = [];
    interBot(:,:,rem) = [];
    interTop(:,:,rem) = [];
end

%% Printout

if meshGen == 1
    
    % Wing
    fid = fopen(outFname,'w');
    fprintf(fid,'%9.0f %9.0f \n\n',nPoints,sum(nSpanStat)-length(nSpanStat)+1);
    for i=1:sum(nSpanStat)-length(nSpanStat)+1
        for j=2:2:nPoints-1
            fprintf(fid,'%9.6f %9.6f %9.6f',inter(j-1,1,i),inter(j-1,2,i),inter(j-1,3,i));
            fprintf(fid,' %9.6f %9.6f %9.6f \n',inter(j,1,i),inter(j,2,i),inter(j,3,i));
        end
        fprintf(fid,'%9.6f %9.6f %9.6f \n',inter(nPoints,1,i),inter(nPoints,2,i),inter(nPoints,3,i));
    end
    fprintf(fid,'\n');
    
    % Tip (may require modifications)
    disp('WARNING: check points (especially LE) if you plan to use tip newtork!');
    fprintf(fid,'%9.0f %9.0f \n\n',(nPoints+1)/2,2);
    for j=(nPoints-1)/2:-2:2
        fprintf(fid,'%9.6f %9.6f %9.6f',interTop(j,1,i),interTop(j,2,i),interTop(j,3,i));
        fprintf(fid,' %9.6f %9.6f %9.6f \n',interTop(j-1,1,i),interTop(j-1,2,i),interTop(j-1,3,i));
    end
    fprintf(fid,'%9.6f %9.6f %9.6f \n',inter((nPoints+1)/2,1,i),inter((nPoints+1)/2,2,i),inter((nPoints+1)/2,3,i));
    fprintf(fid,'\n');
    for j=(nPoints-1)/2:-2:2
        fprintf(fid,'%9.6f %9.6f %9.6f',interBot(j,1,i),interBot(j,2,i),interBot(j,3,i));
        fprintf(fid,' %9.6f %9.6f %9.6f \n',interBot(j-1,1,i),interBot(j-1,2,i),interBot(j-1,3,i));
    end
    fprintf(fid,'%9.6f %9.6f %9.6f \n',inter((nPoints+1)/2,1,i),inter((nPoints+1)/2,2,i),inter((nPoints+1)/2,3,i));
    fid = fclose(fid);
    
    disp('');
    disp('Panair mesh file generated!')
    disp('');

elseif meshGen==2
    
    % Wing
    fid = fopen(outFname,'w');
    fprintf(fid,'Mesh File\n');
    fprintf(fid,'$size\n');
    fprintf(fid,'%9.0f %9.0f\n',nPoints,sum(nSpanStat)-length(nSpanStat)+1);
    fprintf(fid,'$points\n');
    for i=1:sum(nSpanStat)-length(nSpanStat)+1
        for j=1:nPoints
            fprintf(fid,' %9.6f %9.6f %9.6f \n',inter(j,1,i),inter(j,2,i),inter(j,3,i));
        end
    end
    fid = fclose(fid);
    
    disp('');
    disp('FPM mesh file generated!')
    disp('');
else
    disp('');
    disp('No mesh file generated!')
    disp('');
end

%% Display

figure
hold on
% axis([0,1,-1.5,1])
% set(gca,'XTick',0:0.1:1);
% set(gca,'YTick',-1.5:0.5:1);
% set(gca,'YDir','reverse');
view([0,1,0])

for i=1:sum(nSpanStat)-length(nSpanStat)+1
    plot3(inter(:,1,i),inter(:,2,i),inter(:,3,i),'Linewidth', 2, 'Marker', 'o');
end
xlabel('$x$', 'Fontsize', 16, 'Interpreter', 'Latex');
ylabel('$y$', 'Fontsize', 16, 'Interpreter', 'Latex');
zlabel('$z$', 'Fontsize', 16, 'Interpreter', 'Latex');
title('Mesh Points', 'Fontsize', 20, 'Interpreter', 'Latex');