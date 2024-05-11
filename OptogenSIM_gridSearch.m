function OptogenSIM
% Version of main OptogenSIM.m without any GUI functionality 
% for simplification of grid search purposes

clear all; home;

if (~isdeployed)
    addpath('./SliceBrowser');
    addpath('.');
end

%%
COMP = ismac;  % 0 = Windows, 1 = Mac
%% locate the mc_VOL.bin, mc_T.bin,mc_default.mci directory

[ctfpath VOLname] = fileparts(which('mc_VOL.bin'));

librdir = fullfile(pwd,'OGSworking');
if ~exist(librdir,'dir')
  mkdir(librdir)
end


if ~exist(fullfile(librdir,'mc.mci'),'file')
    
    copyfile(fullfile(ctfpath,'mc_default.mci'), fullfile(librdir,'mc.mci'))
 
end

if ~exist(fullfile(librdir,'mc_VOL.bin'),'file')
    
    copyfile(fullfile(ctfpath,'mc_VOL.bin'), fullfile(librdir,'mc_VOL.bin'))
 
end

if ~exist(fullfile(librdir,'mc_T.bin'),'file')
    
    copyfile(fullfile(ctfpath,'mc_T.bin'), fullfile(librdir,'mc_T.bin'))
 
end

if COMP
    if ~exist(fullfile(librdir,'gomcxyzOGS'),'file')
        copyfile(fullfile(ctfpath,'gomcxyzOGS'), fullfile(librdir,'gomcxyzOGS'))
    end
    
else
    if ~exist(fullfile(librdir,'gomcxyzOGS.exe'),'file')
        copyfile(fullfile(ctfpath,'gomcxyzOGS.exe'), fullfile(librdir,'gomcxyzOGS.exe'))
    end
    
end


%% Default parameters
timeFLAG = 1;      % 1: specify simulation time or photons; 2: specify the number of photons in simulation 
simFLAG     = 1  ; % 1: default or previous simulation parameters; 0: updated simulations
flagOUTGO   = 1; % 1 = on, 0 = off
conVALUE = 5;  % 5mw/mm-2, YL: default value for the specified contour, 
conFLAG = 0; % YL: 0: use built-in four contours; 1: use the specified contour
VOL2 = [];  % YL
T   = [];
VOL = [];
Nx  = [];
Ny  = [];
Nz = [];
dx = [];
dy = [];
dz = [];
xs = [];
ys = [];
zs = [];
mcflag = [];  
mcflagtxt = {'De-focused Gaussian(1/e^2 radius)','De-focused Gaussian(1/e radius)','Defocused flat beam',...
    'Isotropic point source','Focused Gaussian(1/e^2 radius)','Focused Gaussian(1/e radius)','Focused flat beam'};
currentmcflag = [];
originalname ='';

SIMname = 'mc';
[T VOL H SIMname] = loadmc(SIMname,librdir); % Loads last mc simulation (mc.mci, mc_T.bin, mc_VOL.bin)
if H(16) ==  1.0000e+12
    H(16) = inf;
end
time_photons = H(1); % time or photons  of simulation
Nx      = H(2);
Ny      = H(3);
Nz      = H(4);
dx      = H(5);
dy      = H(6);
dz      = H(7);
mcflag  = H(8);   %0: focused flat; 1:focused Gaussian 1/e; 2: isotropic point; 3:defocused flat ...
                  %4: defocused Gaussian 1/e; 5: focused Gaussian 1/e2; 6: de-focused Gaussian 1/e2;
launchflag = H(9);
boundaryflag = H(10);
xs      = H(11);
ys      = H(12);
zs      = H(13);
xfocus  = H(14);
yfocus  = H(15);
zfocus  = H(16);
ux0     = H(17);
uy0     = H(18);
uz0     = H(19);
radius  = H(20);
waist   = H(21);
nm = H(22);    % light wavelength
pwr     = H(23);       % power of the light source
NA =     H(24);        % numeric aperture of the fiber
thmax =  H(25);
Nt      = H(26);
timeFLAG = H(27);      % 1: specify simulation time or photons; 2: specify the number of photons in simulation 

tissue      = struct('name',[],'mua',[],'mus',[],'g',[]);
tissue      = makeTissueList_OGS(nm);

%% Set beam choice

H(8) = 6;     % 6 is "De-focused Gaussian(1/e^2 radius)"

%% Set contour
% Using default contours
% Go to original source code lines 477-490 to see how to use specified
% contour

defaultcontours = sprintf('%4s %4s %4s %4s','0.1', '1','10','100');
conFLAG = 0;

%% Set SIMname

SIMname = 'mc';     % Default for OptogenSIM

%% Power

pwr = 10;     % mW
H(23) = pwr;

%% Set time

H(1) = 1;     % Time in minutes
H(27) = 1;    % Flag uses time instead of photons

%% Set fiber radius

H(20) = 0.01;    % Fiber radius in cm

%% Fiber NA

NA = 0;
H(24) = NA;
H(25) = asin(NA/1.36);

%% XS
% Change default numbers later
% Figure out how to run only over grey matter

xs = 0.0164;
H(11) = xs;    % Source X in cm

%% YS

ys = -0.1434;
H(12) = ys;   % Source Y in cm

%% ZS

zs = 0.2891;
H(13) = zs;    % Source X in cm


%% Section of code that grid search will loop over to save outputs

% initialize starting wavelength value
wavelength = 300;

for i = 1:3     % change 3 to number of times you want to run
    %Set wavelength
    H(22) = wavelength;

    %% Load mc_H.mci and mc_T.bin
    %  	 Assumes mc_T.bin and mc_H.mci exist. If not, quits.
    %       exist() yields 2 if file exists. 0 if not.
    flagH = (exist(fullfile(librdir,'mc.mci'))>0);
    flagT = (exist(fullfile(librdir,'mc_T.bin'))>0);
    flagV = (exist(fullfile(librdir,'mc_VOL.bin'))>0);
    flagF = (exist(fullfile(librdir,'mc_F.bin'))>0);
    if ~flagH, disp('missing mc.mci');end
    if ~flagT, disp('missing mc_T.bin. Load a simulation.');end
    if ~flagV, disp('missing mc_VOL.bin. Uses mc_T.bin instead.');end
    if ~flagF,
        disp('missing mc_F.bin. Be sure to run GO before using OUTPUT.');
        set(callGO,'backgroundColor','r','String',sprintf('G O\nnot yet run'))
        set(callOUTPUT,'Enable','off')
    end

    %% Save .mci
    saveHmci('mc', librdir, H, nm, tissue,originalname); % create .mci file
    c = round(clock);
    switch COMP
        case 1 % 'Mac'
            bkupname = sprintf('mc_%d-%d_%d_%d_%d',c(4),c(5),c(2),c(3),c(1))% SLJ
        case 0 % 'Win'
            bkupname = sprintf('mc_%d-%d_%d_%d_%d',c(4),c(5),c(2),c(3),c(1)) %
    end
    saveHmci(bkupname, librdir, H,nm, tissue,originalname); % <-------- create .mci file -------

    %% Output

    % Load Fluence rate F(y, x, z)
    mcFname = fullfile(librdir,'mc_F.bin');
    tic
    fid = fopen(mcFname, 'rb');
    [Data count] = fread(fid, Ny*Nx*Nz, 'float');
    fclose(fid);
    toc
    F = reshape(Data,Ny,Nx,Nz);
    bj = 3; % choose to view fluence contours

    % YL: add the absolute threshold for the desired light power density
    [ix,iy,iz] = lookOGSoutput(F,VOL,T,bj,conFLAG,conVALUE,librdir);

    newdir = librdir + "\" + gridsearch_nm; % Create new dir with nm appended
    mkdir(newdir);
    [ix,iy,iz] = lookOGSoutput(F,VOL,T,bj,conFLAG,conVALUE,newdir); % Save output under new dir

    
    wavelength = wavelength + 5;    % change 5 to desired wavelength increment
end

end
