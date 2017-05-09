%% Display
%
% Adrien Crovato
% ULg 2016-2017
%
% Visual post-processing for Celia

%% Initilization

clear; 
close all;

% Colors
red = [255,10,10]./255;
darkRed = [200,10,10]./255;
green = [10,140,70]./255;
lightGreen = [10,180,70]./255;
blue = [10,10,255]./255;
lightBlue = [10,100,202]./255;
indigo = [75,0,130]./255;
violet = [148,0,211]./255;
orange = [255,127,0]./255;

% Flow and geometry
Mach = 0.6; % Mach number
AOA = 2; % angle of attack
span = 1.196; % half-span of wing
rootChord = 0.8059; % root chord
taper = 0.56; % taper ratio
sweep = 30*pi/180; % LE sweep
MAC = 0.6461; % mean aerodynamic chord
Cp_star = 1/(0.7*Mach^2)*((1/1.2+(Mach^2)/6)^3.5-1);

% Mesh
n = 100;
m = 10;
x = 41;
y = 7;
z = 19;
boundX = [-0.5, 1.6];
boundY = [0, 1.3];
boundZ = [-0.5, 0.5];

% Path
fpath = ''; % folder where data are
fname1 = 'm6s41X7X19';
fname2 = 'm6f41X7X19';
file1 = [fpath fname1 '.dat'];
file2 = [fpath fname2 '.dat'];

% Tranair Cp (at MAC) for comparison
tranData = [9.999999E-01  3.416026E-01
9.999367E-01  3.420147E-01
9.997509E-01  3.421271E-01
9.994387E-01  3.407735E-01
9.990034E-01  3.390408E-01
9.984422E-01  3.367126E-01
9.977580E-01  3.337867E-01
9.969493E-01  3.304496E-01
9.960176E-01  3.265308E-01
9.949610E-01  3.218942E-01
9.937818E-01  3.161748E-01
9.924793E-01  3.088508E-01
9.910542E-01  2.996938E-01
9.895069E-01  2.886930E-01
9.878374E-01  2.763424E-01
9.860461E-01  2.632296E-01
9.841337E-01  2.497632E-01
9.821007E-01  2.366495E-01
9.799482E-01  2.247064E-01
9.776754E-01  2.141579E-01
9.752836E-01  2.043291E-01
9.727726E-01  1.946466E-01
9.701451E-01  1.852735E-01
9.674006E-01  1.763376E-01
9.645384E-01  1.675766E-01
9.615614E-01  1.588385E-01
9.584687E-01  1.502402E-01
9.552616E-01  1.417852E-01
9.519410E-01  1.332972E-01
9.485080E-01  1.247337E-01
9.449647E-01  1.161355E-01
9.413088E-01  1.074623E-01
9.375431E-01  9.872467E-02
9.336687E-01  8.993997E-02
9.296863E-01  8.102377E-02
9.255957E-01  7.187203E-02
9.214009E-01  6.247268E-02
9.170998E-01  5.295017E-02
9.126949E-01  4.359137E-02
9.081879E-01  3.475999E-02
9.035777E-01  2.675039E-02
8.988679E-01  1.959814E-02
8.940591E-01  1.314965E-02
8.891516E-01  7.211460E-03
8.841466E-01  1.558045E-03
8.790467E-01 -3.952852E-03
8.738524E-01 -9.346612E-03
8.685644E-01 -1.458969E-02
8.631847E-01 -1.972690E-02
8.577141E-01 -2.484623E-02
8.521547E-01 -2.994894E-02
8.465075E-01 -3.500606E-02
8.407739E-01 -4.002740E-02
8.349563E-01 -4.503329E-02
8.290539E-01 -4.999242E-02
8.230694E-01 -5.485913E-02
8.170058E-01 -5.965345E-02
8.108621E-01 -6.441033E-02
8.046409E-01 -6.910985E-02
7.983449E-01 -7.375338E-02
7.919738E-01 -7.838273E-02
7.855298E-01 -8.300467E-02
7.790144E-01 -8.761554E-02
7.724298E-01 -9.224444E-02
7.657768E-01 -9.691208E-02
7.590587E-01 -1.016148E-01
7.522743E-01 -1.063613E-01
7.454283E-01 -1.111685E-01
7.385213E-01 -1.160559E-01
7.315542E-01 -1.210318E-01
7.245301E-01 -1.261079E-01
7.174489E-01 -1.313212E-01
7.103145E-01 -1.366730E-01
7.031276E-01 -1.421269E-01
6.958899E-01 -1.476723E-01
6.886029E-01 -1.532997E-01
6.812692E-01 -1.590021E-01
6.738909E-01 -1.648034E-01
6.664689E-01 -1.707017E-01
6.590053E-01 -1.766849E-01
6.515023E-01 -1.827637E-01
6.439617E-01 -1.889207E-01
6.363847E-01 -1.951396E-01
6.287743E-01 -2.014255E-01
6.211311E-01 -2.077634E-01
6.134574E-01 -2.141353E-01
6.057567E-01 -2.205380E-01
5.980287E-01 -2.269547E-01
5.902772E-01 -2.333501E-01
5.825025E-01 -2.397171E-01
5.747062E-01 -2.460490E-01
5.668930E-01 -2.523038E-01
5.590624E-01 -2.584913E-01
5.512169E-01 -2.645923E-01
5.433585E-01 -2.705397E-01
5.354909E-01 -2.763727E-01
5.276126E-01 -2.820712E-01
5.197285E-01 -2.875671E-01
5.118389E-01 -2.929148E-01
5.039466E-01 -2.980718E-01
4.960533E-01 -3.030214E-01
4.881610E-01 -3.077913E-01
4.802714E-01 -3.123044E-01
4.723873E-01 -3.166159E-01
4.645090E-01 -3.206741E-01
4.566403E-01 -3.243284E-01
4.487819E-01 -3.277523E-01
4.409375E-01 -3.310194E-01
4.331069E-01 -3.340861E-01
4.252926E-01 -3.369858E-01
4.174974E-01 -3.396102E-01
4.097227E-01 -3.419818E-01
4.019712E-01 -3.441501E-01
3.942432E-01 -3.459858E-01
3.865414E-01 -3.475905E-01
3.788689E-01 -3.491305E-01
3.712256E-01 -3.506160E-01
3.636152E-01 -3.520716E-01
3.560382E-01 -3.534449E-01
3.484972E-01 -3.546760E-01
3.409935E-01 -3.558329E-01
3.335299E-01 -3.568835E-01
3.261079E-01 -3.577786E-01
3.187296E-01 -3.586146E-01
3.113959E-01 -3.594614E-01
3.041100E-01 -3.603096E-01
2.968723E-01 -3.611830E-01
2.896854E-01 -3.621287E-01
2.825499E-01 -3.631505E-01
2.754698E-01 -3.642376E-01
2.684457E-01 -3.654317E-01
2.614786E-01 -3.667401E-01
2.545705E-01 -3.681314E-01
2.477244E-01 -3.696520E-01
2.409412E-01 -3.713214E-01
2.342220E-01 -3.730867E-01
2.275702E-01 -3.749980E-01
2.209855E-01 -3.771120E-01
2.144701E-01 -3.793347E-01
2.080261E-01 -3.816727E-01
2.016550E-01 -3.842229E-01
1.953579E-01 -3.869677E-01
1.891378E-01 -3.898558E-01
1.829941E-01 -3.928863E-01
1.769294E-01 -3.960654E-01
1.709460E-01 -3.993895E-01
1.650436E-01 -4.028388E-01
1.592249E-01 -4.063891E-01
1.534924E-01 -4.100543E-01
1.478452E-01 -4.139056E-01
1.422847E-01 -4.180117E-01
1.368152E-01 -4.223345E-01
1.314355E-01 -4.267376E-01
1.261475E-01 -4.311459E-01
1.209532E-01 -4.356480E-01
1.158522E-01 -4.403632E-01
1.108483E-01 -4.453171E-01
1.059408E-01 -4.504761E-01
1.011309E-01 -4.558385E-01
9.642112E-02 -4.614492E-01
9.181198E-02 -4.674174E-01
8.730499E-02 -4.738667E-01
8.290015E-02 -4.808608E-01
7.859901E-02 -4.884946E-01
7.440311E-02 -4.970084E-01
7.031356E-02 -5.066688E-01
6.633122E-02 -5.178317E-01
6.245677E-02 -5.312447E-01
5.869107E-02 -5.479705E-01
5.503523E-02 -5.690736E-01
5.149079E-02 -5.956618E-01
4.805775E-02 -6.282455E-01
4.473724E-02 -6.668327E-01
4.153121E-02 -7.103826E-01
3.843855E-02 -7.566098E-01
3.546149E-02 -8.012901E-01
3.259935E-02 -8.406817E-01
2.985477E-02 -8.726853E-01
2.722621E-02 -8.973400E-01
2.471635E-02 -9.216326E-01
2.232448E-02 -9.534629E-01
2.005172E-02 -9.878360E-01
1.789806E-02 -1.001857E+00
1.586618E-02 -9.842816E-01
1.395383E-02 -9.403147E-01
1.216255E-02 -8.681248E-01
1.049305E-02 -7.676828E-01
8.945737E-03 -6.248134E-01
7.519499E-03 -4.372711E-01
6.218114E-03 -2.513611E-01
5.037802E-03 -9.226103E-02
3.982344E-03  4.467917E-02
3.049502E-03  1.598130E-01
2.241935E-03  2.661596E-01
1.556561E-03  3.790513E-01
9.964624E-04  4.853633E-01
5.600984E-04  5.711626E-01
2.490099E-04  6.376294E-01
6.207850E-05  6.886544E-01
0.000000E+00  7.159074E-01
0.000000E+00  7.260201E-01
0.000000E+00  7.355866E-01
6.207850E-05  7.551654E-01
2.490099E-04  7.762503E-01
5.600984E-04  7.867794E-01
9.964624E-04  7.852415E-01
1.556561E-03  7.700414E-01
2.241935E-03  7.429880E-01
3.049502E-03  7.076319E-01
3.982344E-03  6.599380E-01
5.037802E-03  5.939070E-01
6.218114E-03  5.063455E-01
7.519499E-03  3.922463E-01
8.945737E-03  2.673375E-01
1.049305E-02  1.667324E-01
1.216255E-02  9.136314E-02
1.395383E-02  2.945198E-02
1.586618E-02 -1.870855E-02
1.789806E-02 -5.302520E-02
2.005172E-02 -7.114954E-02
2.232448E-02 -7.965068E-02
2.471635E-02 -8.896133E-02
2.722621E-02 -1.003542E-01
2.985477E-02 -1.101026E-01
3.259935E-02 -1.150294E-01
3.546149E-02 -1.149428E-01
3.843855E-02 -1.104630E-01
4.153121E-02 -1.034008E-01
4.473724E-02 -9.585433E-02
4.805775E-02 -8.907414E-02
5.149079E-02 -8.400794E-02
5.503523E-02 -8.092010E-02
5.869107E-02 -7.984806E-02
6.245677E-02 -8.049764E-02
6.633122E-02 -8.238485E-02
7.031356E-02 -8.499460E-02
7.440311E-02 -8.799645E-02
7.859901E-02 -9.121581E-02
8.290015E-02 -9.450951E-02
8.730499E-02 -9.781138E-02
9.181198E-02 -1.011291E-01
9.642112E-02 -1.044523E-01
1.011309E-01 -1.077208E-01
1.059408E-01 -1.108833E-01
1.108483E-01 -1.139207E-01
1.158522E-01 -1.168534E-01
1.209532E-01 -1.197270E-01
1.261475E-01 -1.225282E-01
1.314355E-01 -1.251978E-01
1.368152E-01 -1.276881E-01
1.422847E-01 -1.300529E-01
1.478452E-01 -1.323889E-01
1.534924E-01 -1.347092E-01
1.592249E-01 -1.369700E-01
1.650436E-01 -1.391859E-01
1.709460E-01 -1.413813E-01
1.769294E-01 -1.434932E-01
1.829941E-01 -1.454576E-01
1.891378E-01 -1.473116E-01
1.953579E-01 -1.491486E-01
2.016550E-01 -1.510364E-01
2.080261E-01 -1.529860E-01
2.144701E-01 -1.549834E-01
2.209855E-01 -1.570330E-01
2.275702E-01 -1.591263E-01
2.342220E-01 -1.612451E-01
2.409412E-01 -1.633918E-01
2.477244E-01 -1.655706E-01
2.545705E-01 -1.677661E-01
2.614786E-01 -1.699838E-01
2.684457E-01 -1.722292E-01
2.754698E-01 -1.744859E-01
2.825499E-01 -1.767507E-01
2.896854E-01 -1.790291E-01
2.968723E-01 -1.813000E-01
3.041100E-01 -1.835606E-01
3.113959E-01 -1.858010E-01
3.187296E-01 -1.879775E-01
3.261079E-01 -1.901032E-01
3.335299E-01 -1.921423E-01
3.409935E-01 -1.939915E-01
3.484972E-01 -1.957066E-01
3.560382E-01 -1.973219E-01
3.636152E-01 -1.987496E-01
3.712256E-01 -2.000680E-01
3.788689E-01 -2.013278E-01
3.865414E-01 -2.024772E-01
3.942432E-01 -2.035361E-01
4.019712E-01 -2.043627E-01
4.097227E-01 -2.048496E-01
4.174974E-01 -2.051267E-01
4.252926E-01 -2.051429E-01
4.331069E-01 -2.048765E-01
4.409375E-01 -2.044293E-01
4.487819E-01 -2.037550E-01
4.566403E-01 -2.029046E-01
4.645090E-01 -2.018148E-01
4.723873E-01 -2.003118E-01
4.802714E-01 -1.985454E-01
4.881610E-01 -1.965591E-01
4.960533E-01 -1.942904E-01
5.039466E-01 -1.918208E-01
5.118389E-01 -1.891170E-01
5.197285E-01 -1.861965E-01
5.276126E-01 -1.831022E-01
5.354909E-01 -1.797712E-01
5.433585E-01 -1.762765E-01
5.512169E-01 -1.726369E-01
5.590624E-01 -1.688038E-01
5.668930E-01 -1.648517E-01
5.747062E-01 -1.607936E-01
5.825025E-01 -1.566168E-01
5.902772E-01 -1.523688E-01
5.980287E-01 -1.480470E-01
6.057567E-01 -1.436623E-01
6.134574E-01 -1.392527E-01
6.211311E-01 -1.348244E-01
6.287743E-01 -1.303890E-01
6.363847E-01 -1.259647E-01
6.439617E-01 -1.215578E-01
6.515023E-01 -1.171734E-01
6.590053E-01 -1.128266E-01
6.664689E-01 -1.085289E-01
6.738909E-01 -1.042794E-01
6.812692E-01 -1.000901E-01
6.886029E-01 -9.596319E-02
6.958899E-01 -9.188418E-02
7.031276E-01 -8.785770E-02
7.103145E-01 -8.389001E-02
7.174489E-01 -7.999125E-02
7.245301E-01 -7.618695E-02
7.315542E-01 -7.247485E-02
7.385213E-01 -6.882942E-02
7.454283E-01 -6.524009E-02
7.522743E-01 -6.169981E-02
7.590587E-01 -5.819547E-02
7.657768E-01 -5.470983E-02
7.724298E-01 -5.123331E-02
7.790144E-01 -4.777299E-02
7.855298E-01 -4.431040E-02
7.919738E-01 -4.081589E-02
7.983449E-01 -3.729454E-02
8.046409E-01 -3.374165E-02
8.108621E-01 -3.011715E-02
8.170058E-01 -2.642009E-02
8.230694E-01 -2.267027E-02
8.290539E-01 -1.882740E-02
8.349563E-01 -1.486141E-02
8.407739E-01 -1.081765E-02
8.465075E-01 -6.733366E-03
8.521547E-01 -2.591115E-03
8.577141E-01  1.615146E-03
8.631847E-01  5.842932E-03
8.685644E-01  1.005871E-02
8.738524E-01  1.430265E-02
8.790467E-01  1.862645E-02
8.841466E-01  2.310251E-02
8.891516E-01  2.782963E-02
8.940591E-01  3.289906E-02
8.988679E-01  3.847440E-02
9.035777E-01  4.473683E-02
9.081879E-01  5.183891E-02
9.126949E-01  5.975033E-02
9.170998E-01  6.818987E-02
9.214009E-01  7.680731E-02
9.255957E-01  8.531936E-02
9.296863E-01  9.359404E-02
9.336687E-01  1.016416E-01
9.375431E-01  1.095757E-01
9.413088E-01  1.174864E-01
9.449647E-01  1.253432E-01
9.485080E-01  1.331216E-01
9.519410E-01  1.408655E-01
9.552616E-01  1.485496E-01
9.584687E-01  1.562093E-01
9.615614E-01  1.639854E-01
9.645384E-01  1.718557E-01
9.674006E-01  1.797395E-01
9.701451E-01  1.879138E-01
9.727726E-01  1.967841E-01
9.752836E-01  2.062378E-01
9.776754E-01  2.159892E-01
9.799482E-01  2.263733E-01
9.821007E-01  2.378616E-01
9.841337E-01  2.502562E-01
9.860461E-01  2.628762E-01
9.878374E-01  2.751076E-01
9.895069E-01  2.865729E-01
9.910542E-01  2.967041E-01
9.924793E-01  3.049949E-01
9.937818E-01  3.114169E-01
9.949610E-01  3.161995E-01
9.960176E-01  3.199159E-01
9.969493E-01  3.229918E-01
9.977580E-01  3.255720E-01
9.984422E-01  3.278329E-01
9.990034E-01  3.296118E-01
9.994387E-01  3.309222E-01
9.997509E-01  3.319675E-01
9.999367E-01  3.318149E-01
9.999999E-01  3.314395E-01];


%% Scan

% Surface
fid = fopen(file1); % open file
header = fgetl(fid);
surfData = fscanf(fid,'%g %g %g %g', [4 inf]); % scan data
surfData = surfData'; % transpose to correct format
fid =fclose(fid); % close file

% Field
fid = fopen(file2); % open file
header = fgetl(fid);
fieldData = fscanf(fid,'%g %g %g %g %g %g %g %g %g %g %g', [8 inf]); % scan data
fieldData = fieldData'; % transpose to correct format
fid =fclose(fid); % close file

%% Post-processing

% Surface
surfData_X = reshape(surfData(:,1),n,m);
surfData_Y = reshape(surfData(:,2),n,m);
surfData_Z = reshape(surfData(:,3),n,m);
surfData_Cp = reshape(surfData(:,4),n,m);

% Field
xMax=max(fieldData(:,1));
yMax=max(fieldData(:,2));
zMax=max(fieldData(:,3));
xMin=min(fieldData(:,1));
yMin=min(fieldData(:,2));
zMin=min(fieldData(:,3));
[fieldData_X,fieldData_Y,fieldData_Z] = meshgrid(linspace(xMin,xMax,x),linspace(yMin,yMax,y),linspace(zMin,zMax,z));
fieldData_ux = griddata(fieldData(:,1),fieldData(:,2),fieldData(:,3),fieldData(:,4),fieldData_X,fieldData_Y,fieldData_Z);
fieldData_uy = griddata(fieldData(:,1),fieldData(:,2),fieldData(:,3),fieldData(:,5),fieldData_X,fieldData_Y,fieldData_Z);
fieldData_uz = griddata(fieldData(:,1),fieldData(:,2),fieldData(:,3),fieldData(:,6),fieldData_X,fieldData_Y,fieldData_Z);
fieldData_mach = griddata(fieldData(:,1),fieldData(:,2),fieldData(:,3),fieldData(:,7),fieldData_X,fieldData_Y,fieldData_Z);
fieldData_sigma = griddata(fieldData(:,1),fieldData(:,2),fieldData(:,3),fieldData(:,8),fieldData_X,fieldData_Y,fieldData_Z);

% Min Velocity
minU = min(fieldData_ux(fieldData_ux > 0.01));
maxU = max(fieldData_ux(fieldData_ux > 0.01));
minW = min(fieldData_uz(fieldData_ux > 0.01));
maxW = max(fieldData_uz(fieldData_ux > 0.01));

% Min Mach number
minMach = min(fieldData_mach(fieldData_mach > 0.01));
maxMach = max(fieldData_mach(fieldData_mach > 0.01));

%% Display

% 3D mesh
figure
hold on
surf(surfData_X,surfData_Y,surfData_Z,ones(n,m),'EdgeColor','k');
for i=1:y
    surf(reshape(fieldData_X(i,:,:),x,z),reshape(fieldData_Y(i,:,:),x,z),reshape(fieldData_Z(i,:,:),x,z),ones(x,z),'FaceColor','none','EdgeColor','k');
end
colormap([200,200,200]./255);
axis([boundX(1),boundX(2),boundY(1),boundY(2),boundZ(1),boundZ(2)]);
daspect([1 1 1])
view(3);
camlight;
lighting gouraud;
xlabel('$x$', 'Fontsize', 16, 'Interpreter', 'Latex');
ylabel('$y$', 'Fontsize', 16, 'Interpreter', 'Latex');
zlabel('$z$', 'Fontsize', 16, 'Interpreter', 'Latex');
title('Mesh', 'Fontsize', 20, 'Interpreter', 'Latex');

% 3D surface pressure
figure
hold on
surf(surfData_X,surfData_Y,surfData_Z,surfData_Cp,'EdgeColor','none');
daspect([1 1 1])
view(3);
xlabel('$x$', 'Fontsize', 16, 'Interpreter', 'Latex');
ylabel('$y$', 'Fontsize', 16, 'Interpreter', 'Latex');
zlabel('$z$', 'Fontsize', 16, 'Interpreter', 'Latex');
title('Pressure distribution', 'Fontsize', 20, 'Interpreter', 'Latex');

% 2D plot velocity
figure;
hold on
plot(surfData(1:n,1),surfData(1:n,3),'Color','k','LineWidth',2);
quiver(reshape(fieldData_X(1,:,:),x,z),reshape(fieldData_Z(1,:,:),x,z),reshape(fieldData_ux(1,:,:),x,z),reshape(fieldData_uz(1,:,:),x,z));
xlabel('$x$', 'Fontsize', 16, 'Interpreter', 'Latex');
ylabel('$z$', 'Fontsize', 16, 'Interpreter', 'Latex');
title('Velocity', 'Fontsize', 20, 'Interpreter', 'Latex');

% 2D countour slice (x-velocity)
figure
hold on
plot(surfData(1:n,1),surfData(1:n,3),'Color','k','LineWidth',2);
contour(reshape(fieldData_X(1,:,:),x,z),reshape(fieldData_Z(1,:,:),x,z),reshape(fieldData_ux(1,:,:),x,z),'LevelList', linspace(minU,maxU,10));
colorbar;
caxis([minU maxU]);
xlabel('$x$', 'Fontsize', 16, 'Interpreter', 'Latex');
ylabel('$z$', 'Fontsize', 16, 'Interpreter', 'Latex');
title('X-velocity contour at station 1', 'Fontsize', 20, 'Interpreter', 'Latex');

% 2D countour slice (z-velocity)
figure
hold on
plot(surfData(1:n,1),surfData(1:n,3),'Color','k','LineWidth',2);
contour(reshape(fieldData_X(1,:,:),x,z),reshape(fieldData_Z(1,:,:),x,z),reshape(fieldData_uz(1,:,:),x,z));
colorbar;
caxis([minW maxW]);
xlabel('$x$', 'Fontsize', 16, 'Interpreter', 'Latex');
ylabel('$z$', 'Fontsize', 16, 'Interpreter', 'Latex');
title('Z-velocity contour at station 1', 'Fontsize', 20, 'Interpreter', 'Latex');

% 3D slice (mach)
figure
hold on
surf(reshape(fieldData_X(1,:,:),x,z),reshape(fieldData_Y(1,:,:),x,z),reshape(fieldData_Z(1,:,:),x,z),reshape(fieldData_mach(1,:,:),x,z),'EdgeColor','none');
surf(reshape(fieldData_X(floor(y/2),:,:),x,z),reshape(fieldData_Y(floor(y/2),:,:),x,z),reshape(fieldData_Z(floor(y/2),:,:),x,z),reshape(fieldData_mach(floor(y/2),:,:),x,z),'EdgeColor','none');
surf(reshape(fieldData_X(y-2,:,:),x,z),reshape(fieldData_Y(y-2,:,:),x,z),reshape(fieldData_Z(y-2,:,:),x,z),reshape(fieldData_mach(y-2,:,:),x,z),'EdgeColor','none');
colorbar;
caxis([minMach maxMach]);
axis([boundX(1),boundX(2),boundY(1),boundY(2),boundZ(1),boundZ(2)]);
daspect([1 1 1])
view(3);
xlabel('$x$', 'Fontsize', 16, 'Interpreter', 'Latex');
ylabel('$y$', 'Fontsize', 16, 'Interpreter', 'Latex');
zlabel('$z$', 'Fontsize', 16, 'Interpreter', 'Latex');
title('Mach number', 'Fontsize', 20, 'Interpreter', 'Latex');

% 3D slice (sigma)
figure
hold on
surf(reshape(fieldData_X(1,:,:),x,z),reshape(fieldData_Y(1,:,:),x,z),reshape(fieldData_Z(1,:,:),x,z),reshape(fieldData_sigma(1,:,:),x,z),'EdgeColor','none');
surf(reshape(fieldData_X(floor(y/2),:,:),x,z),reshape(fieldData_Y(floor(y/2),:,:),x,z),reshape(fieldData_Z(floor(y/2),:,:),x,z),reshape(fieldData_sigma(floor(y/2),:,:),x,z),'EdgeColor','none');
surf(reshape(fieldData_X(y-2,:,:),x,z),reshape(fieldData_Y(y-2,:,:),x,z),reshape(fieldData_Z(y-2,:,:),x,z),reshape(fieldData_sigma(y-2,:,:),x,z),'EdgeColor','none');
colorbar;
axis([boundX(1),boundX(2),boundY(1),boundY(2),boundZ(1),boundZ(2)]);
daspect([1 1 1])
view(3);
xlabel('$x$', 'Fontsize', 16, 'Interpreter', 'Latex');
ylabel('$y$', 'Fontsize', 16, 'Interpreter', 'Latex');
zlabel('$z$', 'Fontsize', 16, 'Interpreter', 'Latex');
title('Field source', 'Fontsize', 20, 'Interpreter', 'Latex');

% 2D countour slice (mach)
figure
hold on
plot(surfData(1:n,1),surfData(1:n,3),'Color','k','LineWidth',2);
contour(reshape(fieldData_X(1,:,:),x,z),reshape(fieldData_Z(1,:,:),x,z),reshape(fieldData_mach(1,:,:),x,z), 'LevelList', linspace(minMach,maxMach,10));
colorbar;
caxis([minMach maxMach]);
xlabel('$x$', 'Fontsize', 16, 'Interpreter', 'Latex');
ylabel('$z$', 'Fontsize', 16, 'Interpreter', 'Latex');
title('Mach countour at station 1', 'Fontsize', 20, 'Interpreter', 'Latex');

% 2D countour slice (sigma)
figure
hold on
plot(surfData(1:n,1),surfData(1:n,3),'Color','k','LineWidth',2);
contour(reshape(fieldData_X(1,:,:),x,z),reshape(fieldData_Z(1,:,:),x,z),reshape(fieldData_sigma(1,:,:),x,z));
colorbar;
xlabel('$x$', 'Fontsize', 16, 'Interpreter', 'Latex');
ylabel('$z$', 'Fontsize', 16, 'Interpreter', 'Latex');
title('Field source contour at station 1', 'Fontsize', 20, 'Interpreter', 'Latex');

% Pressure at MAC
figure
hold on
set(gca,'YDir','reverse');
set(gca,'Fontsize', 17,'Fontname', 'Times', 'LineWidth', 0.5);
plot(tranData(:,1),tranData(:,2),'Color',lightBlue,'LineWidth',2);
plot((surfData_X(:,5)-min(surfData_X(:,5)))/(max(surfData_X(:,5))-min(surfData_X(:,5))),surfData_Cp(:,5),'Color',darkRed,'LineWidth',2);
plot([0 1],[Cp_star Cp_star],'Color','k','LineWidth',1.2,'LineStyle','--');
leg=legend('Tranair', 'Code', '$C_{p}^{\star}$'); legend BOXOFF;
leg.FontSize = 16;
set(leg,'Interpreter','Latex')
xlabel('$x$', 'Interpreter', 'Latex');
ylabel('$C_p$', 'Interpreter', 'Latex');
title('$C_p$ at MAC', 'Interpreter', 'Latex');