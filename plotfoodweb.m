function plotfoodweb(A,varargin)
%PLOTFOODWEB(A) plots a foodweb from adjacency matrix A
% A must be square and contain 0 and 1. By covention A(x,y)= 1 means 
% that species 'y' eats species 'x'. NaNs and negative values will be 
% treated as 0, while positive values will be treated as 1. 
%
% By default each species is represented by a colored square (see below 
% high quality plot) with coordinates X, Y, Z. X-Y coordinates are 
% calculated automatically as two circles (one for plants and another 
% for animals). Z coordinates, if not specified, are the respective 
% trophic level of each species calculated with function 'trophiclevel' 
% (included).
%
% Optional arguments can be abreviated to 2 or 3 letters.
% 
% PLOTFOODWEB(...,'quality','high') plots in high quality (much slower).
% Species will be represented by semi-transparent spheres connected with 
% semi-transparent tapered cilinders (wider in the prey size). Default is
% 'low'
% 
% PLOTFOODWEB(...,'stransparency',ST) will plot the spheres with transparency 
% 'ST' (from 0 to 1, 0 means invisible). Default 0.6;
% 
% PLOTFOODWEB(...,'size',S) uses vector 'S' as sizes for the squares or 
% spheres. The default is 50 for all. In case of high quality plot (see 
% above) the default is 2.
% 
% PLOTFOODWEB(...,'color',C) assigns the colors of the squares. C must be a 
% matrix of N by 3 containing RBG values (between 0 and 1). If N is 2, the
% colors will be interpolated between C(1,:) (for plants) 
% and C(2,:) for top predators. If C is not provided, green will be used
% for plants while animal's color will be interplated between pale and dark blue
% according to their trophic level. If N is larger than 2 it must be 
% equal to the number of species, specifying a color for each species.
% Cannibals (i.e. A(i,i)=1) will be plotted in red.
% 
% PLOTFOODWEB(...,'lcolor',L) where L is a vector of 3 values fom 0 to 1, 
% assigns color L to the lines that connect predators and prey. Default 
% is [.7 .7 .7].
%
% PLOTFOODWEB(...,'ltransparency',LT) will plot lines that connect species
% with transparency 'LT' (from 0 to 1, 0 means invisible). Default 0.5. It
% also applies to the cylinders used in high quality plot.
% 
% PLOTFOODWEB(...,'lmethod',method) uses 'method' for rendering. Default is 
% 'flat'. The only other accepted value is 'gouraud' (slower, higher quality).
% 
% PLOTFOODWEB(...,'tl',T) uses vector T as trophic level and Z coordinate.
% 
% PLOTFOODWEB(...,'numbers','on') plots the number of the species close to the
% square (only works in low quality). Default 'off'.
%
% PLOTFOODWEB(...,'names','on') plots the name of the species close to the
% square. Names must be provided in a cell array. Default 'off'. If set to
% 'on', then numbers will be set to 'off'.
%
% PLOTFOODWEB(...,'arrange',method) specifies the method to arrange species 
% in space.  The default is 'circle' in which species's are in circles in the 
% XY plane and trophic level is used as the Z coordinate. Other possibilities 
% are:
%	'flat': the food web will be ploted in 2D. Plants and herbivores will 
%	be placed centered and uniformly in the X axis, at Z coordinates 0 and 1
%	respectively. For the rest of species, the Z coordinate will be the
%	trophic level and the X coordinate will be random.
% 
%	'mds' in this case a multi-dimensional scaling is performed based on the 
%	similitude of predators and prey for each species and the resulting 
%	coordinates used to arrange the species in the XY plane. The vertical 
%	value (Z) is again trophic level (requires statistical toolbox). 
%	Species will cluster based on predators and prey shared.
% 
% Requires (included in the zip):
% functions: trophiclevel, fwcylinder, fwmds, fwrotation
%
% Author: Fancisco de Castro (2015)


%-- Some checking of input
if nargin < 1					error('Need at least 1 argument: square adjacency matrix'); end
if size(A,1) ~= size(A,2)	error('Adjacency matrix must be square'); end
if ~isnumeric(A)				error('Adjacency matrix must be numeric'); end
A(isnan(A))= 0;
A(A < 0)= 0;
A(A > 0)= 1;


%-- Defaults, etc.
ns= size(A,1);
nprey= sum(A,1);
plants=  find(nprey == 0);
animals= find(nprey > 0);
bsize= [];
linecolor= [.7 .7 .7];
linetrnsp= 0.5;
spheretrnsp= 0.6;
xc= zeros(ns,1);
yc= zeros(ns,1);
zc= [];
color= [];
quality= 'low';
lmethod= 'flat';
numbers= 'off';
names= {};
arrange= 'circle';


%-- Optional arguments
for j= 1:2:length(varargin)
	string= lower(varargin{j});
	switch string(1:min(3,length(string)))
		case 'tl'
			zc= varargin{j+1};
		case 'col'
			color= varargin{j+1};
		case 'lco'
			linecolor= varargin{j+1};
		case 'str'
			spheretrnsp= varargin{j+1};
		case 'ltr'
			linetrnsp= varargin{j+1};
		case 'siz'
			bsize= varargin{j+1};
		case 'qua'
			quality= varargin{j+1};
		case 'lme'
			lmethod= varargin{j+1};
		case 'num'
			numbers= varargin{j+1};
			if numel(numbers) ~= ns
				error('Size of numbers array does not match the number of species');
			end
		case 'nam'
			names= varargin{j+1};
			numbers= 'off';
			names= names(:);
			if numel(names) ~= ns
				error('Size of names array does not match the number of species');
			end
		case 'arr'
			arrange= varargin{j+1};
		otherwise
			error('Unknown parameter name');
	end
end


%-- No trophic level assigned
if isempty(zc) zc= trophiclevel(A); end


%-- Cannibals & weird cannibals
cann= find(diag(A));
weird= find(zc(cann)== 1);
zc(cann(weird))= zc(cann(weird))+1;
maxzc= max(zc);


%-- X-Y coordinates
nplants=  length(plants);
nanimals= length(animals);
herbs=  find(zc==2);
carns=  find(zc> 2);
nherbs= numel(herbs);
ncarns= numel(carns);

switch arrange(1:3)
	case 'cir'	%-- Default X-Y Coordinates plants (circle)
		alfa= linspace(0,6.1,nplants)';
		xc(plants)= 0.5*maxzc*cos(alfa);
		yc(plants)= 0.5*maxzc*sin(alfa);
		% Animals
		alfa= linspace(0,6.2,nanimals)';
		xc(animals)= 0.5*maxzc*cos(alfa);
		yc(animals)= 0.5*maxzc*sin(alfa);
	case 'mds'
		[xc,yc,zc]= fwmds(A);
	case 'fla'
		centr= round(max([nplants,nherbs])/2);
		minxpl=  centr-floor((nplants-1)/2);
		maxxpl=  centr+floor(nplants/2);
		xc(plants)= minxpl:maxxpl;
		zc(plants)= zeros(nplants,1);

		minxhe=  centr-floor((nherbs-1)/2);
		maxxhe=  centr+floor(nherbs/2);
		xc(herbs)= minxhe:maxxhe;
		
		minx= min(minxpl,minxhe);
		maxx= max(maxxpl,maxxhe);
		xc(carns)= minx+ (maxx-minx)* rand(ncarns,1);
		yc= zeros(ns,1);
	otherwise
		error('Unknown arrange value');
end



%-- Colors
switch size(color,1)
	case 0 % No color provided: gradient light blue to dark blue
		color(plants,:)= repmat([0 1 0],nplants,1); % Plants all green
		slope= -1/(maxzc-1);
		color(animals,1:3)= [zeros(nanimals,1), 1+slope*(zc(animals)-1), ones(nanimals,1)];

	case 2 % 2 colors provided: gradient between them
		color(plants,:)= repmat(color(1,1:3),nplants,1); % Plants all base color
		slope= [(color(2,1)-color(1,1))/maxzc; ...
				  (color(2,2)-color(1,2))/maxzc; ...
				  (color(2,3)-color(1,3))/maxzc];
		for j= 1:nanimals
			color(animals(j),1:3)= [color(1,1)+slope(1)*(zc(animals(j))-2), ...
										   color(1,2)+slope(2)*(zc(animals(j))-2), ...
										   color(1,3)+slope(3)*(zc(animals(j))-2)];
		end

	otherwise % All colors provided: check size	
		if size(color,1) ~= ns error('Wrong number of rows in color array'); end
end

color(color > 1)= 1;
color(color < 0)= 0;


%########################
%####### Plotting
%########################

hold on
if strcmp(quality(1:3),'low') % Low quality

	%-- Sizes (if not provided)
	if isempty(bsize) bsize= 50*ones(ns,1); end

	%-- Dots
	scatter3(xc(:),yc(:),zc(:),bsize(:),color(:,1:3),'s','filled','MarkerEdgeColor','k');

	%-- Cannibals: red dot
	if ~isempty(cann)
% 		if size(bsize,1) cannsize= bsize 
		cannsize= bsize(cann);
% 		end
		scatter3(xc(cann),yc(cann),zc(cann),cannsize,'r','s','filled','MarkerEdgeColor','k');
	end


	%-- Lines
	for p= 1:length(animals)
		x1= xc(animals(p)); y1= yc(animals(p)); z1= zc(animals(p));
		presas= find(A(:,animals(p))== 1);
		for j= 1:length(presas)
			x2= xc(presas(j)); y2= yc(presas(j)); z2= zc(presas(j));
			patchline([x1 x2],[y1 y2],[z1 z2],'edgecolor',linecolor,'linewidth',0.1,'edgealpha',linetrnsp);
		end
	end

else % High quality, then

	%-- Sizes (if not provided)
	if isempty(bsize)		bsize= max([xc;yc;zc])*0.02*ones(ns,1); end
	if size(bsize,1)==1	bsize= bsize*0.02*ones(ns,1); end


	%-- Spheres
	[sx,sy,sz]= sphere(30);
	for j= 1:length(xc)
		surf(sx*bsize(j)+xc(j), sy*bsize(j)+yc(j), sz*bsize(j)+zc(j),...
		'LineStyle','none','AmbientStrength',1,'FaceColor',color(j,:),...
		'SpecularStrength',1,'DiffuseStrength',1,'FaceAlpha',spheretrnsp);
	end
	if ~isempty(cann)
		for j= 1:length(cann)
		surf(sx*bsize(cann(j))+xc(cann(j)), sy*bsize(cann(j))+yc(cann(j)), sz*bsize(cann(j))+zc(cann(j)),...
		'LineStyle','none','AmbientStrength',1,'FaceColor',[1 0.4 0],...
		'SpecularStrength',1,'DiffuseStrength',1,'FaceAlpha',spheretrnsp);
		end
	end

	%-- Lines
	diam= max(bsize)/8;
	for p= 1:length(animals)
		x1= xc(animals(p)); y1= yc(animals(p)); z1= zc(animals(p));
		presas= find(A(:,animals(p))==1);
		for j= 1:length(presas)
			x2= xc(presas(j)); y2= yc(presas(j)); z2= zc(presas(j));
			fwcylinder([x1,y1,z1],[x2,y2,z2],diam,20,[.6 .6 .6],linetrnsp);
		end
	end

	light('Position',[0 0 1],'Style','infinit','Color',[1 1 1]);
	light('Position',[0 -1 0],'Style','infinit','Color',[1 1 1]);
	material shiny
	
	switch lmethod
	case 'flat'
		lighting flat
	case 'gouraud'
		lighting gouraud
	otherwise
		error('Wrong value for lmethod')
	end

	clear sx sy sz x1 y1 z1 x2 y2 z2
end


%-- Species numbers or names
offset= max([xc;yc;zc])/50;
if isempty(names)
	if strcmp(numbers(1:2),'on') 
		for j= 1:ns
			text(offset+xc(j),offset+yc(j),offset+zc(j),num2str(j),'FontSize',8,'color',[0 0 0])
		end
	end
else
	for j= 1:ns
		text(offset+xc(j),offset+yc(j),offset+zc(j),names(j,:),'FontSize',8,'color',[0 0 0])
	end
end

axis equal
axis tight
axis off
view(0,0);
hold off



