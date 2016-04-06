function [root_mult_array] = o_roots_BezierClipping(fx)

% Pseudocode from 'Computing roots of polynomials by quadratic clipping'
% Barton and Juttler
%
% Algorithm bezclip(p,[a,b])
% where p polynomial
% 1. if length of interval [a,b] > e then
% 2.    C <- Convex hull of control points of p with respect to [a,b]
% 3.    if C intersects t axis then
% 4.        Find [a',b'] by intersecting C with the t axis
% 5.        if |a'-b'|< 1/2*|a-b| then
% 6.            return (bezclip(p,[a',b'])
% 7.        else
% 8.            return (bezclip(p,[a,1/2(a+b)]) union bezclip(p,1/2(a+b,b))
% 9.        end if
% 10.   else
% 11.   return empty set
% 12.   end if
% 13. else
% 14. return ([a,b])
% 15. end if


% bool_preproc = '';
% low_rank_approx_method = '';
% apf_method = '';
% SetGlobalVariables(bool_preproc,low_rank_approx_method,apf_method);

% Set interval size
global xmin
global xmax
global ymin
global ymax


% max error
epsilon = 1e-5;

% Set iteration number to one
outer_loop_ite = 1;

% Initialise a vector to store calculated roots
root_mult_array = [];

% While f(x) is not a constant
while GetDegree(fx) >= 1
    
    if GetDegree(fx) == 1;
        
        
        % Curve is a straight line.
        m = 1;
        x = linspace(0,1,m+1)';
        CP = [x,fx];
        
        xmin = min(CP(:,1));
        xmax = max(CP(:,1));
        ymin = min(CP(:,2));
        ymax = max(CP(:,2));
        
        x_intercept_new = GetXIntercept(CP);
        fprintf('x intercept : %2.4f',x_intercept_new)
    else
        
        % Get degree of polynomial f(x)
        m = GetDegree(fx);
        
        % Define the set of control points of f(x)
        % Split the unit interval into m+1 parts
        x = linspace(0,1,m+1)';
        CP = [x,fx];
        
        xmin = min(CP(:,1));
        xmax = max(CP(:,1));
        ymin = min(CP(:,2));
        ymax = max(CP(:,2));
        
        x_intercept_old = 0;
        x_intercept_new = 1;
        inner_loop_ite = 1;
        
        while abs(x_intercept_new - x_intercept_old) > epsilon
            
            % Get the convex hull of the control points of f(x)
            k = convhull(CP(:,1),CP(:,2));
            
            convex_hull_vertices = flipud(CP(k,:));
            
            % for each line of the convex hull, check if it intersects
            % Plot the control points
            Plot_fx(CP,k)
            
            % Get the first point at which the convex hull crosses the x axis.
            x_intercept_old = x_intercept_new;
            x_intercept_new = GetXIntercept(convex_hull_vertices);
           
            if(x_intercept_new == -1)
                break;
            end
            
            % Evaluate f(x) at the point.
            Bernstein_Evaluate(fx,x_intercept_new);
            
            % Subdivide at the point t1.
            [~, Pk_right] = BezierSubdivide(CP,m,x_intercept_new);
            
            CP = Pk_right;
            inner_loop_ite = inner_loop_ite + 1;
        end
    end
    % Save all plots
    
    if (x_intercept_new == -1)
        break;
    end
    
    %dir_name = sprintf('outputs-%i',outer_loop_ite);
    %save_all_figures_to_directory(dir_name);
    
    % Close all open plots
    close all;
    
    
    % Get the root.
    root = x_intercept_new;
    
    % Add the root to the root array.
    root_mult_array = [root_mult_array ; root 1 ];
    
    % Get polynomial of the factor (x-r)
    tx =[...
        -root;
        1-root;
        ];
    
    % Deconvolve root from polynomial
    fx = Bernstein_Deconvolve(fx,tx);
    
    % Increment iteration number
    outer_loop_ite = outer_loop_ite + 1;
    
end


% % Print out roots
PrintoutRoots('CLIPPING',root_mult_array);


end


function [new_left_pk, new_right_pk] = BezierSubdivide(Pk,degree,c)
% Given a set of control points, subdivide such that two new sets of
% control points are obtained.
%
% Inputs.
%
%
% Pk :  Set of original control points.
%
% degree :  degree of control points
%
% c :    point of subdivision.
%
%                           Outputs.
%
%
% new_left_pk : - The set of control points on the left of c.
%
% new_right_pk : - The set of control points on the right of c.
%

% obtain each set of control points
Pk_array{1} = Pk;

for i = 2:1:degree+1
    Pk_array{i} = deCasteljau(c,Pk_array{i-1})  ;
end

% From the generated sets of control points obtain set of control points
% for sub-curve P[a,c]
new_left_pk = zeros(degree+1,2);
for i = 1:1:degree+1
    new_left_pk(i,:) = Pk_array{i}(1,:);
end

% obtain set of control points for sub-curve P[c,b]
new_right_pk = zeros(degree+1,2);
for i=1:1:degree+1
    new_right_pk(i,:) = Pk_array{i}(end,:);
end

new_right_pk(i,:) = Pk_array{i}(end,:);

new_right_pk = sortrows(new_right_pk,1);

end

function [Pk_1] = deCasteljau(c, Pk_0)
% Obtain the next set of control points, given a set of control points Pk_0
% in the deCasteljau subdivision process.
% Uses deCasteljau algorithm found in 'Applied Geometry for computer
% graphics and CAD' - Page 151

% c = subdivision point on horizontal axis
% Pk_0 = original control points

% get first control poitn
a = Pk_0(1,1);

% get last control point
b = Pk_0(end,1);


T = (c-a)./(b-a);


% Build set of control points
Pk_1 = [];
% for
for i = 1:1:length(Pk_0)-1
    %Pk_1(i)
    aa =  (1-T).*Pk_0(i,:) + T.*Pk_0(i+1,:) ;
    Pk_1 = [Pk_1 ; aa];
    
end



end


function [x_intercept] = GetXIntercept(convex_hull_vertices)
% Get the intersections between the convex hull and the x axis, return the
% first intercept.
%
% Inputs.
%
%
% convex_hull_vertices : Set of [x,y] paris of vertices of the convex hull.
%                        where the first vertex is included twice.
%
% Outputs.
%
% x_intercept : the first point at which the convex hull crosses the x
%               axis.
%



% Get the number of vertices in the convex hull
[nVerticesCH,~] = size(convex_hull_vertices);

% Note that convex_hull_vertices contains the first vertex twice
nVerticesCH = nVerticesCH - 1;

% Get the set of edges which make up the convex hull
nEdges = nVerticesCH;

% Initialise a vector to store xaxis intercepts
v_x_intercept = [];

for i =1:1:nEdges
    x0 = convex_hull_vertices(i,1);
    y0 = convex_hull_vertices(i,2);
    x1 = convex_hull_vertices(i+1,1);
    y1 = convex_hull_vertices(i+1,2);
    
    if (y0*y1) > 0 % No change of sign
        % Do Nothing
    else % Change of sign
        
        
        % Get x intercept
        m = (y1-y0)./ (x1-x0);
        
        % Get location of intercept on x axis
        location = x1 - (y1./m);
        
        v_x_intercept = ...
            [
            v_x_intercept ;
            location
            ];
        
        
    end
end

if size(v_x_intercept,1)<1
    x_intercept = -1;
    return;
end

% Get all x intercepts in order
v_x_intercept = sort(v_x_intercept,'ascend');

% Get the first intercept
x_intercept = v_x_intercept(1);
end


function save_all_figures_to_directory(FolderName)

% First, get the name of the folder you're using.
% For example if your folder is 'D:\photos\Info',
% parentFolder  would = 'D:\photos, and deepestFolder
% would = 'Info'.

myFolder = pwd;

[parentFolder,deepestFolder] = fileparts(myFolder);

% Next, create a name for a subfolder within that.
% For example 'D:\photos\Info\DATA-Info'
newSubFolder = sprintf('%s/Outputs/DATA-%s', myFolder, FolderName);

% Finally, create the folder if it doesn't exist already.
if ~exist(newSubFolder, 'dir')
    mkdir(newSubFolder);
end



figlist=findobj('type','figure');

% for each figure in the figure list
for i=1:numel(figlist)
    file_name = sprintf('Figure-%i',i);
    ext = sprintf('.fig');
    ext_jpeg = sprintf('.jpg');
    saveas(figlist(i),fullfile(newSubFolder,[file_name,ext]));
    saveas(figlist(i),fullfile(newSubFolder,[file_name,ext_jpeg]));
end
end