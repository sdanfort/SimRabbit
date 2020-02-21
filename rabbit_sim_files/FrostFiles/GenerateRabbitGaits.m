function gait = GenerateRabbitGaits(stepLengths,avgVels,wRefVal,COMPILE)

%% Setup
cur = pwd;
addpath(genpath(cur));

% addpath('../../');
% frost_addpath;
export_path = fullfile(cur,'gen/');
% if load_path is empty, it will not load any expression.
% if non-empty and assigned valid directory, then symbolic expressions will
% be loaded  from the MX binary files from the given directory.
load_path = [];%fullfile(cur, 'gen/sym');
delay_set = false;
if nargin < 2
    COMPILE = true;
end

% Load model
rabbit = RABBIT('urdf/five_link_walker.urdf');
if isempty(load_path)
    rabbit.configureDynamics('DelayCoriolisSet',delay_set);
else
    % load symbolic expression for the dynamics equations
    rabbit.loadDynamics(load_path, delay_set);
end


% Define domains
r_stance = RightStance(rabbit, load_path);
% l_stance = LeftStance(rabbit, load_path);
r_impact = RightImpact(r_stance, load_path);
% l_impact = LeftImpact(l_stance, load_path);

% Define hybrid system
rabbit_1step = HybridSystem('Rabbit_1step');
rabbit_1step = addVertex(rabbit_1step, 'RightStance', 'Domain', r_stance);

srcs = {'RightStance'};
tars = {'RightStance'};

rabbit_1step = addEdge(rabbit_1step, srcs, tars);
rabbit_1step = setEdgeProperties(rabbit_1step, srcs, tars, ...
    'Guard', {r_impact});

%% Define User Constraints
r_stance.UserNlpConstraint = str2func('right_stance_constraints');
r_impact.UserNlpConstraint = str2func('right_impact_constraints');

%% Define User Costs
u = r_stance.Inputs.Control.u;
ucost = tovector(sum(u.^2));
ucost_fun = SymFunction(['ucost_' r_stance.Name],ucost,{u});

%% Create optimization problem
num_grid.RightStance = 10;
num_grid.LeftStance = 10;
nlp = HybridTrajectoryOptimization('Rabbit_1step',rabbit_1step,num_grid,...
    [],'EqualityConstraintBoundary',1e-6);

% Configure bounds 
bounds = setBounds(rabbit);

% load some optimization related expressions here
if ~isempty(load_path)
    nlp.configure(bounds, 'LoadPath',load_path);
else
    nlp.configure(bounds);
end

nlp.Phase(1).OptVarTable.wref(1).setBoundary(0,0);
nlp.Phase(1).OptVarTable.xref(1).setBoundary(0,0);
 
% Add costs
addRunningCost(nlp.Phase(getPhaseIndex(nlp,'RightStance')),ucost_fun,{'u'});

% Update
nlp.update;

x = r_stance.States.x;
xref = r_stance.Params.xref;
wref = r_stance.Params.wref;
for i = 1:nlp.Phase(1).NumNode
    xcost_ = tovector(sum((x-xref(7*(i-1) + (1:7))).^2)*wref);
    xcost_fun_ = SymFunction(['xcost' num2str(i) '_' r_stance.Name],xcost_,{x,xref,wref});
    
    addNodeCost(nlp.Phase(1),xcost_fun_,{'x','xref','wref'},i);
end

nlp.update;

% save expressions after you run the optimization. It will save all required
% expressions
% do not need to save expressions if the model configuration is not
% changed. Adding custom constraints does not require saving any
% expressions.
% load_path = fullfile(cur, 'gen/sym');
% rabbit_1step.saveExpression(load_path);
%% Compile
if COMPILE
    if ~exist([export_path, 'opt/'])
        mkdir([export_path, 'opt/'])
    end
    rabbit.ExportKinematics([export_path,'kinematics/']);
    compileConstraint(nlp,[],[],[export_path, 'opt/']);
    compileObjective(nlp,[],[],[export_path, 'opt/']);
end

% Example constraint removal
% removeConstraint(nlp.Phase(1),'u_friction_cone_RightToe');

%% Create Ipopt solver for each step length and solve
addpath(genpath(export_path));
options.tol = 1e-6;
options.acceptable_tol = 1e-6;
options.constr_viol_tol = 1e-6;
options.acceptable_constr_viol_tol = 1e-6;

gait = cell(0,3);

for i = 1:length(stepLengths)
    updateConstrProp(nlp.Phase(2),'stepLengthRightImpact','first','lb',...
        stepLengths(i),'ub',stepLengths(i));
    updateConstrProp(nlp.Phase(1),'average_velocity','last','lb',...
        avgVels(i),'ub',avgVels(i));
    
    nlp.update;
    solver = IpoptApplication(nlp,options);

    % Run Optimization
    tic
    % old = load('x0');
    % [sol, info] = optimize(solver, old.sol);
    [sol, info] = optimize(solver);
    toc
    if ~info.status
        [tspan, states, inputs, params] = exportSolution(nlp, sol);
        gait{end+1,1} = stepLengths(i);
        %seems like here is where the alpha values are generated:
        gait{end,2} = avgVels(i);
        gait{end,3} = struct(...
                        'tspan',tspan,...
                        'states',states,...
                        'inputs',inputs,...
                        'params',params);
                
        % Animate
        q_log_R = states{1}.x; % Right stance
        
        q_log_L = q_log_R([1:3,6:7,4:5],:); % symmetric Left stance
        
        q_log_L(1:3,:) = q_log_L(1:3,:) + repmat((q_log_R(1:3,end)-q_log_R(1:3,1)),1,21);

        t_log_R = tspan{1};
        t_log_L = t_log_R + t_log_R(end);

        q_log = [q_log_R, q_log_L];
        t_log = [t_log_R, t_log_L];
        
%         pToe = 0*t_log;
%         for j = 1:length(pToe)
%             pToe_ = p_LeftToe(q_log(:,j));
%             pToe(j) = pToe_(3);
%         end
%         figure,plot(t_log,pToe)
%         
%         close all
%         anim = Animator.FiveLinkAnimator(t_log, q_log);
%         anim.pov = Animator.AnimatorPointOfView.West;
%         anim.isLooping = true;
%         anim.updateWorldPosition = false;
%         anim.endTime = 2;
%         
%         for j = 1:20
%             anim.Animate
%             pause(0.05)
%         end
        
        for j = 1:length(nlp.VariableArray)
            nlp.VariableArray(j).setInitialValue(sol(nlp.VariableArray(j).Indices));
        end
        
        xRefVal = zeros(7*21,1);
        for j = 1:nlp.Phase(1).NumNode
            xRefVal(7*(j-1) + (1:7)) = sol(nlp.Phase(1).OptVarTable.x(j).Indices);
        end
        
        nlp.Phase(1).OptVarTable.wref(1).setBoundary(wRefVal,wRefVal);
        nlp.Phase(1).OptVarTable.xref(1).setBoundary(xRefVal,xRefVal);
        
        for j = 1:nlp.Phase(1).NumNode
            nlp.Phase(1).OptVarTable.wref(j).setInitialValue(wRefVal);
            nlp.Phase(1).OptVarTable.xref(j).setInitialValue(xRefVal);
        end
    
        nlp.update;
    else
        warning('Optimization did not converge')
    end
end

%% Remove outliers

% gaitNoOut = gait(1,:);
% 
% for i = 2:length(gait)
%     if sum(sum((gait{i,2}(1).states.dx - gaitNoOut{end,2}(1).states.dx).^2)) < 0.1
%         gaitNoOut{end+1,:} = gait{i,:};
%     end
% end
% 
% gait = gaitNoOut;

%% Visualize
nNodes = length(gait{1,3}(1).tspan);

if size(gait,1) > 1
    for i = 1:7
        figure(i),clf
        for j = 1:size(gait,1)
            x_(j,:) = gait{j,1}*ones(1,nNodes);
            y_(j,:) = gait{j,3}(1).tspan;
            z_(j,:) = gait{j,3}(1).states.x(i,:);
        end
        surf(x_,y_,z_,'linestyle','none')
    end
end

end
