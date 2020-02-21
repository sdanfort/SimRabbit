function right_impact_constraints(nlp, src, tar, bounds, varargin)
    
    plant = nlp.Plant;
    
    % fist call the class method
    plant.rigidImpactConstraint(nlp, src, tar, bounds, varargin{:});
    
    % Don't need time continuity constraint
    removeConstraint(nlp,'tContDomain');

    % the relabeling of joint coordiante is no longer valid
    removeConstraint(nlp,'xDiscreteMapRightImpact');
    
    % the configuration only depends on the relabeling matrix
    R = plant.R;
    x = plant.States.x;
    xn = plant.States.xn;
    x_diff = R*x-xn;
%     footHeight = 1000*(plant.ImpactConstraints.RightToe.ConstrExpr(3) + plant.ImpactConstraints.RightToe.Param(3));
    x_map = SymFunction(['xDiscreteMap' plant.Name],x_diff(2:end),{x,xn});
    sl_fun = SymFunction(['stepLength' plant.Name],x_diff(1),{x,xn});
%     swfh_fun = SymFunction(['swingFootHeight' plant.Name],footHeight,{x});
%     stfh_fun = SymFunction(['stanceFootHeight' plant.Name],subs(footHeight,x,xn),{xn});
    
    addNodeConstraint(nlp, x_map, {'x','xn'}, 'first', 0, 0, 'Linear');
    addNodeConstraint(nlp, sl_fun, {'x','xn'}, 'first', 0.5, 0.5, 'Linear');
%     addNodeConstraint(nlp, swfh_fun, {'x'}, 'first', 0, 0);
%     addNodeConstraint(nlp, stfh_fun, {'xn'}, 'first', 0, 0);
end