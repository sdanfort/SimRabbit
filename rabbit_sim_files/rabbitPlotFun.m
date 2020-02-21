function rabbitPlotFun(x,v,~,~,~,polyMod)
    nModes = length(x);
    
    figure(1), clf, hold on
    dalpha_ = 0;
    
    xSym = sym('x',[4,1]);
    
    thetaLim = linspace(-1,1,nModes+1);

    for i = 1:nModes
        [alpha_,th_,dth_] = ndgrid(linspace(-1,1,40),linspace(thetaLim(i),thetaLim(i+1),ceil(50/nModes)),linspace(-1,1,50));
        xSC_ = [th_(:),dth_(:),alpha_(:),0*alpha_(:)+dalpha_];
        vfn_ = fn(v{i},x{i});
        vfn_ = matlabFunction(vfn_(xSym),'vars',{xSym});
        v_ = vfn_(xSC_.');
    %     v_ = full(msubs(v_0{i},polyMod.x{i},[phi_(:),th_(:),dphi_+0*th_(:),dth_(:)].'));
        v_ = reshape(v_,size(th_));
        
        patch_ = patch(isosurface(th_,alpha_,dth_,v_,0));

        patch_.FaceColor = [55,200,100]/256;%/256;
        patch_.EdgeColor = 'none';

        isonormals(th_,alpha_,dth_,v_,patch_);
    end
    ax_ = gca;
    
   
    xlabel('$\theta$','interpreter','latex')
    ylabel('$\alpha$','interpreter','latex')
    zlabel('$\dot{\theta}$','interpreter','latex')

    view([34.1800, 32.5600])
    camlight HEADLIGHT; lighting phong
    drawnow
    
%     phiRange_ = ax_.YLim;
%     dphiRange_ = ax_.ZLim;
%     
%     pg = patch(pi/12*[ 1; 1; 1; 1; 1],...
%                 phiRange_([1,1,2,2,1]),...
%                 dphiRange_([1,2,2,1,1]),[84,52,83]/100);
%         
%     set(pg,'facealpha',0.5)
    box on
end