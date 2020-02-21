function uaFN = GetSafeAlphaController(polyMod,Tmax,vThresh)
nModes = length(polyMod.x);

thetaLim = linspace(-1,1,nModes+1);
thetaLim(1) = -inf;
thetaLim(end) = inf;

uHat_FN = maskInputInterpHybrid(polyMod.invSetOut.u, polyMod.invSetOut.v, polyMod.x, vThresh);

function ua = UnscaledMaskedController(xHat,ua0)
    xHatSC = polyMod.scaleParam.x_SC_Fun{1}(xHat);
    thetaSC = xHatSC(1);
    m = find(thetaSC < thetaLim,1)-1;
    
    uHat0 = UalphaToUhat(ua0,xHat,Tmax);
    
    uHat = uHat_FN{m}(xHatSC,uHat0);
    
    ua = UhatToUalpha(uHat,xHat,Tmax);
end

uaFN = @UnscaledMaskedController;

end