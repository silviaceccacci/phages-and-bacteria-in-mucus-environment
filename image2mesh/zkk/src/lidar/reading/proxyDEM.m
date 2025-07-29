function [ demProxy ] = proxyDEM(dem)

demProxy = dem.elevationProxy_I_J;

demProxy = cell2mat(demProxy);

end

