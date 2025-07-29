function [  ] = showDEMProxy( dem )

demProxy_i_j = proxyDEM(dem);
demProxy_j_i = demProxy_i_j';
demProxy_jInv_i = demProxy_j_i(end:-1:1,:);

demProxyPlot = demProxy_jInv_i;

isValue = demProxyPlot > -9999;

maxDem = max(demProxyPlot(isValue));
minDem = min(demProxyPlot(isValue));

epsilon = 1e-6;

% normalized = ((demProxy_jInv_i) - minDem) / (maxDem - minDem);
logi = log10(abs(demProxyPlot)+epsilon);

maxi = max(logi(isValue));
mini = min(logi(isValue));

noValue = -0.1;

normalized = mini + noValue + 0*demProxyPlot;

normalized(isValue) = (logi(isValue) - mini) / (maxi - mini);

imshow(normalized,[noValue, 1]);

end

