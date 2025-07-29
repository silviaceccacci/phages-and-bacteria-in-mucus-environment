

function []=writeNodeStatistics(nodalDistortion,nameFile,numDegIni,numDegGeo)


nameTextFile = ['./results/exportMeshes/' nameFile '.txt'];
fopen(nameTextFile);
fid = fopen(nameTextFile, 'wt');

fprintf(fid, 'Optimized quality:\n');
fprintf(fid, '\\begin{table}[H] \n');
fprintf(fid, '\\centering \n');
% fprintf(fid, '\\begin{tabular}{c|c|c|c|c|c|} \n');
fprintf(fid, '\\begin{tabular}{cccccc} \n');
fprintf(fid, '\\hline \n' );
fprintf(fid, '\\cline{2-6}      & Min. Q.   & Max. Q. & Mean Q. & Std. Dev. & $\\#$inv. el. \\\\ \n' );
fprintf(fid, '\\hline \n' );
% Initial mesh statistics IniIdeal
qualityValue = nodalDistortion(:,1);
fprintf(fid, '\\multicolumn{1}{|c|}{Ideal initial}       & %1.2f    &  %1.2f & %1.2f  & %1.2f    & %7.0f \\\\ \n',... 
    min(qualityValue), max(qualityValue) , sum(qualityValue)/length(qualityValue),...
    std(qualityValue), numDegIni);

% Initial mesh statistics IniIdeal
qualityValue = nodalDistortion(:,2);
fprintf(fid, '\\multicolumn{1}{|c|}{Ideal geometric}      & %1.2f    &  %1.2f & %1.2f  & %1.2f    & %7.0f \\\\ \n',... 
    min(qualityValue), max(qualityValue) , sum(qualityValue)/length(qualityValue),...
    std(qualityValue), numDegGeo);

fprintf(fid, '\\hline \n' );
fprintf(fid, '\\end{tabular} \n' );
fprintf(fid, '\\label{tab:%s} \n',nameFile );
fprintf(fid, '\\end{table} \n' );

end


