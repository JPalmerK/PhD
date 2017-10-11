% Plot Pdet for 10 sites

close all; clear all; clc
csvdir=['W:\KJP PHD\6-Bayesian Habitat Use\Pdet at Time\2013\*PdetVals.csv']
csvfnames=dir(csvdir)

figure
subplot(2,1,1)
cc=jet(10);
order=[6 5 2 7 4 3 10 1 8 9];
nms={};
pdet=table();
for ii=1:length(csvfnames)

    pdetfname=['W:\KJP PHD\6-Bayesian Habitat Use\Pdet at Time\2013\', csvfnames(order(ii)).name];
    mm=readtable(pdetfname);
    mm=mm(mm.MedianPdet>0,:);
    mm=mm(mm.MatlabDate<7.3554*10^5,:);
    
    nms(ii)=cellstr(csvfnames(order(ii)).name(1:6));
    mm.DepLoc=repmat(nms{ii}, height(mm), 1);
    
    hold on
    plot(mm.MatlabDate, mm.MedianPdet, 'color', cc(ii,:), 'LineWidth', 1.5)
    pdet=[pdet; mm];
end

datetick('x','dd-mmm-yy','keepticks','keeplimits')
ylabel('Effective Detection Probability')

% Effective because it's the probability of detecting 50% of the clicks
legend(nms, 'Interpreter', 'none')
ax.XTickLabelRotation=45;
% set(findall(gcf,'type','text'),'FontSize',12,'FontName','TimesNewRoman')


subplot(2,1,2)
boxplot(pdet.MedianPdet,pdet.DepLoc)
ylabel('Effective Detection Probability', 'Interpreter', 'tex')

% Save the file as PNG
 pname={'Effective Detection Probability SM units 2013'};
 %print(cell2mat(pname),'-dpng','-r250')
 

 
 Pdet=readtable('W:\KJP PHD\3-Detection Function\Propagation Model\Area_vs_NL\Pdet 41 kHz oct.txt',...
     'ReadRowNames', 1);
figure
NLvals=[0:0.1:70];


for ii=1:10
    
    % dummy variable to place in correct order
    mm=order(ii);
    
    idxs=[3*mm-2:mm*3];
    subplot(3,4,ii)
    plot(NLvals, table2array(Pdet([idxs(3),idxs(1),idxs(2)],:))', 'LineWidth', 1.5)
    legend(Pdet.Properties.RowNames([idxs(3),idxs(1),idxs(2)]),...
        'Interpreter', 'none', 'Location','northeast', 'FontSize',9,'FontName','TimesNewRoman')
    title(names(ii))
    
    xlim([40 60])
    
    xlabel('NL_a_s_l', 'Interpreter', 'tex')
    
    ylim([0 0.1])
        if ii==1|| ii==5 || ii==9
        ylabel('Detection Probability')
        end
    
    %set(findall(gcf,'type','text'),'FontSize',12,'FontName','TimesNewRoman')
    set(findall(gcf,'FontSize',10,'FontName','TimesNewRoman'));
    % Save the file as PNG
    %pname=kk(ii);
    %print(cell2mat(pname),'-dpng','-r500')
end



%Plot Noise Levels

csvdir=['W:/KJP PHD/3-Detection Function/Propagation Model/Hourly NL 41khz thrdOctBand/*.csv']
csvfnames=dir(csvdir)

figure

subplot(2,1,1)
cc=jet(10);
NLs=table();
nms={};
order=[9 5 10 7 2 4 8 3 6 1];
names={'Latheron' 'Helmsdale', 'Cromarty', 'Spey Bay', 'Fraserburgh',...
    'Cruden Bay', 'Stonehaven', 'Arbroath', 'St Andrews', 'St Abbs'}

for ii=1:length(csvfnames)

    NLfname=[csvdir(1:end-5), csvfnames(order(ii)).name];
    mm=readtable(NLfname);
    mm=mm(mm.MatlabDate<7.3556*10^5,:);
    NLs=[NLs; mm];
    nms(ii)=mm.DepLoc(1)
    subplot(2,1,1)
    plot(mm.MatlabDate, mm.NLasl, 'color', cc(ii,:), 'LineWidth', 1)
    hold on

end
xlim([min(NLs.MatlabDate)-1 max(NLs.MatlabDate)+1])
datetick('x','dd-mmm-yy','keepticks','keeplimits')
ylabel('dB_A_S_L', 'Interpreter', 'tex')
legend(nms, 'Interpreter', 'none')
hold off
subplot(2,1,2)


% Make a boxplot of the NL values
boxplot(NLs.NLasl,NLs.DepLoc)
ylabel('dB_A_S_L', 'Interpreter', 'tex')

 pname={'Noise Levels from SM units 2013'};
 %print(cell2mat(pname),'-dpng','-r250')


