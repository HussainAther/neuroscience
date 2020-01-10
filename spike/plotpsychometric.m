function h=plot_psychometric(D,varargin)
% h=plot_psychometric(D,varargin)
%
% Plots psychometric data
% D structure contains data: strength, choice, rt, corr, lcorr [correct from
% logistic]
% this can be real or predicted data
% 
% The expected varargin consist of 
% data      (0 [1]) plot as data i.e. points siwth s.e
% pred      (0 [1]) plot as predictions i.e. line
% folded    (0 [1]) plot folded
% rt        (0 [1]) plot rt
% correct   ([0] 1) plot only correct rt
% choice    (0 [1]) plot choice
% split0     ([0] 1) split0 rt at zero to allow different tnd
% xsplit    (0 [1]) split0 x axis near 0 
% col        "k"    plot color

pSet = inputParser;
gms=get(groot,"defaultLineMarkerSize");
glw=get(groot,"defaultlinelinewidth");

addParameter(pSet,"data",1); % plot  predictions 
addParameter(pSet,"pred",1);  % plot  predictions 
addParameter(pSet,"folded",1);  % plot folded
addParameter(pSet,"rt",1);  %  plot rt
addParameter(pSet,"choice",1);  % plot choice
addParameter(pSet,"split0",0);  % split0 rt at zero to allow differnet tnd
addParameter(pSet,"col","k"); % color for plot
addParameter(pSet,"correct",0); % color for plot
addParameter(pSet,"xsplit",1); % color for plot

parse(pSet,varargin{:});
v2struct(pSet.Results);


%done unpacking
D.uastrength=unique(abs(D.strength));
x1=unique(abs(D.strength));
x1(x1==0)=min(x1(x1>0))/2;
x2=unique(abs(D.strength));
x2(x2==0)=min(x1);
dd=1.2;
D.ustrength=unique(D.strength);

if rt & choice
    twoplot=1;
else
    twoplot=0;
end

if folded
    if twoplot
        subplot(1,2,1)
    end
    S=bindata(abs(D.strength),D.corr);
    
    if choice
        hold on
        if data
            S.m(1,1)=0.5;
            errbar(x1,S.m,S.se,"b","linewidth",1,"color",col)
            h=plot(x1,S.m,"o","MarkerFace",col,"MarkerSize",gms,"MarkerEdgeColor",col);
            set(h,"linewidth",glw);
        end
        
        if pred
            plot([x2(1)/dd dd*x2(1)],[S.m(1) S.m(1)],col)
            h=plot(x2(2:end),S.m(2:end),"-","Color",col);
        end
        
        xlabel("Motion strength (%coherence)")
        ylabel("Proportion correct")
        set(gca,"Xscale","log")
        set(gca,"XTick",x1)
        set(gca,"XMinorTick","off")
        set(gca,"XTickLabel",num2cell(D.uastrength*100))
        if  xsplit
            xbreaklog(0.023,1.05)
        end
        
        box off
    end
    if twoplot
        
        subplot(1,2,2)
    end
    if rt
        if correct
            S=bindata(abs(D.strength(D.lcorr>0)),D.rt(D.lcorr>0));
        else
            S=bindata(abs(D.strength),D.rt);     
        end
        hold on
        
        if data
            errbar(x1,S.m,S.se,"b","linewidth",1,"Color",col)
            h=plot(x1,S.m,"o","MarkerFace",col,"MarkerSize",gms,"MarkerEdgeColor",col);
            set(h,"linewidth",glw);
        end    
        
        if pred
            plot([x2(1)/dd dd*x2(1)],[S.m(1) S.m(1)],col)
            h= plot(x2(2:end),S.m(2:end),"-","Color",col);
        end
        
        xlabel("Motion strength (%coherence)")
        ylabel("Reaction time (s)")
        set(gca,"Xscale","log")
        set(gca,"XTick",x1)
        set(gca,"XMinorTick","off")
        
        set(gca,"XTickLabel",num2cell(D.uastrength*100))
        if  xsplit
            xbreaklog(0.023,1.05)
        end
        box off
    end
    
    
else
    if choice
        D.ustrength=unique([D.uastrength -D.uastrength]);
        if twoplot
            subplot(1,2,1)
        end
        hold on
        S=bindata(D.strength,D.choice);
        
        if data% here xxxxx
            errbar(S.x,S.m,S.se,"b","linewidth",1,"color",col)
            h=plot(S.x,S.m,"o","MarkerFace",col,"MarkerSize",gms,"Color",col,"MarkerEdgeColor",col);

        end
        
        if pred
            h=plot(S.x,S.m,"-","Color",col);

            
        end
        
        xlabel("Motion strength (%strength)")
        ylabel("P_{right}")
        set(gca,"XMinorTick","off")
        set(gca,"XTick",[min(D.strength) 0 max(D.strength)])
        set(gca,"XTickLabel",num2cell([min(D.strength) 0 max(D.strength)]*100))
        box off
    end
    
    
    if twoplot
        subplot(1,2,2)
    end
    if rt
        hold on
        if correct
            S=bindata(D.strength(D.lcorr>0),D.rt(D.lcorr>0));
        else
            S=bindata(D.strength,D.rt);     
        end
        if data
            errbar(S.x,S.m,S.se,"b","linewidth",1,"Color",col,"Marker","o")
            h=plot(S.x,S.m,"o","MarkerFace",col,"MarkerSize",gms,"Color",col,"MarkerEdgeColor",col);
            set(h,"linewidth",glw)
        end
        
        if pred
            if split0
                i=find(S.x==0)
                S.m(i-1)=NaN;
            end
            h=plot(S.x,S.m,"-","Color",col);
        end
        xlabel("Motion strength (%strength)")
        ylabel("Reaction time (s)")
        set(gca,"XTick",[min(D.strength) 0 max(D.strength)])
        set(gca,"XMinorTick","off")
        set(gca,"XTickLabel",num2cell([min(D.strength) 0 max(D.strength)]*100))
        box off
    end
    
end