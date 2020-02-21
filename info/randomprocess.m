%{Brownian motion, Gaussian stochastic process
Three sample paths of a Wiener process. These paths were obtained by summing up and scaling by delta_t the white noise paths in C. B. 
Sample mean and mean ± one standard deviation of 100 Wiener sample paths. Note the time dependence of the standard deviation, implying that 
the Wiener process is not stationary. C. Three sample paths of white noise sampled at a time step deltat = 1 ms. 

For clarity, two of the paths have been shifted above and below the horizontal zero line (dashed). D. Mean and mean ± one standard deviation of 
100 white noise samples (solid and dotted lines, respectively).
%}

% Generate sample paths of white noise and Wiener process
dt = 1e-3; %time step in ms

% Corresponding standard deviation of a white noise sample
s_wn = 1/sqrt(dt);

%generate 100 white noise samples at scale dt, 100*dt s long
wn_paths = normrnd(0,s_wn,100,100);

%corresponding wiener sample paths. Because the white noise has a variance of 1/dt, 
%we need to multiply by dt to obtain a variance of dt.
wiener_paths = [zeros(100,1) cumsum(wn_paths,2)*dt];

colors = {[1 0 0]; [0.5 0 0]; [0 0 0]};

h_f1 = figure;
h_a1 = subplot(4,1,1);
h_a2 = subplot(4,1,2);
h_a3 = subplot(4,1,3);
h_a4 = subplot(4,1,4);

h_f2 = figure; 
h_a5 = subplot(2,1,1);
h_a6 = subplot(2,1,2);

for i=1:3
    line('Parent', eval(sprintf('h_a%i',i)), ...
        'XData',[1:100]*dt,'YData',wn_paths(i,:),'Color',colors{i});
    
    line('Parent', h_a5,'XData',[0:100]*dt,'YData',wiener_paths(i,:),'Color',colors{i});
end

axes(h_a4);
xlabel('time (s)');
axes(h_a6);
xlabel('time (s)');

%%
mwn_path = mean(wn_paths,1);
swn_path = std(wn_paths,1,1);

mwiener_path = mean(wiener_paths,1);
swiener_path = std(wiener_paths,1,1); 

line('Parent', h_a4,'XData',[1:100]*dt,'YData',mwn_path);
line('Parent', h_a4,'XData',[1:100]*dt,'YData',mwn_path+swn_path,...
     'LineStyle',':');
line('Parent', h_a4,'XData',[1:100]*dt,'YData',mwn_path-swn_path,...
     'LineStyle',':');

set(h_a4,'YLim',[-1e2 1e2]);
set(h_a5,'YLim',[-0.8 0.8]);

line('Parent', h_a6,'XData',[0:100]*dt,'YData',mwiener_path);
line('Parent', h_a6,'XData',[0:100]*dt,'YData',mwiener_path+swiener_path,...
     'LineStyle',':');
line('Parent', h_a6,'XData',[0:100]*dt,'YData',mwiener_path-swiener_path,...
     'LineStyle',':');
 
%print(handles.figure1,'-depsc','figure1.eps'); 

