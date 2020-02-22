% Poisson spike train statistics
h_f = figure; 
h_a = axes;

m_t = 20e-3; %in sec

for i = 1:10
    t_exp = exprnd(m_t*ones(1,50));
    t_vect = cumsum(t_exp);
    i_vect = find(t_vect < 0.5); %plot 500 ms
    t_plot = t_vect(i_vect);

    line([t_plot; t_plot],[i i+0.9]'*ones(size(t_plot)),'Parent',h_a,'Color','k');
end;

axes(h_a);
set(h_a,'XLim',[-0.01 0.51],'XTick',[0 0.1 0.2 0.3 0.4 0.5]);
xlabel('time (s)');

