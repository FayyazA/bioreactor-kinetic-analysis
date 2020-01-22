function exchange_fitting_1D_plot(spc_config, spc_params, ...
    S_data, S_fit, varargin)

% ---- set parameters ----
% addpath_matlab

Scolor = [37 64 97; 55 96 146;...
          152 72 7; 228 108 10;...
          79 98 40; 119 147 60]/256;

      
varargs_length = 4;
if length(varargin) > varargs_length
    error('myfuns:exchange_fitting_1D_plot:TooManyInputs', ...
        'requires at most %d optional inputs',varargs_length);
end

optargs = {1, 5};
% optargs = {plot_scale,t_offset}
% model_params is model-specific parameters
for i_args = 1:length(varargin),
    if ~isempty(varargin{i_args})
        optargs{i_args} = varargin{i_args};
    end
end

% ---- plot data and fit ----
plot_scale = optargs{1};
t_offset = optargs{2};
t = spc_params.timepoint_dt*(1:size(S_data,2))+t_offset;

f_h = figure,
hold on
for i_plot = 1:size(S_data,1),
    if i_plot == 1
        plot(t,S_data(i_plot,:)*plot_scale,'-','Color',Scolor(i_plot*2-1,:),'linewidth',2);
%         plot(t,S_fit(i_plot,:)*plot_scale,':','Color',Scolor(i_plot*2,:),'linewidth',2);
    else
        plot(t,S_data(i_plot,:),'-','Color',Scolor(i_plot*2-1,:),'linewidth',2);
        plot(t,S_fit(i_plot,:),':','Color',Scolor(i_plot*2,:),'linewidth',2);
    end
end
% title('1D model fit','FontName','Arial','fontsize',14);
xlabel('Time Post-injection (s)','FontName','Arial','fontsize',14)
ylabel('HP-^{13}C signal (AU)','FontName','Arial','fontsize',14);
if plot_scale == 1
    pyr_leg = 'Pyruvate (data)';
else
    pyr_leg = sprintf('Pyruvate (data)x%.1f',plot_scale);
end

    
legend(pyr_leg,'Lactate (data)','Lactate (fit)');
xlim([0 spc_params.timepoint_dt*size(S_data,2)+2]+t_offset);
ylim([-0.05*plot_scale*max(S_data(:)) 1.1*plot_scale*max(S_data(:))]);

set(f_h, 'position', [100 100 400 250])
% if length(varargin) >= 2,
% for i_text = 2:length(varargin),
%         text(spc_params.timepoint_dt,max(S_data(:))*(1.2-0.1*i_text),...
%             sprintf('%s = %.5f',optargs.name{i_text},eval(sprintf('%s',optargs.name{i_text}))),...
%             'FontName','Arial','fontsize',14);
% end
% end
