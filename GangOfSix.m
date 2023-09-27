classdef GangOfSix
    %GANGOFSIX Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sltuner_objs slTuner
        disturbance_input_ap (1, 1) string
        disturbance_output_ap (1, 1) string
        reference_ap (1, 1) string
        noise_ap (1, 1) string
        plant_input_ap (1, 1) string
        system_output_ap (1, 1) string
    end

    properties(Constant)
        freq_to_plot = logspace(-2, 3, 100);
    end
    
    methods
        function obj = GangOfSix(st, di, do, r, n, pi, so)
            obj.sltuner_objs = st;
            obj.disturbance_input_ap = di;
            obj.disturbance_output_ap = do;
            obj.reference_ap = r;
            obj.noise_ap = n;
            obj.plant_input_ap = pi;
            obj.system_output_ap = so;
        end
        
        function plot(obj)
            fig = figure();
            %t = tiledlayout(fig, 2, 3); $ Tiledlayout does not work with
            %sigma plot
            %t.TileSpacing = 'compact';

            % So plot
            ax1 = subplot(2, 3, 1);
            plot_So(obj, ax1);

            ax2 = subplot(2, 3, 2);
            plot_KSo(obj, ax2);

            ax3 = subplot(2, 3, 3);
            plot_To(obj, ax3);

            ax4 = subplot(2, 3, 4);
            plot_SoG(obj, ax4);

            ax5 = subplot(2, 3, 5);
            plot_KSoF(obj, ax5);

            ax6 = subplot(2, 3, 6);
            plot_ToF(obj, ax6);
            

        end

        function ax = plot_So(obj, ax)
            wrap_get_SO = @(tuner_obj)obj.get_tf(tuner_obj, obj.disturbance_output_ap, obj.system_output_ap, 1);
            So = arrayfun(wrap_get_SO, obj.sltuner_objs, UniformOutput=false);
            
            ax = obj.plot_tf(ax, So);
            title(ax, "So");
            grid(ax, "on");
        end

        function ax = plot_KSo(obj, ax)
            wrap_get_KSO = @(tuner_obj)obj.get_tf(tuner_obj, obj.disturbance_output_ap, obj.plant_input_ap, -1);
            So = arrayfun(wrap_get_KSO, obj.sltuner_objs, UniformOutput=false);
            
            ax = obj.plot_tf(ax, So);
            title(ax, "KSo");
            grid(ax, "on");
        end

        function ax = plot_To(obj, ax)
            wrap_get_To = @(tuner_obj)obj.get_tf(tuner_obj, obj.noise_ap, obj.system_output_ap, -1);
            So = arrayfun(wrap_get_To, obj.sltuner_objs, UniformOutput=false);
            
            ax = obj.plot_tf(ax, So);
            title(ax, "To");
            grid(ax, "on");
        end

        function ax = plot_SoG(obj, ax)
            wrap_get_SoG = @(tuner_obj)obj.get_tf(tuner_obj, obj.disturbance_input_ap, obj.system_output_ap, -1);
            So = arrayfun(wrap_get_SoG, obj.sltuner_objs, UniformOutput=false);
            
            ax = obj.plot_tf(ax, So);
            title(ax, "SoG");
            grid(ax, "on");
        end

        function ax = plot_KSoF(obj, ax)
            wrap_get_KSoF = @(tuner_obj)obj.get_tf(tuner_obj, obj.reference_ap, obj.plant_input_ap, -1);
            So = arrayfun(wrap_get_KSoF, obj.sltuner_objs, UniformOutput=false);
            
            ax = obj.plot_tf(ax, So);
            title(ax, "KSoF");
            grid(ax, "on");
        end

        function ax = plot_ToF(obj, ax)
            wrap_get_ToF = @(tuner_obj)obj.get_tf(tuner_obj, obj.reference_ap, obj.system_output_ap, -1);
            So = arrayfun(wrap_get_ToF, obj.sltuner_objs, UniformOutput=false);
            
            ax = obj.plot_tf(ax, So);
            title(ax, "ToF");
            grid(ax, "on");
        end

    end
 
    methods(Static)
        function ax = plot_tf(ax, tfs)
            hold on
            wrap_sigma_plot = @(sys)sigmaplot(ax, genss(sys));
            cellfun(wrap_sigma_plot, tfs);
            hold off
            
            % wrap_freqresp = @(sys)sigma(genss(sys), GangOfSix.freq_to_plot);
            % H = cellfun(wrap_freqresp, tfs, UniformOutput=false);
            % 
            % H = squeeze(vertcat(H{:}));
            % H = mag2db(H);
            % semilogx(ax, GangOfSix.freq_to_plot, H, 'b');
            % 
           
        end

        function transfer_function = get_tf(sltuner_obj, input, output, sign)
            transfer_function = getIOTransfer(sltuner_obj, input, output);
            transfer_function = sign*tf(transfer_function);
        end
    end
end

