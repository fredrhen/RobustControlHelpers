classdef OpenLoopAnalysis
    %OPENLOOPANALYSIS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sltuner_objs slTuner
        plant_input_ap (1, 1) string
        plant_output_ap (1, 1) string
    end
    
    methods
        function obj = OpenLoopAnalysis(st, pi, po)
            obj.sltuner_objs = st;
            obj.plant_input_ap = pi;
            obj.plant_output_ap = po;
        end
        
        function plot(obj)
            fig = figure();
            t = tiledlayout(fig, 2, 3);
            t.TileSpacing = 'compact';

            % 1. Nyquist plot at the plant input
            ax1 = nexttile;
            obj.my_nyquist_plot(ax1, obj.plant_input_ap);
            title(ax1, "Nyquist Plot @ Plant Input");
            % 2. Plot the minimum gain margin at plant input
            ax2 = nexttile; 
            obj.my_gain_margin_plot(ax2, obj.plant_input_ap);
            title(ax2, "Gain Margin @ Plant Input");

            % 3. Plot the minimum phase margin at plant input
            ax3 = nexttile;
            obj.my_phase_margin_plot(ax3, obj.plant_input_ap);
            title(ax3, "Phase Margin @ Plant Input");

            % 4. Nyquist plot at the plant output
            ax4 = nexttile;
            obj.my_nyquist_plot(ax4, obj.plant_output_ap);
            title(ax4, "Nyquist Plot @ Plant Output");

            % 5. Plot the minimum gain margin at plant output
            ax5 = nexttile;
            obj.my_gain_margin_plot(ax5, obj.plant_output_ap);
            title(ax5, "Gain Margin @ Plant Output");

            % 6. Plot the minimum phase margin at plant output
            ax6 = nexttile;
            obj.my_phase_margin_plot(ax6, obj.plant_output_ap);
            title(ax6, "Phase Margin @ Plant Output");

        end

        function my_nyquist_plot(obj, ax, location)
            wrap_get_senitivity = @(st)getSensitivity(st, location);
            sensitivity = arrayfun(wrap_get_senitivity, obj.sltuner_objs, UniformOutput=false);


            
            wrap_plot_nyquist = @(sys)obj.plot_nyquist(ax, sys);
            cellfun(wrap_plot_nyquist, sensitivity, UniformOutput=false);

            % Plot the unit circle
            [x, y] = obj.circle(-1, 0, 1);
            plot(ax, x, y, 'black--', LineWidth=1);
            
            plot(ax, -1, 0, 'r+');
            
            %Plot the lowest gain margin point
            [Gm, Pm, Wcg, Wcp] = cellfun(@margin, sensitivity);
            sensitivity = cellfun(@ss, sensitivity, UniformOutput=false);
            Gm = squeeze(reshape(Gm, 1, []));

            [minGm, minGmI] = min(abs(mag2db(Gm)));
            freq = Wcg(minGmI);
            
            corosponing_sens = sensitivity(minGmI);
            H = freqresp(corosponing_sens{:}, freq);
            X = real(H);
            Y = imag(H);
            
            plot(ax, X, Y, 'g.', MarkerSize=12);

            % Plot the lowest phase margin
            Pm = squeeze(reshape(Pm, 1, []));

            [minPm, minPmI] = min(abs(Pm));
            freq = Wcp(minPmI);
            
            corosponing_sens = sensitivity(minPmI);
            H = freqresp(corosponing_sens{:}, freq);
            X = real(H);
            Y = imag(H);
            
            plot(ax, X, Y, 'm.', MarkerSize=12);
            
            graph_text = sprintf("Min Gain Margin = %0.2f dB\nMin Phase Margin = %0.1f deg", minGm, minPm);
            text(ax, 0.5, -0.25, graph_text);

            hold off
            grid(ax, 'on');
        end

        function my_gain_margin_plot(obj, ax, location)
            wrap_get_senitivity = @(st)getSensitivity(st, location);
            sensitivity = arrayfun(wrap_get_senitivity, obj.sltuner_objs, UniformOutput=false);

            [Gm, ~, ~, ~] = cellfun(@margin, sensitivity);
            Gm = abs(mag2db(Gm));
            
            surf(ax, Gm);

            xlabel(ax, "Altitude [m]");
            ylabel(ax, "Mach []");
            zlabel(ax, "Gain Margin [dB]");
        end

        function my_phase_margin_plot(obj, ax, location)
            wrap_get_senitivity = @(st)getSensitivity(st, location);
            sensitivity = arrayfun(wrap_get_senitivity, obj.sltuner_objs, UniformOutput=false);

            [~, Pm, ~, ~] = cellfun(@margin, sensitivity);
            Pm = abs(Pm);
            surf(ax, Pm);
            
            xlabel(ax, "Altitude [m]");
            ylabel(ax, "Mach []");
            zlabel(ax, "Phase Margin [deg]");

        end
    end

    methods(Static)
        function plot_nyquist(ax, system)
            [re, im] = nyquist(system);
            re = squeeze(re);
            im = squeeze(im);

            plot(ax, re, im, 'b');
            hold on
        end

        function [xp, yp] = circle(x,y,r)
            %x and y are the coordinates of the center of the circle
            %r is the radius of the circle
            %0.01 is the angle step, bigger values will draw the circle faster but
            %you might notice imperfections (not very smooth)
            ang=0:0.01:2*pi; 
            xp=r*cos(ang) + x;
            yp=r*sin(ang) + y;
        end
    end 
end

