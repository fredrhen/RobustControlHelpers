classdef OpenLoopAnalysis
    %OPENLOOPANALYSIS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sltuner_objs slTuner
        plant_input_ap (1, 1) string
        plant_output_ap (1, 1) string
    end

    properties(Constant)
        freq_to_plot = logspace(-2, 3, 100);
    end
    
    methods
        function obj = OpenLoopAnalysis(st, pi, po)
            obj.sltuner_objs = st;
            obj.plant_input_ap = pi;
            obj.plant_output_ap = po;
        end
        
        function plot(obj)
            fig = figure();
            t = tiledlayout(fig, 2, 4);
            t.TileSpacing = 'compact';

            % 1. Nyquist plot at the plant input
            ax1 = nexttile;
            obj.my_nyquist_plot(ax1, obj.plant_input_ap);
            title(ax1, "Nyquist Plot @ Plant Input");
            % 2. Plot the minimum gain margin at plant input
            ax2 = nexttile; 
            obj.my_disk_gain_margin_plot(ax2, obj.plant_input_ap);
            title(ax2, "Min Disk Gain Margin @ Plant Input");

            % 3. Plot the minimum phase margin at plant input
            ax3 = nexttile;
            obj.my_disk_phase_margin_plot(ax3, obj.plant_input_ap);
            title(ax3, "Min Disk Phase Margin @ Plant Input");
            
            ax4 = nexttile;
            obj.loop_gain(ax4, obj.plant_input_ap);
            title(ax4, "Input Loop Transfer Function");

            % 4. Nyquist plot at the plant output
            ax5 = nexttile;
            obj.my_nyquist_plot(ax5, obj.plant_output_ap);
            title(ax5, "Nyquist Plot @ Plant Output");

            % 5. Plot the minimum gain margin at plant output
            ax6 = nexttile;
            obj.my_disk_gain_margin_plot(ax6, obj.plant_output_ap);
            title(ax6, "Min Disk Gain Margin @ Plant Output");

            % 6. Plot the minimum phase margin at plant output
            ax7 = nexttile;
            obj.my_disk_phase_margin_plot(ax7, obj.plant_output_ap);
            title(ax7, "Min Disk Phase Margin @ Plant Output");

            ax8 = nexttile;
            obj.loop_gain(ax8, obj.plant_input_ap);
            title(ax8, "Output Loop Transfer Function");

        end

        function plot_margins(obj)
            fig = figure();
            t = tiledlayout(fig, 2, 2);
            t.TileSpacing = 'compact';

            ax1 = nexttile;
            obj.my_gain_margin_plot(ax1, obj.plant_input_ap);
            title(ax1, "Min Gain Margin @ Plant Input");
            
            ax2 = nexttile; 
            obj.my_disk_gain_margin_plot(ax2, obj.plant_input_ap);
            title(ax2, "Min Disk Gain Margin @ Plant Input");

            ax3 = nexttile;
            obj.my_phase_margin_plot(ax3, obj.plant_input_ap);
            title(ax3, "Min Phase Margin @ Plant Input");

            ax4 = nexttile;
            obj.my_disk_phase_margin_plot(ax4, obj.plant_input_ap);
            title(ax4, "Min Disk Phase Margin @ Plant Input");
        end

        function plot_margin_save(obj, folder)
            fig = figure();
            ax1 = axes(fig);
            obj.my_gain_margin_plot(ax1, obj.plant_input_ap);
            exportgraphics(fig, folder+"gain_margin.pdf", Resolution=450);
            
            fig = figure();
            ax1 = axes(fig);
            obj.my_disk_gain_margin_plot(ax1, obj.plant_input_ap);
            exportgraphics(fig, folder+"disk_gain_margin.pdf", Resolution=450);

            fig = figure();
            ax1 = axes(fig);
            obj.my_phase_margin_plot(ax1, obj.plant_input_ap);
            exportgraphics(fig, folder+"phase_margin.pdf", Resolution=450);

            fig = figure();
            ax1 = axes(fig);
            obj.my_disk_phase_margin_plot(ax1, obj.plant_input_ap);
            exportgraphics(fig, folder+"disk_phase_margin.pdf", Resolution=450);
        end



        function my_nyquist_plot(obj, ax, location)
            % Solve for the size of the analysis point
            temp_loop_transfer = getLoopTransfer(obj.sltuner_objs(1), location, -1);
            ap_dim = size(temp_loop_transfer, 1);
            
            for i=1:ap_dim
                sub_location = location + sprintf("(%i)", i);
                wrap_get_LoopTransfer = @(st)getLoopTransfer(st, sub_location, -1);
                loop_transfer{i} = arrayfun(wrap_get_LoopTransfer, obj.sltuner_objs, UniformOutput=false);

                wrap_plot_nyquist = @(sys)obj.plot_nyquist(ax, sys);
                cellfun(wrap_plot_nyquist, loop_transfer{i}, UniformOutput=false);
            end
            
            xlim(ax, [-5, 5]);
            ylim(ax, [-5, 5]);
            % Plot the unit circle
            [x, y] = obj.circle(0, 0, 1);
            plot(ax, x, y, 'black--', LineWidth=1);
            
            plot(ax, -1, 0, 'r+');
            
            %Plot the lowest gain margin point
            % [Gm, Pm, Wcg, Wcp] = cellfun(@margin, sensitivity);
            % sensitivity = cellfun(@ss, sensitivity, UniformOutput=false);
            % Gm = squeeze(reshape(Gm, 1, []));
            % 
            % [minGm, minGmI] = min(abs(mag2db(Gm)));
            % freq = Wcg(minGmI);
            % 
            % corosponing_sens = sensitivity(minGmI);
            % H = freqresp(corosponing_sens{:}, freq);
            % X = real(H);
            % Y = imag(H);
            % 
            % plot(ax, X, Y, 'g.', MarkerSize=12);
            % 
            % % Plot the lowest phase margin
            % Pm = squeeze(reshape(Pm, 1, []));
            % 
            % [minPm, minPmI] = min(abs(Pm));
            % freq = Wcp(minPmI);
            % 
            % corosponing_sens = sensitivity(minPmI);
            % H = freqresp(corosponing_sens{:}, freq);
            % X = real(H);
            % Y = imag(H);
            % 
            % plot(ax, X, Y, 'm.', MarkerSize=12);
            % 
            % graph_text = sprintf("Min Gain Margin = %0.2f dB\nMin Phase Margin = %0.1f deg", minGm, minPm);
            % text(ax, 0.5, -0.25, graph_text);

            hold off
            grid(ax, 'on');
        end

        function my_gain_margin_plot(obj, ax, location)
            temp_loop_transfer = getLoopTransfer(obj.sltuner_objs(1), location, -1);
            ap_dim = size(temp_loop_transfer, 1);
            domain = temp_loop_transfer.SamplingGrid;
            
            for i=1:ap_dim
                sub_location = location + sprintf("(%i)", i);
                wrap_get_LoopTransfer = @(st)getLoopTransfer(st, sub_location, -1);
                loop_transfer{i} = arrayfun(wrap_get_LoopTransfer, obj.sltuner_objs, UniformOutput=false);

                [Gm(i, :, :), ~, ~, ~] = cellfun(@margin, loop_transfer{i});
                Gm(i, :, :) = abs(mag2db(Gm(i, :, :)));
            end

            wrap_get_domain = @(lt)lt.SamplingGrid;

            domain = cellfun(wrap_get_domain, loop_transfer{1});

            domain = cellfun(@(x)x.SamplingGrid, loop_transfer{1});
            alt_grid = arrayfun(@(x)x.altitude, domain);
            mach_grid = arrayfun(@(x)x.mach, domain);
            
            min_Gm = min(Gm, [], 1);
            min_Gm = squeeze(min_Gm);
            surf(alt_grid, mach_grid, min_Gm);

            xlabel(ax, "Altitude [m]");
            ylabel(ax, "Mach []");
            zlabel(ax, "Gain Margin [dB]");
        end

        function my_phase_margin_plot(obj, ax, location)
            temp_loop_transfer = getLoopTransfer(obj.sltuner_objs(1), location, -1);
            ap_dim = size(temp_loop_transfer, 1);
            
            for i=1:ap_dim
                sub_location = location + sprintf("(%i)", i);
                wrap_get_LoopTransfer = @(st)getLoopTransfer(st, sub_location, -1);
                loop_transfer{i} = arrayfun(wrap_get_LoopTransfer, obj.sltuner_objs, UniformOutput=false);

                [~, Pm(i, :, :), ~, ~] = cellfun(@margin, loop_transfer{i});
                Pm(i, :, :) = abs(Pm(i, :, :));
            end
            domain = cellfun(@(x)x.SamplingGrid, loop_transfer{1});
            alt_grid = arrayfun(@(x)x.altitude, domain);
            mach_grid = arrayfun(@(x)x.mach, domain);

            min_Pm = min(Pm, [], 1);
            min_Pm = squeeze(min_Pm);
            surf(ax, alt_grid, mach_grid, min_Pm);
            
            xlabel(ax, "Altitude [m]");
            ylabel(ax, "Mach []");
            zlabel(ax, "Phase Margin [deg]");

        end
        
        function my_disk_gain_margin_plot(obj, ax, location)
            temp_loop_transfer = getLoopTransfer(obj.sltuner_objs(1), location, -1);
            ap_dim = size(temp_loop_transfer, 1);
            
            for i=1:ap_dim
                sub_location = location + sprintf("(%i)", i);
                wrap_get_LoopTransfer = @(st)getLoopTransfer(st, sub_location, -1);
                loop_transfer{i} = arrayfun(wrap_get_LoopTransfer, obj.sltuner_objs, UniformOutput=false);
                
                [DM(i, :, :), MM(i, :, :)] = cellfun(@diskmargin, loop_transfer{i});
               
            end
            
            GM = arrayfun(@(st)st.GainMargin, DM,UniformOutput=false);
            GM = cellfun(@(element)min(abs(mag2db(element))), GM);
            GM = min(GM, [], 1);
            GM = squeeze(GM);

            % Get the sampling grid points
            domain = cellfun(@(x)x.SamplingGrid, loop_transfer{1});
            alt_grid = arrayfun(@(x)x.altitude, domain);
            mach_grid = arrayfun(@(x)x.mach, domain);

            surf(ax, alt_grid, mach_grid, GM);

            xlabel(ax, "Altitude [m]");
            ylabel(ax, "Mach []");
            zlabel(ax, "Gain Margin [dB]");
        end

        function my_disk_phase_margin_plot(obj, ax, location)
            temp_loop_transfer = getLoopTransfer(obj.sltuner_objs(1), location, -1);
            ap_dim = size(temp_loop_transfer, 1);
            
            for i=1:ap_dim
                sub_location = location + sprintf("(%i)", i);
                wrap_get_LoopTransfer = @(st)getLoopTransfer(st, sub_location, -1);
                loop_transfer{i} = arrayfun(wrap_get_LoopTransfer, obj.sltuner_objs, UniformOutput=false);

                [DM(i, :, :), MM(i, :, :)] = cellfun(@diskmargin, loop_transfer{i});
            end
            
            PM = arrayfun(@(st)st.PhaseMargin, DM,UniformOutput=false);
            PM = cellfun(@(element)min(abs(element)), PM);
            PM = min(PM, [], 1);
            PM = squeeze(PM);
            
            % Get the sampling grid points
            domain = cellfun(@(x)x.SamplingGrid, loop_transfer{1});
            alt_grid = arrayfun(@(x)x.altitude, domain);
            mach_grid = arrayfun(@(x)x.mach, domain);

            surf(ax, alt_grid, mach_grid, PM);
            
            xlabel(ax, "Altitude [m]");
            ylabel(ax, "Mach []");
            zlabel(ax, "Phase Margin [deg]");
        end

        function loop_gain(obj, ax, location)
            temp_loop_transfer = getLoopTransfer(obj.sltuner_objs(1), location, -1);
            ap_dim = size(temp_loop_transfer, 1);

            for i=1:ap_dim
                sub_location = location + sprintf("(%i)", i);
                wrap_get_LoopTransfer = @(st)getLoopTransfer(st, sub_location, -1);
                loop_transfer{i} = arrayfun(wrap_get_LoopTransfer, obj.sltuner_objs, UniformOutput=false);
                obj.plot_tf(ax, loop_transfer{i});
            end
            grid on
        end
    
        function my_nichols_plot(obj, ax, location)
            temp_loop_transfer = getLoopTransfer(obj.sltuner_objs(1), location, -1);
            ap_dim = size(temp_loop_transfer, 1);
            
            for i=1:ap_dim
                sub_location = location + sprintf("(%i)", i);
                wrap_get_LoopTransfer = @(st)getLoopTransfer(st, sub_location, -1);
                loop_transfer{i} = arrayfun(wrap_get_LoopTransfer, obj.sltuner_objs, UniformOutput=false);

                wrap_plot_nichols = @(sys)obj.plot_nichols(ax, sys);
                cellfun(wrap_plot_nichols, loop_transfer{i}, UniformOutput=false);
            end

            % Plot the exlusion zone
            exclusion = polyshape([-145, -180, -180, -145], [3, 6, -6, -3]);
            
            exclusion_plot = plot(ax, exclusion);
            exclusion_plot.FaceAlpha = 0;
            xlim(ax, [-200 -100]);
            ylim(ax, [-10, 10]);

            grid(ax, 'on');
            
            txt = {'Nichols' 'Exclusion'};
            t = text(ax, -163, 0, txt, HorizontalAlignment="center");
            t.FontSize = 14;

            ylabel(ax, "Gain [dB]");
            xlabel(ax, "Phase [deg]");

        end

    end



    methods(Static)
        function plot_nyquist(ax, system)
            for i =1:size(system, 1)
                for j=1:size(system, 2)
                    [re, im, ~] = nyquist(system(:, :, i, j));
                    temp_re = squeeze(re);
                    temp_im = squeeze(im);
                    plot(ax, temp_re, temp_im, 'b');
                    plot(ax, temp_re, -temp_im, 'b');
                    hold on
                end
            end
        end

        function plot_nichols(ax, system)
            for i =1:size(system, 1)
                for j=1:size(system, 2)
                    [mag, phase, ~] = nichols(system(:, :, i, j));
                    temp_mag = squeeze(mag);
                    temp_phase = squeeze(phase);
                    plot(ax, temp_phase, mag2db(temp_mag), 'b');
                    hold on
                end
            end
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

        function ax = plot_tf(ax, tfs)            
            wrap_freqresp = @(sys)sigma(sys, GangOfSix.freq_to_plot);
            H = cellfun(wrap_freqresp, tfs, UniformOutput=false);

            H = squeeze(vertcat(H{:}));
            H = mag2db(H);
            semilogx(ax, GangOfSix.freq_to_plot, H, 'b');
        end 
    end 
end

