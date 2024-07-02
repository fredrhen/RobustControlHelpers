classdef ClosedLoopResponse
    %CLOSELOOPRESPONSE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ref_to_output 
        ref_to_plant_input
    end

    properties(Constant)
        g = 9.80665;
        g_imperial = 32.174;
        t = -0.5:0.01:8;
    end
    
    methods
        function obj = ClosedLoopResponse(st, r, pi, so)
            if nargin == 4
                wrap_ref_to_output = @(tuner_obj)obj.get_tf(tuner_obj, r, so);
                ref_to_output = arrayfun(wrap_ref_to_output, st, UniformOutput=false);
        
                wrap_ref_to_plant_input = @(tuner_obj)obj.get_tf(tuner_obj, r, pi);
                ref_to_plant_input = arrayfun(wrap_ref_to_plant_input, st, UniformOutput=false);
                
                obj.ref_to_output = cell2MultiModel(ref_to_output);
                obj.ref_to_plant_input = cell2MultiModel(ref_to_plant_input);
            else
                obj.ref_to_output = ss();
                obj.ref_to_plant_input = ss();
            end
        end
        
        function plot(obj, ref_model)
            fig = figure();
            t = tiledlayout(fig, 2, 2);
            t.TileSpacing = 'compact';

            % 1.Plot the step response of all the closed loop systems for a
            % 9g pull up
            
            ax1 = nexttile;
            obj.plot_pull_up(ax1, ref_model);
            
            % 2. Plot the input deflection in degrees with a relevant step
            % at the input

            ax2 = nexttile;
            obj.plot_elevon_deflection(ax2);

            ax3 = nexttile;
            obj.plot_elevon_deflection_rate(ax3);



            % 3. Plot CAP graph to show complience with the requirements


        end

        function ax = plot_pull_up(obj, ax, ref_model)
            required_q = obj.convert_to_load_q(obj.ref_to_output, 9);
            
            one_array = ones(size(required_q));

            wrap_gen_input = @(req_q)obj.gen_input(obj.t, req_q);
            u = arrayfun(wrap_gen_input, one_array, UniformOutput=false);
            
            ref_to_output_cell = multiModel2Cell(obj.ref_to_output);
            u_cell = cell(u);
            wrap_lsim = @(sys, u)lsim(sys, u, obj.t);
            time_response = cellfun(wrap_lsim, ref_to_output_cell, u_cell, UniformOutput=false);
            
            ref_model_response = lsim(ref_model, u{1}, obj.t);

            wrap_plot = @(response)obj.plot_time_response(ax, response, obj.t, '-b');
            cellfun(wrap_plot, time_response);
            obj.plot_time_response(ax, ref_model_response, obj.t, '-r')
            hold off
            grid(ax, "on");

            title(ax, "Pitch Rate Reponse with n=9");
            xlabel(ax, "Time [s]");
            ylabel(ax, "Pitch rate [rad/s]")
        end

        function ax = plot_elevon_deflection(obj, ax)
            required_q = obj.convert_to_load_q(obj.ref_to_plant_input, 9);
            
            wrap_gen_input = @(req_q)obj.gen_input(obj.t, req_q);
            u = arrayfun(wrap_gen_input, required_q, UniformOutput=false);
            u = cell(u);
            
            ref_to_input_cell = multiModel2Cell(obj.ref_to_plant_input);

            wrap_lsim = @(sys, u)lsim(sys, u, obj.t);
            time_response = cellfun(wrap_lsim, ref_to_input_cell, u, UniformOutput=false);
            
            time_response_deg = cellfun(@rad2deg, time_response, UniformOutput=false);


            wrap_plot = @(response)obj.plot_time_response(ax, response, obj.t, '-b');
            cellfun(wrap_plot, time_response_deg);
            hold off
            grid(ax, "on");


            title(ax, "Commanded Elevon deflection with n=9");
            xlabel(ax, "Time [s]");
            ylabel(ax, "Commanded Elevon deflection [deg]");

            % Add limits on deflections
            yline(ax, 25, 'black--', LineWidth=2);
            yline(ax, -30, 'black--', LineWidth=2);

        end

        function ax = plot_elevon_deflection_rate(obj, ax)
            required_q = obj.convert_to_load_q(obj.ref_to_plant_input, 9);
            
            wrap_gen_input = @(req_q)obj.gen_input(obj.t, req_q);
            u = arrayfun(wrap_gen_input, required_q, UniformOutput=false);
            u = cell(u);
           
            s = tf('s');
            ref_to_input_rate = obj.ref_to_plant_input * s;

            ref_to_input_cell = multiModel2Cell(ref_to_input_rate);
            
            wrap_lsim = @(sys, u)lsim(sys, u, obj.t);
            time_response = cellfun(wrap_lsim, ref_to_input_cell, u, UniformOutput=false);
            
            time_response_deg = cellfun(@rad2deg, time_response, UniformOutput=false);


            wrap_plot = @(response)obj.plot_time_response(ax, response, obj.t, '-b');
            cellfun(wrap_plot, time_response_deg);
            hold off
            grid(ax, "on");


            title(ax, "Commanded Elevon deflection with n=9");
            xlabel(ax, "Time [s]");
            ylabel(ax, "Commanded Elevon deflection rate [deg/s]");

            % Add limits on deflections
            yline(ax, 50, 'black--', LineWidth=2);
            yline(ax, -50, 'black--', LineWidth=2);

        end
        


    end

    methods(Static)

        function q = pitch_rate(velocity, load_factor)
            q = ClosedLoopResponse.g ./ velocity .* (load_factor - 1);
        end

        function q = convert_to_load_q(sys, load_factor)
            grid_point = sys.SamplingGrid;
            
            [~, a, ~, ~, ~] = atmosisa(grid_point.altitude);
            velocity = a .* grid_point.mach;
            
            q = ClosedLoopResponse.pitch_rate(velocity, load_factor);

        end

        function transfer_function = get_tf(sltuner_obj, input, output)
            transfer_function = getIOTransfer(sltuner_obj, input, output);
            transfer_function = tf(transfer_function);
        end

        function plot_time_response(ax, response, t, color)
            colors = {'b', 'r', 'g'};
            for i=1:size(response, 2)
                plot(ax, t, response(:, i), color);
            end
            
            hold on
        end

        function u = gen_input(t, req_q)
            u = zeros(size(t));
            u(t>=0) = req_q;
        end
        
    end
end

