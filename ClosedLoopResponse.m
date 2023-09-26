classdef ClosedLoopResponse
    %CLOSELOOPRESPONSE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sltuner_objs slTuner
        reference_ap (1, 1) string
        plant_input_ap (1, 1) string
        system_output_ap (1, 1) string
    end

    properties(Constant)
        g = 9.80665;
        g_imperial = 32.174;
        t = -0.5:0.001:4;
    end
    
    methods
        function obj = ClosedLoopResponse(st, r, pi, so)
            obj.sltuner_objs = st;
            obj.reference_ap = r;
            obj.plant_input_ap = pi;
            obj.system_output_ap = so;
        end
        
        function plot(obj)
            fig = figure();
            t = tiledlayout(fig, 1, 3);
            t.TileSpacing = 'compact';

            % 1.Plot the step response of all the closed loop systems for a
            % 9g pull up
            
            ax1 = nexttile;
            obj.plot_pull_up(ax1);
            
            % 2. Plot the input deflection in degrees with a relevant step
            % at the input

            ax2 = nexttile;
            obj.plot_elevon_deflection(ax2);

            % 3. Plot CAP graph to show complience with the requirements
            wrap_get_response = @(tuner_obj)obj.get_tf(tuner_obj, obj.reference_ap, obj.system_output_ap);
            response_tf = arrayfun(wrap_get_response, obj.sltuner_objs, UniformOutput=false);
            ax3 = nexttile;
            obj.plot_cap(ax3, response_tf);


        end

        function ax = plot_pull_up(obj, ax)
            wrap_get_response = @(tuner_obj)obj.get_tf(tuner_obj, obj.reference_ap, obj.system_output_ap);
            response_tf = arrayfun(wrap_get_response, obj.sltuner_objs, UniformOutput=false);
            
            laod_factor_to_q = @(sys)obj.convert_to_load_q(sys, 9);
            required_q = cellfun(laod_factor_to_q, response_tf);
            
            wrap_gen_input = @(req_q)obj.gen_input(obj.t, req_q);
            u = arrayfun(wrap_gen_input, required_q, UniformOutput=false);

            wrap_lsim = @(sys, u)lsim(sys, u, obj.t);
            time_response = cellfun(wrap_lsim, response_tf, u, UniformOutput=false);
            
            wrap_plot = @(response)obj.plot_time_response(ax, response, obj.t);
            cellfun(wrap_plot, time_response);
            hold off
            grid(ax, "on");

            title(ax, "Pitch Rate Reponse with n=9");
            xlabel(ax, "Time [s]");
            ylabel(ax, "Pitch rate [rad/s]")
        end

        function ax = plot_elevon_deflection(obj, ax)
            wrap_get_response = @(tuner_obj)obj.get_tf(tuner_obj, obj.reference_ap, obj.plant_input_ap);
            response_tf = arrayfun(wrap_get_response, obj.sltuner_objs, UniformOutput=false);

            laod_factor_to_q = @(sys)obj.convert_to_load_q(sys, 9);
            required_q = cellfun(laod_factor_to_q, response_tf);
            
            wrap_gen_input = @(req_q)obj.gen_input(obj.t, req_q);
            u = arrayfun(wrap_gen_input, required_q, UniformOutput=false);

            wrap_lsim = @(sys, u)lsim(sys, u, obj.t);
            time_response = cellfun(wrap_lsim, response_tf, u, UniformOutput=false);
            
            time_response_deg = cellfun(@rad2deg, time_response, UniformOutput=false);


            wrap_plot = @(response)obj.plot_time_response(ax, response, obj.t);
            cellfun(wrap_plot, time_response_deg);
            hold off
            grid(ax, "on");


            title(ax, "Commanded Elevon deflection with n=9");
            xlabel(ax, "Time [s]");
            ylabel(ax, "Commanded Elevon deflection [rad]");

            % Add limits on deflections
            yline(ax, 25, 'black--', LineWidth=2);
            yline(ax, -30, 'black--', LineWidth=2);

        end
        


    end

    methods(Static)
        
        function ax = plot_cap(ax, response_tf)
            LOES = cellfun(@modelOrderReduction, response_tf, UniformOutput=false);
            
            zpk_loes = cellfun(@zpk, LOES, UniformOutput=false);

            [CAP, n_alpha, omega_sp, zeta_sp, t_theta_2] = cellfun(@ClosedLoopResponse.get_cap_param, zpk_loes);

            shortPeriodCategoryAPlot(ax, n_alpha, omega_sp, "bo");

        end

        function q = pitch_rate(velocity, load_factor)
            q = ClosedLoopResponse.g / velocity * (load_factor - 1);
        end

        function q = convert_to_load_q(sys, load_factor)
            grid_point = sys.SamplingGrid;
            
            [~, a, ~, ~, ~] = atmosisa(grid_point.altitude);
            velocity = a * grid_point.mach;
            
            q = ClosedLoopResponse.pitch_rate(velocity, load_factor);

        end

        function transfer_function = get_tf(sltuner_obj, input, output)
            transfer_function = getIOTransfer(sltuner_obj, input, output);
            transfer_function = tf(transfer_function);
        end

        function plot_time_response(ax, response, t)
            colors = {'b', 'r', 'g'};
            for i=1:size(response, 2)
                plot(ax, t, response(:, i), colors{i});
            end
            
            hold on
        end

        function u = gen_input(t, req_q)
            u = zeros(size(t));
            u(t>=0) = req_q;
        end
        
        function [CAP, n_alpha, omega_sp, zeta_sp, t_theta_2] = get_cap_param(loes_zpk)
            t_theta_2 = -1 / loes_zpk.Z{1};
            [omega_sp, zeta_sp] = damp(loes_zpk);
             
            if(omega_sp(1) ~= omega_sp(2))
                warning("Omega does not equal each other")
            end

            omega_sp = omega_sp(1);
            zeta_sp = zeta_sp(1);


            grid_point = loes_zpk.SamplingGrid;
            [~, a, ~, ~, ~] = atmosisa(grid_point.altitude);
            velocity = a * grid_point.mach;

            n_alpha = (velocity/ClosedLoopResponse.g) * (1/t_theta_2);
            CAP = omega_sp.^2 / n_alpha;

        end
    end
end

