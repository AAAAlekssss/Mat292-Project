%Parameters
Blood_Fat_Exchange = 0.195;
Blood_Reaction_Exchange = 0.085;
Blood_Losses = 0.0875;
Reaction_Losses = 0.07;

%initial Value Problem Characteristics
Initial_Value = 0.5;
End_Time = 20;
Time_Divisions = 400;

outdir = "Simulation_Results";
if ~exist(outdir, 'dir')
    mkdir(outdir);
end


Filename = fullfile(outdir, "diffusion_vary_3.pdf"); 
firstPage = true;
%% Sweep for all parameter values
for bf = Blood_Fat_Exchange
    for br = Blood_Reaction_Exchange
        for bl = Blood_Losses
            for rl = Reaction_Losses

                % Derivative Matrix
                Derivative_matrix = [-(bl + bf + br),  bf,br;
                     bf, -bf, 0;
                     br,   0, -(br + rl)];
                % RK Solver Fcn
                [SOL, TimeStep] = Runge_Kutta_3_Loop(Derivative_matrix, Initial_Value, End_Time, Time_Divisions);
                figure;
                plot(TimeStep, SOL(1,:), TimeStep, SOL(2,:), TimeStep, SOL(3,:), TimeStep, SOL(1,:)+SOL(2,:) +SOL(3,:));
                    legend("Bloodstream","Storage","Target", "Total"); xlabel('Time (t, Hours)'); ylabel('Concentration (mg/kg)'); grid on;
                    title(sprintf('Runge-Kutta method approximation of three comparment concentrations with N = %d', Time_Divisions));
                %File management
                subtitle(sprintf('BF=%.3f, BR=%.3f, BL=%.3f, RL=%.3f', bf, br, bl, rl));
                if firstPage
                    exportgraphics(gcf, Filename, "Append", false);
                    firstPage = false;
                else
                    exportgraphics(gcf, Filename, "Append", true);
                end

                close(gcf);
            end
        end
    end
end
