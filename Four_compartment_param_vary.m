%Parameters - CHange Me!
Digestion_Blood_Exchange = 0.475;
Digestion_Losses = 0.6;
Blood_Fat_Exchange = 0.195;
Blood_Reaction_Exchange = 0.085;
Blood_Losses = 0.0875;
Reaction_Losses = 0.07;
%initial Value Problem
Initial_Value = 0.5;
End_Time = 20; %Time of simulation - change as needed!
Time_Divisions = 400; %I found that 20* whatever you had as End Time was best

outdir = "Simulation_Results";
if ~exist(outdir, 'dir')
    mkdir(outdir);
end


Filename = fullfile(outdir, "Reactions_Testing.pdf"); 
firstPage = true;
%%
%[text] ## Sweep for all parameter values
for db = Digestion_Blood_Exchange
    for dl = Digestion_Losses
        for bf = Blood_Fat_Exchange
            for br = Blood_Reaction_Exchange
                for bl = Blood_Losses
                    for rl = Reaction_Losses
        
                        % Derivative Matrix
                        Derivative_matrix = [-(db+dl), 0, 0,0;
                            db, -(bl + bf + br),  bf,br;
                             0, bf, -bf, 0;
                             0, br,   0, -(br + rl)];
                        % RK Solver Fcn
                        [SOL, TimeStep] = Runge_Kutta_4_Loop(Derivative_matrix, Initial_Value, End_Time, Time_Divisions);
                        figure;
                        plot(TimeStep, SOL(1,:), TimeStep, SOL(2,:), TimeStep, SOL(3,:), TimeStep, SOL(4,:), TimeStep, SOL(2,:)+SOL(3,:) +SOL(4,:));
                    legend("Digestion","Bloodstream","Storage","Target", "Body Total"); xlabel('Time (t, Hours)'); ylabel('Concentration (mg/kg)'); grid on;
                            title(sprintf('Runge-Kutta method approximation of four comparment concentrations with N = %d', Time_Divisions));
                        %File management
                        subtitle(sprintf('DB=%.3f, DL=%.3f, BF=%.3f, BR=%.3f, BL=%.4f, RL=%.3f', db, dl, bf, br, bl, rl));
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
    end
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":21.9}
%---
