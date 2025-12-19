%Parameters
Digestion_Blood_Exchange = 0.475;
Digestion_Losses = 0.6;
Blood_Fat_Exchange = 0.195;
Blood_Reaction_Exchange = 0.085;
Blood_Losses = 0.0875;
Reaction_Losses = 0.07;
%parameters tested by seeng what made sense based on literature values - If you need to change these Please use the Param Vary version, then input here
%initial Value Problem
End_Time = 32;
Time_Divisions = 640;
outdir = "Simulation_Results";
if ~exist(outdir, 'dir')
    mkdir(outdir);
end


Filename = fullfile(outdir, "Mass_.pdf"); 
firstPage = true;

Mass Inputs (EDIT HERE!)
Mass_body = [56,66,76]; %in kg - variable dependent on patient mass (other factors are variable too)
%note that it is in a different format to the rest of the mass variables as 
%to ensure it is well differentiated.
%can set ranges as needed
fat_percent = [14,15];
digestive_percent = [8,10];
target_percent = [1];
blood_percent = [8];
Initial_Amount = [10,15,20]; %mg
%do not change anything after this


Create Graph and run iterations
for mb = Mass_body
    for fp = fat_percent
        for dp = digestive_percent
            for tp = target_percent
                for bp = blood_percent
                    for ia = Initial_Amount
                        Digestive_Mass = dp/100*mb;
                        Target_Mass = tp/100*mb;
                        Fat_Mass = fp/100*mb;
                        Blood_Mass = bp/100*mb;
                        %remainder of body - muscles, bones, etc - is assumed to have no effect nor storage capacity.
                        %convert to initial concentration as to begin reaction
                        Initial_Concentration = ia/Digestive_Mass; %mg/kg
                        % Derivative Matrix
                        Derivative_matrix = [-(Digestion_Blood_Exchange+Digestion_Losses), 0, 0,0;
                            Digestion_Blood_Exchange, -(Blood_Losses + Blood_Fat_Exchange + Blood_Reaction_Exchange),  Blood_Fat_Exchange,Blood_Reaction_Exchange;
                             0, Blood_Fat_Exchange, -Blood_Fat_Exchange, 0;
                             0, Blood_Reaction_Exchange,   0, -(Blood_Reaction_Exchange + Reaction_Losses)];
                        % RK Solver Fcn
                        [SOL, TimeStep] = Runge_Kutta_4_Loop(Derivative_matrix, Initial_Concentration, End_Time, Time_Divisions);
                        figure;
                        %return to mass from concentrations
                        SOL(1,:) = Digestive_Mass .* SOL(1,:);
                        SOL(2,:) = Blood_Mass .* SOL(2,:);
                        SOL(3,:) = Fat_Mass .* SOL(3,:);
                        SOL(4,:) = Target_Mass .* SOL(4,:);
                        %plotting graphs
                        plot(TimeStep, SOL(1,:), TimeStep, SOL(2,:), TimeStep, SOL(3,:), TimeStep, SOL(4,:), TimeStep, (SOL(2,:)+SOL(3,:)+SOL(4,:)));
                            legend("Digestion","Bloodstream","Storage","Target", "Body Total"); 
                            xlabel('Time (t, Hours)'); ylabel('mg of medication in each body system'); grid on;
                            title(sprintf('Runge-Kutta approximation of four comparment mg of medication with N = %d', Time_Divisions));
                        %File management shenanigans to make one neat graph
                        subtitle(sprintf('BM=%.3f, FP=%.3f, DP=%.3f, TP=%.3f, BP=%.4f, IA=%.3f', bm, fp, dp, tp, bp, ia));
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
