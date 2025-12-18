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
%%
%[text] ## Mass Inputs (EDIT HERE!)
Mass_body = 76; %in kg - variable dependent on patient mass (other factors are variable too)
%note that it is in a different format to the rest of the mass variables as 
%to ensure it is well differentiated.
fat_percent = 14;
digestive_percent = 8;
target_percent = 1;
blood_percent = 8;
Initial_Amount = 10; %mg
%do not change anything after this
Digestive_Mass = digestive_percent/100*Mass_body;
Target_Mass = target_percent/100*Mass_body;
Fat_Mass = fat_percent/100*Mass_body;
Blood_Mass = blood_percent/100*Mass_body;
%remainder of body - muscles, bones, etc - is assumed to have no effect nor storage capacity.
%convert to initial concentration as to begin reaction
Initial_Concentration = Initial_Amount/Digestive_Mass; %mg/kg

%%
%[text] ## Create Graph
% Derivative Matrix
Derivative_matrix = [-(db+dl), 0, 0,0;
    db, -(bl + bf + br),  bf,br;
     0, bf, -bf, 0;
     0, br,   0, -(br + rl)];
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
    title(sprintf('Runge-Kutta method approximation of four comparment concentrations with N = %d', Time_Divisions));
%File management shenanigans to make one neat graph
subtitle(sprintf('DB=%.3f, DL=%.3f, BF=%.3f, BR=%.3f, BL=%.4f, RL=%.3f', db, dl, bf, br, bl, rl));
if firstPage
    exportgraphics(gcf, Filename, "Append", false);
    firstPage = false;
else
    exportgraphics(gcf, Filename, "Append", true);
end

close(gcf);

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":21.9}
%---
