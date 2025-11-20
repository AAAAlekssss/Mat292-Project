function [SOL, TimeStep] = Runge_Kutta_3_Loop(Derivative_matrix, Init_Value, Final_Time, Div_num)
    Concentrations=@(t) Derivative_matrix*t; %initial value is only in the blood.
    dt = floor(Final_Time/Div_num);
    time = 0:dt:Final_Time;
    SOL = NaN(3,length(time));
    SOL(:,1) = [Init_Value, 0, 0];
    h = floor(Final_Time/Divisions);
    TimeStep = linspace(0, Final_Time, Divisions+1);
    for n = 1:Divisions
        k1 = Concentrations(SOL(:,n));
        k2 = Concentrations(SOL(:,n) + 0.5*h*k1);
        k3 = Concentrations(SOL(:,n) + 0.5*h*k2);
        k4 = Concentrations(SOL(:,n) + h*k3);
        SOL(:,n+1) = SOL(:,n) + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
    end

    plot(TimeStep, SOL(1,:), TimeStep, SOL(2,:), TimeStep, SOL(3,:));
        legend("Bloodstream","Storage","Target"); xlabel('Time (t, Hours)'); ylabel('Concentration (mg/kg'); grid on;
        title('Runge-Kutta method approximation of three comparment concentrations with N = %d' + Divisions);
end
