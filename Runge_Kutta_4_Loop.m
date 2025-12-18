function [SOL, TimeStep] = Runge_Kutta_4_Loop(Derivative_matrix, Init_Value, Final_Time, Div_num)
    Concentrations=@(t) Derivative_matrix*t; %initial value is only in the blood.
    dt = Final_Time/Div_num;
    time = 0:dt:Final_Time;
    SOL = NaN(4,length(time));
    SOL(:,1) = [Init_Value, 0, 0, 0];
    TimeStep = linspace(0, Final_Time, Div_num+1);
    for n = 1:Div_num
        k1 = Concentrations(SOL(:,n));
        k2 = Concentrations(SOL(:,n) + 0.5*dt*k1);
        k3 = Concentrations(SOL(:,n) + 0.5*dt*k2);
        k4 = Concentrations(SOL(:,n) + dt*k3);
        SOL(:,n+1) = SOL(:,n) + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    end
end