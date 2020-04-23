classdef PIDController < handle
    properties (GetAccess='public', SetAccess='private')
        kp;
        ki;
        kd;
        t;
        integral;
        prevError;
    end
    
    methods
        % CONSTRUCTOR
        function obj = PIDController(kp, ki, kd)
            obj.kp = kp;
            obj.ki = ki;
            obj.kd = kd;
            obj.t = [];
            obj.integral = 0;
            obj.prevError = 0;
        end
        
        % GET THE CONTROL OUTPUT BASED ON THE CURRENT ERROR AND TIME
        function u = control(obj, t, e)
            if isempty(obj.t)
                obj.t = t;
                obj.prevError = e;
                u = 0;
            else
                dt = t - obj.t;
                obj.t = t;
                obj.integral = obj.integral + (e + obj.prevError) / 2 * dt;
                de = (e - obj.prevError) / dt;
                if isnan(de) || de == inf
                    de = 0;
                end
                u = obj.kp * e + obj.ki * obj.integral + obj.kd * de;
                obj.prevError = e;
            end
        end
        
        % SET ONE OR ALL GAINS
        function setGains(obj, kp, ki, kd)
            if ischar(kp)
                switch kp
                    case 'p'
                        obj.kp = ki;
                    case 'i'
                        obj.ki = ki;
                    case 'd'
                        obj.kd = ki;
                end
                return;
            end
            obj.kp = kp;
            obj.ki = ki;
            obj.kd = kd;
            obj.integral = 0;
            obj.prevError = 0;
        end
        
        % RESET THE CONTROLLER
        function reset(obj)
            obj.t = [];
            obj.integral = 0;
            obj.prevError = 0;
        end
    end
end