%function [outputArg1,outputArg2] = untitled(inputArg1,inputArg2)
function stop = outputFunction_2(optimValues, state)
persistent foundLocal
stop = false;
switch state
    case 'init'
        foundLocal = []; % initialized as empty
    case 'iter'
        newf = optimValues.localsolution.Fval;
        eflag = optimValues.localsolution.Exitflag;
        %temp_param = [newf, optimValues.bestx];
        temp_param = [newf, optimValues.localsolution.X];
        % Now check if the exit flag is positive and
        % the new value differs from all others by at least 1e-4
        % If so, add the new value to the newf list
        if  newf <= optimValues.bestfval %eflag >= 0 &&
            foundLocal = [foundLocal;newf];
            %dlmwrite('test.dat',newf,'delimiter', '\t','-append');
            dlmwrite('Results\Glob_sol_All.dat',temp_param,'delimiter', '\t','newline', 'pc','-append');
            %keyboard
            % Now check if the latest value added to foundLocal
            % is less than 1/2
            % Also check if there are 5 local minima in foundLocal
            % If so, then stop
            if foundLocal(end) < 0.005 %|| length(foundLocal) >= 5
                stop = true;
            end
        end
        %keyboard
%outputArg1 = inputArg1;
%outputArg2 = inputArg2;
end

