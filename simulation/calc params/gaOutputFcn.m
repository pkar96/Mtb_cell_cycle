% gaOutputFcn.m
function [state, options, optchanged] = gaOutputFcn(options, state, flag)
    optchanged = false; % Indicate whether the optimization options have changed
    persistent bestLoss;
    persistent bestParams;

    if isempty(bestLoss)
        bestLoss = Inf;
        bestParams = [];
    end

    % Ensure state.Best is non-empty and has valid values before accessing
    if ~isempty(state.Best) && isvector(state.Best)
        currentBest = state.Best(end);
        if state.Generation > 0 && currentBest < bestLoss
            bestLoss = currentBest;
            % Find the index of the best individual
            [~, bestIdx] = min(state.Score);
            bestParams = state.Population(bestIdx, :);

            % Save parameters and loss to a file
            fileID = fopen('optimization_log.txt', 'a');
            fprintf(fileID, 'Generation %d: Best Loss = %f, Params = %s\n', ...
                state.Generation, currentBest, mat2str(state.Population(bestIdx, :), 5));
            fclose(fileID);
        end
        % Log the current state
%         fprintf('Generation %d: Best Loss = %f\n', state.Generation, currentBest);
    end

    % Save the best loss value and parameters to a file at the end of optimization
    if strcmp(flag, 'done')
        save('best_loss.mat', 'bestLoss', 'bestParams');
        fprintf('Best Loss through all Generations: %f\n', bestLoss);
        if ~isempty(bestParams)
            fprintf('Parameters corresponding to the Best Loss: %s\n', mat2str(bestParams, 5));
        else
            fprintf('No valid parameters found for the best loss.\n');
        end
    end
end
