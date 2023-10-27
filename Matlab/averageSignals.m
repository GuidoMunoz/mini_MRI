function [output, t] = averageSignals(nameFormat)
    numSignals = 50;
    signalSum = zeros(8192, 1); % Assuming signals have 1000 data points
    
    % Loop through and accumulate the signals
    for signalNum = 1:numSignals
        filename = sprintf(nameFormat, signalNum);
        data = readtable(filename);
        signalSum = signalSum + table2array(data(:,2)) - mean(table2array(data(:,2)));
    end
    
    % Calculate the average
    averagedSignal = signalSum / numSignals;
    
    % Extract the filename without extension
    %[~, baseFilename, ~] = fileparts(nameFormat);
    
    % Save the averaged signal to a new CSV file
    %outputFilename = [baseFilename, '_av.csv'];
    %csvwrite(outputFilename, averagedSignal);
    output = averagedSignal;
    t = table2array(data(:,1));
    
%    plot(output)
    %fprintf('Averaged signal saved to %s\n', outputFilename);
end
