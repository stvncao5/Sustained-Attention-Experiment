function bootstrappedSamples = bootstrapper(numberOfBootstrapSamples,vectorToBeBootstrapped)
    % Takes a vector of subject IDs (but can also be the data vector itself), and
    % bootstraps it n times in order to create n bootstrapped samples.
    % All bootstrapped samples are of the same size as the vector being bootstrapped.
    
    numElements = length( vectorToBeBootstrapped );
    
    % Preallocate output vector.
    bootstrappedSamples = nan( numberOfBootstrapSamples , numElements );
    
    % For each bootstrap sample...
    for i=1:numberOfBootstrapSamples
        % ...generate a random sample (with replacement) of our ID (or data) vector.
        bootstrapIndices = randsample(1:numElements, numElements, true)';
            % (It's transposed because randsample() returns a column vector for some reason)
        bootstrappedSamples(i,:) = vectorToBeBootstrapped(bootstrapIndices);
    end
    
    
end

