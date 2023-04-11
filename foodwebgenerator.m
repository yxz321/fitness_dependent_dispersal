function fw = foodwebgenerator(S,C)
% fw = FoodwebGenerator(S,C)
% A niche model modified from Williams & Martinez, 2000.
% ---
% Input parameters
% S - number of species
% C - network connectance of foodweb
% ---
% Output parameters
% fw - S*S adjacent matrix of predatory relationship

    beta = 1 / (2 * C) - 1;
    niche = zeros(S,1);
    randbeta = zeros(S,1);
    preyCenter = zeros(S,1);
    iRemoved = 1:S;
    nRemoved = S; % #=S are removed means #=0 species are present in the beginning
    fw = zeros(S);
    while(nRemoved > 0) % until a connected and non-duplicate web is generated (see below)
        niche(iRemoved) = rand(nRemoved, 1);
        randbeta(iRemoved) = (1 - (1 - rand(nRemoved, 1)) .^ (1 / beta)); % random x's follow f(x|1,beta) in [0,1]
        range = randbeta .* niche;
        preyCenter(iRemoved) = rand(nRemoved, 1) .* (niche(iRemoved) - range(iRemoved) / 2)...
            + range(iRemoved) / 2;
        fw  = ((preyCenter - range / 2) <= niche') &...
            (niche' <= (preyCenter + range / 2));

        fw = fw - diag(diag(fw)); % cannibalism is not allowed
        % isolated species are also removed.
        [iNoPrey,~] = find(sum(fw,2) == 0); % select species without prey
        [~,iNoConsumer] = find(sum(fw,1) == 0); % select species without predator
        iRemoved = intersect(iNoPrey,iNoConsumer); % select unlinked species (with neither prey nor predator)
        if isempty(iRemoved)
            if isempty(find(all(fw==0,2),1)) % if no plants present in the foodweb
                [~,iRemoved] = min(niche); % remove the species with smallest niche value
            end
        end
        nRemoved = numel(iRemoved);
    end
    [~,isort] = sort(niche);
    fw = fw(isort,isort); % sort the species by niche value
    fw = double(fw);
end
