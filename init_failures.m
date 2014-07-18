%fail F cells randomly
if opt_fail_random
    for i = 1 : F
        f = ceil(rand(1,1)*N);

        % don't fail targets or sources
        if intersect([targets; sources], f)
            continue;
        end
        failed = [failed; f];
    end
end

% set distance and next pointers for failures
if ~isempty(failed)
    for i = 1 : length(failed)
        f = failed(i);
        Cell(f).failed = 1;
        for tt = 1 : NT
            Cell(f).path(tt).prev = [];
            Cell(f).next(tt) = f;
            Cell(f).dist(tt) = inf;
        end
    end
end