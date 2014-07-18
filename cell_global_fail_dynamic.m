% fail a cell after we've started running (for plots)
if k == failround && opt_faildynamic
    failed = [failed, faildynamic];

    % set distance and next pointers for failures
    if ~isempty(failed)
        for i = 1 : length(failed)
            f = failed(i);
            Cell(f).failed = 1;
            for tt = 1 : NT
                Cell(f).prev(tt) = f;
                Cell(f).next(tt) = f;
                Cell(f).dist(tt) = inf;
            end
        end
    end
end