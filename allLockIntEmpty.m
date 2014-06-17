function out = allLockIntEmpty(locks, path)
    out = 1;
    for ti = locks'
        out = out & ((~isempty(locks) && isempty(find(path( ti ).intEmpty == 0))));
    end
end