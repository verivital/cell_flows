% error checking
if sum(targets > N) ~= 0 || sum(sources > N) ~= 0 || sum(failed > N) ~= 0 || sum(faildynamic > N) ~= 0
    'ERROR: bad source, target, or failed IDs entered (greater than N)'
    return
end