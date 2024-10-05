function adj = randomConnectedGraph(n, E)
    adj = spalloc(n, n, E);
    idx = randperm(n * n, E);
    adj(idx) = 1;
    adj = min( adj + adj.', 1);
end