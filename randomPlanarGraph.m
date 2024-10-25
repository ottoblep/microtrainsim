function adj = randomPlanarGraph(n)
    % adj = spalloc(n, n, 6*n*n); % for finite planar graphs the average degree is strictly less than 6
    del = delaunay(randn(n, 2));
    for i_triangle = 1:size(del, 1)
        adj(del(i_triangle, 1),del(i_triangle, 2)) = 1;
        adj(del(i_triangle, 2),del(i_triangle, 3)) = 1;
        adj(del(i_triangle, 3),del(i_triangle, 1)) = 1;
    end
end