function adj = random_planar_graph(n)
    %adj = spalloc(n, n);
    del = delaunay(randn(n, 2));
    for i_triangle = 1:size(del, 1)
        adj(del(i_triangle, 1),del(i_triangle, 2)) = 1;
        adj(del(i_triangle, 2),del(i_triangle, 3)) = 1;
        adj(del(i_triangle, 3),del(i_triangle, 1)) = 1;
    end
end