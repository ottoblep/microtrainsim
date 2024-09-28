% Copyright (c) 2009, Dustin Arendt
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%
% FastFloyd - quickly compute the all pairs shortest path matrix
% 
% Uses a vectorized version of the Flyod-Warshall algorithm
% see: http://en.wikipedia.org/wiki/Floyd_Warshall
% 
% USAGE:
%
% D = FastFloyd(A)
% 
% D - the distance (geodesics, all pairs shortest path, etc.) matrix.
% A - the adjacency matrix, where A(i,j) is the cost for moving from vertex i to
%     vertex j.  If vertex i and vertex j are not connected then A(i,j) should
%     be >= the diameter of the network (Inf works fine).
%      
% EXAMPLE:
% 
% Here I create a random binary matrix and convert it to an integer format. Then
% I take the reciprocal of the matrix so that all non-adjacent pairs get a value
% of Inf.  The result is stored in D.
%
% A = int32(rand(100,100) < 0.05);
% D = FastFloyd(1./A)
%

function D = FastFloyd(D)

	n = size(D, 1);	

	for k=1:n
	
		i2k = repmat(D(:,k), 1, n);
		k2j = repmat(D(k,:), n, 1);
		
		D = min(D, i2k+k2j);
	
	end

end