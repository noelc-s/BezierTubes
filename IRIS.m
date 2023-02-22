classdef IRIS
    methods (Static)
        
        function [A,b] = separatingHyperplanes(C, d, obstacle_pts)
            
            persistent mosek_res
            dim = size(C,1);
            infeas_start = false;
            n_obs = size(obstacle_pts, 3);
            pts_per_obs = size(obstacle_pts, 2);
            Cinv = inv(C);
            Cinv2 = (Cinv * Cinv');
            if n_obs == 0 || isempty(obstacle_pts)
                A = zeros(0, dim);
                b = zeros(0, 1);
                infeas_start = false;
                return;
            end
            
            
            uncovered_obstacles = true(n_obs,1);
            planes_to_use = false(n_obs, 1);
            image_pts = reshape(Cinv * bsxfun(@minus, reshape(obstacle_pts, dim, []), d), size(obstacle_pts));
            image_dists = reshape(sum(image_pts.^2, 1), size(obstacle_pts, 2), size(obstacle_pts, 3));
            obs_image_dists = min(image_dists, [], 1);
            [~, obs_sort_idx] = sort(obs_image_dists);
            
            flat_obs_pts = reshape(obstacle_pts, dim, []);
            
            A = zeros(n_obs,dim);
            b = zeros(n_obs,1);
            for i = obs_sort_idx;
                if uncovered_obstacles(i)
                    obs = obstacle_pts(:,:,i);
                    ys = image_pts(:,:,i);
                    dists = image_dists(:,i);
                    [~,idx] = min(dists);
                    xi = obs(:,idx);
                    nhat = 2 * Cinv2 * (xi - d);
                    nhat = nhat / norm(nhat);
                    b0 = nhat' * xi;
                    if all(nhat' * obs - b0 >= 0)
                        % nhat is feasible, so we can skip the optimization
                        A(i,:) = nhat';
                        b(i) = b0;
                    else
                        if isempty(mosek_res)
                            [~,mosek_res] = mosekopt('symbcon echo(0)');
                        end
                        ystar = IRIS.mosek_ldp(ys, mosek_res);
                        
                        if norm(ystar) < 1e-3
                            % d is inside the obstacle. So we'll just reverse nhat to try to push the
                            % ellipsoid out of the obstacle.
                            % warning('IRIS:EllipseCenterInObstacle', 'ellipse center is inside an obstacle.');
                            infeas_start = true;
                            A(i,:) = -nhat';
                            b(i) = -nhat' * xi;
                        else
                            xstar = C*ystar + d;
                            nhat = 2 * Cinv2 * (xstar - d);
                            nhat = nhat / norm(nhat);
                            A(i,:) = nhat;
                            b(i) = nhat' * xstar;
                        end
                    end
                    
                    check = A(i,:) * flat_obs_pts >= b(i);
                    check = reshape(check', pts_per_obs, []);
                    excluded = all(check, 1);
                    uncovered_obstacles(excluded) = false;
                    
                    planes_to_use(i) = true;
                    uncovered_obstacles(i) = false;
                    
                    if ~any(uncovered_obstacles)
                        break
                    end
                end
            end
            A = A(planes_to_use,:);
            b = b(planes_to_use);
        end
        
        function ystar = mosek_ldp(ys, res)
            % Use Mosek to find the closest point in the convex hull of the ys to the
            % origin.
            
            if nargin < 2
                [~, res] = mosekopt('symbcon echo(0)');
            end
            
            dim = size(ys, 1);
            
            nw = size(ys,2);
            nvar= 1 + dim+nw;
            prob.c   = [zeros(1, dim+nw) , 1];
            prob.a   = sparse([ [-eye(dim), ys, zeros(dim,1)];[ zeros(1,dim), ones(1,nw),0] ]);
            prob.blc = [zeros(dim,1);1];
            prob.buc = [zeros(dim,1);1];
            prob.blx = [-inf*ones(dim,1);zeros(nw+1,1)];
            prob.bux = inf*ones(nvar,1);
            
            % Specify the cones.
            prob.cones.type   = res.symbcon.MSK_CT_QUAD;
            prob.cones.sub    = [nvar, 1:dim];
            prob.cones.subptr = 1;
            
            % Optimize the problem.
            [~,solution]=mosekopt('minimize echo(0)',prob);
            %toc
            ystar = solution.sol.itr.xx(1:dim);
            
        end
        
        function [A,b] = separatingHyperplanes_cell(C, d, obstacle_pts)
            
            persistent mosek_res
            dim = size(C,1);
            infeas_start = false;
            n_obs = size(obstacle_pts,2);
            pts_per_obs = cellfun(@(s) size(s,2), obstacle_pts); % now an array
            Cinv = inv(C);
            Cinv2 = (Cinv * Cinv');
            if n_obs == 0 || isempty(obstacle_pts)
                A = zeros(0, dim);
                b = zeros(0, 1);
                infeas_start = false;
                return;
            end
            
            
            uncovered_obstacles = true(n_obs,1);
            planes_to_use = false(n_obs, 1);
            
            image_pts = cellfun(@(s) Cinv*(s - d), obstacle_pts, 'UniformOutput', false);
            image_dists = cellfun(@(s) sum(s.^2), image_pts, 'UniformOutput', false);
            obs_image_dists = cellfun(@(s) min(s), image_dists);
            [~, obs_sort_idx] = sort(obs_image_dists);
            
            flat_obs_pts = [obstacle_pts{:}];
            
            A = zeros(n_obs,dim);
            b = zeros(n_obs,1);
            for i = obs_sort_idx;
                if uncovered_obstacles(i)
                    obs = obstacle_pts{i};
                    ys = image_pts{i};
                    dists = image_dists{i};
                    [~,idx] = min(dists);
                    xi = obs(:,idx);
                    nhat = 2 * Cinv2 * (xi - d);
                    nhat = nhat / norm(nhat);
                    b0 = nhat' * xi;
                    if all(nhat' * obs - b0 >= 0)
                        % nhat is feasible, so we can skip the optimization
                        A(i,:) = nhat';
                        b(i) = b0;
                    else
                        if isempty(mosek_res)
                            [~,mosek_res] = mosekopt('symbcon echo(0)');
                        end
                        ystar = IRIS.mosek_ldp(ys, mosek_res);
                        
                        if norm(ystar) < 1e-3
                            % d is inside the obstacle. So we'll just reverse nhat to try to push the
                            % ellipsoid out of the obstacle.
                            % warning('IRIS:EllipseCenterInObstacle', 'ellipse center is inside an obstacle.');
                            infeas_start = true;
                            A(i,:) = -nhat';
                            b(i) = -nhat' * xi;
                        else
                            xstar = C*ystar + d;
                            nhat = 2 * Cinv2 * (xstar - d);
                            nhat = nhat / norm(nhat);
                            A(i,:) = nhat;
                            b(i) = nhat' * xstar;
                        end
                    end
                    
                    for j = 1:size(obstacle_pts,2)
                        check = A(i,:) * obstacle_pts{j} >= b(i);
                        %             check = reshape(check', pts_per_obs(j), []);
                        uncovered_obstacles(j) = ~all(check);
                    end
                    
                    planes_to_use(i) = true;
                    uncovered_obstacles(i) = false;
                    
                    if ~any(uncovered_obstacles)
                        break
                    end
                end
            end
            A = A(planes_to_use,:);
            b = b(planes_to_use);
        end
        
        function [C, d, volume] = inscribedEllipsoid(A,b)
            
            % poly = iris.Polyhedron(A, b).reduce();
            [Ad, ia] = unique(A,'rows');
            A = Ad;
            b = b(ia);
            
            [C, d] = IRIS.mosek_nofusion(A, b);
            % [C, d] = iris.inner_ellipsoid.mosek_ellipsoid(A, b);
            
            % If Mosek fails for you, you can use CVX with the free SDPT3 solver,
            % but it will be much (about 100X) slower. Just swap the above line for the
            % following:
            % [C, d] = iris.inner_ellipsoid.cvx_ellipsoid(A, b);
            
            volume = det(C);
            
        end
        
        function [C, d] = mosek_nofusion(A, b)
            % Find the largest ellipsoid in the polyhedron defined by Ax <= b. This
            % should return the same result as iris.inner_ellipsoid.mosek_ellipsoid.m.
            % It implements the same algorithm, but does not use the Mosek Fusion
            % symbolic API, and is consequently about 2X faster. You can compare all of
            % the available ellipsoid solvers with iris.test.test_ellipsoid.m.
            
            
            DEBUG = false;
            
            % tic
            [m, n] = size(A);
            l = ceil(log2(n));
            
            num.t = 1;
            num.d = n;
            num.s = 2^l - 1;
            num.sprime = num.s;
            num.z = 2^l;
            num.f = m * n;
            num.g = m;
            
            nvar = 0;
            for v = {'t', 'd', 's', 'sprime', 'z', 'f', 'g'}
                var = v{1};
                ndx.(var) = nvar + (1:num.(var));
                nvar = nvar + num.(var);
            end
            ndx.f = reshape(ndx.f, m, n);
            
            ncon = n * m + m + n + n + (2^l - n) + 1 + (n * (n-1) / 2) + (2^l - 1);
            
            nabar = n * m * n + n + n + (n * (n-1) / 2);
            abar_ptr = 1;
            
            % Mosek suggests running the following so we can look up the cone types:
            % [r, res] = mosekopt('symbcon');
            % but this takes about as long as actually solving the SDP, so we'll just
            % cowboy up and hard-code them.
            MSK_CT_RQUAD = 1;
            MSK_CT_QUAD = 0;
            
            prob.c = zeros(nvar, 1);
            
            % maximize t
            prob.c(ndx.t) = 1;
            
            % Y \in S^{2n}_+
            prob.bardim = 2*n;
            
            prob.a = zeros(ncon, nvar);
            prob.blc = -inf(ncon,1);
            prob.buc = inf(ncon,1);
            
            prob.bara.subi = nan(1, nabar);
            prob.bara.subj = nan(1, nabar);
            prob.bara.subk = nan(1, nabar);
            prob.bara.subl = nan(1, nabar);
            prob.bara.val = nan(1, nabar);
            
            prob.cones.type = [];
            prob.cones.sub = [];
            prob.cones.subptr = [];
            
            con_ndx = 1;
            for i = 1:m
                % a_i^T C = [f_{i,1}, f_{i,2}, ..., f_{i,n}]
                for j = 1:n
                    % (a_i^T C)_j = f_{i,j}
                    prob.bara.subi(abar_ptr:abar_ptr+n-1) = con_ndx;
                    prob.bara.subj(abar_ptr:abar_ptr+n-1) = 1;
                    
                    % Do some silliness because Mosek will fail if we try to specify
                    % elements of Abar above the diagonal.
                    subk = j + zeros(1, n);
                    subl = 1:n;
                    swap_mask = subk < subl;
                    swap = subk(swap_mask);
                    subk(swap_mask) = subl(swap_mask);
                    subl(swap_mask) = swap;
                    
                    prob.bara.subk(abar_ptr:abar_ptr+n-1) = subk;
                    prob.bara.subl(abar_ptr:abar_ptr+n-1) = subl;
                    prob.bara.val(abar_ptr:abar_ptr+n-1) = A(i,:);
                    abar_ptr = abar_ptr + n;
                    
                    prob.a(con_ndx, ndx.f(i,j)) = -1;
                    prob.blc(con_ndx) = 0;
                    prob.buc(con_ndx) = 0;
                    con_ndx = con_ndx + 1;
                end
                prob.a(con_ndx, ndx.d) = A(i,:);
                prob.a(con_ndx, ndx.g(i)) = 1;
                prob.blc(con_ndx) = b(i);
                prob.buc(con_ndx) = b(i);
                con_ndx = con_ndx + 1;
            end
            
            for j = 1:n
                % Xbar_{n+j,n+j} == z_j
                prob.bara.subi(abar_ptr) = con_ndx;
                prob.bara.subj(abar_ptr) = 1;
                prob.bara.subk(abar_ptr) = n+j;
                prob.bara.subl(abar_ptr) = j;
                prob.bara.val(abar_ptr) = 1;
                abar_ptr = abar_ptr + 1;
                
                prob.a(con_ndx, ndx.z(j)) = -1;
                prob.blc(con_ndx) = 0;
                prob.buc(con_ndx) = 0;
                con_ndx = con_ndx + 1;
            end
            
            for j = 1:n
                % Xbar_{n+j,n+j} == z_j
                prob.bara.subi(abar_ptr) = con_ndx;
                prob.bara.subj(abar_ptr) = 1;
                prob.bara.subk(abar_ptr) = n+j;
                prob.bara.subl(abar_ptr) = n+j;
                prob.bara.val(abar_ptr) = 1;
                abar_ptr = abar_ptr + 1;
                
                prob.a(con_ndx, ndx.z(j)) = -1;
                prob.blc(con_ndx) = 0;
                prob.buc(con_ndx) = 0;
                con_ndx = con_ndx + 1;
            end
            
            for j = n+1:2^l
                % z_j == t for j > n
                prob.a(con_ndx, ndx.z(j)) = 1;
                prob.a(con_ndx, ndx.t) = -1;
                prob.blc(con_ndx) = 0;
                prob.buc(con_ndx) = 0;
                con_ndx = con_ndx + 1;
            end
            
            % Off-diagonal elements of Y22 are 0
            for k = n+1:(2*n-1)
                for j = k+1:2*n
                    prob.bara.subi(abar_ptr) = con_ndx;
                    prob.bara.subj(abar_ptr) = 1;
                    prob.bara.subk(abar_ptr) = j;
                    prob.bara.subl(abar_ptr) = k;
                    prob.bara.val(abar_ptr) = 1;
                    abar_ptr = abar_ptr + 1;
                    
                    prob.blc(con_ndx) = 0;
                    prob.buc(con_ndx) = 0;
                    con_ndx = con_ndx + 1;
                end
            end
            
            if DEBUG
                assert(~any(isnan(prob.bara.subi)));
                assert(~any(isnan(prob.bara.subj)));
                assert(~any(isnan(prob.bara.subk)));
                assert(~any(isnan(prob.bara.subl)));
                assert(~any(isnan(prob.bara.val)));
            end
            
            % 2^(l/2)t == s_{2l - 1}
            prob.a(con_ndx, ndx.t) = 2^(l/2);
            prob.a(con_ndx, ndx.s(end)) = -1;
            prob.blc(con_ndx) = 0;
            prob.buc(con_ndx) = 0;
            con_ndx = con_ndx + 1;
            
            
            % s_j == sprime_j
            for j = 1:(2^l - 1)
                prob.a(con_ndx, ndx.s(j)) = 1;
                prob.a(con_ndx, ndx.sprime(j)) = -1;
                prob.blc(con_ndx) = 0;
                prob.buc(con_ndx) = 0;
                con_ndx = con_ndx + 1;
            end
            
            if DEBUG
                assert(con_ndx == ncon + 1);
            end
            
            cone_ptr = 1;
            lhs = [ndx.z, ndx.sprime];
            lhs_ptr = 1;
            for j = 1:(2^l - 1)
                prob.cones.type = [prob.cones.type, MSK_CT_RQUAD];
                prob.cones.sub = [prob.cones.sub, lhs(lhs_ptr:lhs_ptr+1), ndx.s(j)];
                prob.cones.subptr = [prob.cones.subptr, cone_ptr];
                lhs_ptr = lhs_ptr + 2;
                cone_ptr = cone_ptr + 3;
            end
            
            for i = 1:m
                prob.cones.type = [prob.cones.type, MSK_CT_QUAD];
                prob.cones.sub = [prob.cones.sub, ndx.g(i), ndx.f(i,:)];
                prob.cones.subptr = [prob.cones.subptr, cone_ptr];
                cone_ptr = cone_ptr + n + 1;
            end
            
            prob.a = sparse(prob.a);
            
            % Divide all off-diagonal entries of Abar by 2. This is necessary because Abar
            % is assumed by the solver to be a symmetric matrix, but we're only setting
            % its lower triangular part.
            mask = prob.bara.subk ~= prob.bara.subl;
            prob.bara.val(mask) = prob.bara.val(mask) / 2;
            
            % fprintf('setup: %f s\n', toc);
            % tic
            [r, res] = mosekopt('maximize echo(0)', prob);
            if ~isfield(res, 'sol')
                res
                error('IRIS:MosekNoSolution', sprintf('MOSEK was unable to return a solution. %s', res.rmsg));
            end
            % fprintf('solve: %f s\n', toc);
            
            % tic
            Y = zeros(2*n, 2*n);
            flat_ndx = 1;
            for k = 1:2*n
                for j = k:2*n
                    Y(j,k) = res.sol.itr.barx(flat_ndx);
                    flat_ndx = flat_ndx + 1;
                end
            end
            
            % Reflect Y to turn the lower-triangular form into a full matrix
            Y = Y + tril(Y,-1)';
            C = Y(1:n, 1:n);
            d = res.sol.itr.xx(ndx.d);
            % fprintf('extract: %f s\n', toc)
            
        end
    end
end