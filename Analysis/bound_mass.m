function [M_bound, ind_bound] = bound_mass(pos, vel, m, method, L)
% Given cloud of particles return largest gravitationally bound mass.

bigG = 6.67384e-11;

switch method
    case 'kory'
        % In the RF of the particle with lowest potential return the particles
        % with negative total energy
        U = bigG*potential(pos(:,1),pos(:,2),pos(:,3),m);
        [~, ind] = min(U);
        VCM = vel(ind,:);
        ind_bound = false(length(pos),1);
        for j=1:length(pos)
            V = vel(j,:) - VCM;
            K = 0.5*(V*V');
            if K + U(j) < 0, ind_bound(j) = true; end
        end
        
    case 'jutzi'
        ind_bound = true(length(pos),1);
        nbb = -1;
        while nbb ~= sum(ind_bound)
            nbb = sum(ind_bound);
            bU = bigG*potential(pos(:,1),pos(:,2),pos(:,3),m,ind_bound);
            [~, ind] = min(bU);
            VCM = vel(ind,:);
            for j=1:length(pos)
                V = vel(j,:) - VCM;
                K = 0.5*(V*V');
                if K + bU(j) >= 0, ind_bound(j) = false; end
            end
        end
        
    case 'naor'
        ind_bound = false(length(pos),1);
        % seed with deepest particle
        U = bigG*potential(pos(:,1),pos(:,2),pos(:,3),m);
        [~, ind] = min(U);
        ind_bound(ind) = true;
        % build bottom-up
        nbb = -1;
        while nbb ~= sum(ind_bound)
            nbb = sum(ind_bound);
            M = sum(m(ind_bound));            
            cmpos = sum(repmat(m(ind_bound),1,3).*pos(ind_bound,:),1)/M;
            cmvel = sum(repmat(m(ind_bound),1,3).*vel(ind_bound,:),1)/M;            
            for j=1:length(pos)
                dr = norm(pos(j,:) - cmpos) + eps;
                U = -bigG*M/dr;
                V = vel(j,:) - cmvel;
                K = 0.5*(V*V');
                if K + U < 0
                    ind_bound(j) = true;
                else
                    ind_bound(j) = false;
                end
            end
        end
        
    otherwise
        error('unknown method')
end
M_bound = sum(m(ind_bound));


function U = potential(x,y,z,m,mask)
if ~exist('mask','var'), mask = true(size(x)); end
U = zeros(size(x));
for j=1:length(x)
    if mask(j)
        for k=1:j-1
            if mask(k)
                dx = x(j) - x(k);
                dy = y(j) - y(k);
                dz = z(j) - z(k);
                dr = sqrt(dx*dx + dy*dy + dz*dz);
                U(j) = U(j) - m(k)/dr;
                U(k) = U(k) - m(j)/dr;
            end
        end
    end
end


function labels = fast_clumps(pos, L)
% Partition a cloud of point masses into distinct clumps based on proximity.

nb_parts = length(pos);
L2 = L^2;
labels = NaN(nb_parts,1);

% The NR algorithm (with the book's exact, and useless, comments)
labels(1)=1;
for j=2:length(pos) % Loop over first element of all pairs.
    labels(j)=j;
    pj = pos(j,:);
    for k=1:j-1      % Loop over second element of all pairs.
        labels(k)=labels(labels(k)); % Sweep it up this much.
        pk = pos(k,:);
        pjk = pj - pk;
        if  pjk(1)^2 + pjk(2)^2 + pjk(3)^2 < L2
            labels(labels(labels(k)))=j;
        end
    end
end
for j=1:length(pos)
    labels(j)=labels(labels(j)); % Only this much sweeping is needed finally.
end
