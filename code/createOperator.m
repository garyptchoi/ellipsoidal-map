function Operator = createOperator(v,f)
% -------------------------------------------------------
% Create common discrete operators
% -------------------------------------------------------
% Input : 
% 	v - vertex (3xn)
% 	f - connection (3xm)
%
% Output :
% 	f2v - Face to vertex operator
%		  (Face-valued to vertex-valued function)
%	v2f - Vertex to face operator
%		  (Vertex-valued to face-valued function)
% -------------------------------------------------------
%
% If you use this code in your own work, please cite the following paper:
% [1] P. T. Choi, K. C. Lam, and L. M. Lui, 
%     "FLASH: Fast Landmark Aligned Spherical Harmonic Parameterization for Genus-0 Closed Brain Surfaces."
%     SIAM Journal on Imaging Sciences, vol. 8, no. 1, pp. 67-94, 2015.
%
% Copyright (c) 2013-2022, Gary Pui-Tung Choi
% https://math.mit.edu/~ptchoi/

Operator.f2v = F2V(v,f);
Operator.v2f = V2F(v,f);

end




% -------------------------------------------------------
function S = F2V(v,f)
    ring = vertexAttachments(TriRep(f',v'));
    nv = length(v); nf = length(f);
    II = cellfun(@times,ring,num2cell(zeros(nv,1)),'UniformOutput',0);
    II = cell2mat(cellfun(@plus,II,num2cell(1:nv)','UniformOutput',0)')';
    JJ = cell2mat(ring')';
    avg = cellfun(@length,ring);
    S = sparse(II,JJ,ones(length(JJ),1),nv,nf);
    S = sparse(1:nv,1:nv,1./avg)*S;
end


% -------------------------------------------------------
function S = V2F(v,f)
    nv = length(v); nf = length(f);
    II = reshape(repmat(1:nf,3,1),3*nf,1);
    JJ = f(:);
    S = sparse(II,JJ,ones(length(JJ),1),nf,nv)./3;
end

